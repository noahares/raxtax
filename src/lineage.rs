use crate::tree::{Node, Tree};
use itertools::Itertools;
use logging_timer::time;

use crate::utils;

#[derive(Debug, Clone)]
pub struct EvaluationResult<'a, 'b> {
    pub query_label: &'b String,
    pub lineage: &'a String,
    pub confidence_values: Vec<f64>,
    pub local_signal: f64,
    pub global_signal: f64,
}

impl EvaluationResult<'_, '_> {
    pub fn get_output_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{:.5}\t{:.5}",
            self.query_label,
            self.lineage,
            self.confidence_values
                .iter()
                .map(|v| format!("{1:.0$}", utils::F64_OUTPUT_ACCURACY as usize, v))
                .join(","),
            self.local_signal,
            self.global_signal
        )
    }

    pub fn get_tsv_string(&self, sequence: &String) -> String {
        format!(
            "{}\t{}\t{:.5}\t{:.5}\t{}",
            self.query_label,
            self.lineage
                .split(',')
                .map(std::string::ToString::to_string)
                .interleave(self.confidence_values.iter().map(|v| format!(
                    "{1:.0$}",
                    utils::F64_OUTPUT_ACCURACY as usize,
                    v
                )))
                .join("\t"),
            self.local_signal,
            self.global_signal,
            sequence
        )
    }
}

pub struct Lineage<'a, 'b> {
    query_label: &'b String,
    tree: &'a Tree,
    confidence_values: Vec<f64>,
    confidence_prefix_sum: Vec<f64>,
    confidence_vectors: Vec<(usize, Vec<f64>, Vec<f64>)>,
    rounding_factor: f64,
}

impl<'a, 'b> Lineage<'a, 'b> {
    pub fn new(query_label: &'b String, tree: &'a Tree, confidence_values: Vec<f64>) -> Self {
        let mut confidence_prefix_sum = vec![0.0];
        confidence_prefix_sum.extend(confidence_values.iter().scan(0.0, |sum, i| {
            *sum += i;
            Some(*sum)
        }));
        let rounding_factor = f64::from(10_u32.pow(utils::F64_OUTPUT_ACCURACY));
        let expected_num_results = rounding_factor as usize / 2;
        Self {
            query_label,
            tree,
            confidence_values,
            confidence_prefix_sum,
            confidence_vectors: Vec::with_capacity(expected_num_results),
            rounding_factor,
        }
    }

    #[time("debug")]
    pub fn evaluate(mut self) -> Vec<EvaluationResult<'a, 'b>> {
        self.eval_recurse(&self.tree.root, &[], &[]);
        // NOTE: This would be the correct maximum leaf confidence and ideally we would normalize with this.
        // However, because this is already 0.99 for 100 tips, it is not worth it, as it is
        // basically 1 for any reasonable reference lineage.
        // let max_leaf_confidence = ((1.0 - 1.0 / self.tree.num_tips as f64).powi(2) + ((self.tree.num_tips as f64 - 1.0) / (self.tree.num_tips as f64).powi(2))).sqrt();
        let leaf_confidence = utils::euclidean_norm(
            self.confidence_values
                .iter()
                .map(|&v| (v - 1.0 / self.tree.num_tips as f64)),
        );
        self.confidence_vectors
            .into_iter()
            .sorted_by(|a, b| b.1.iter().partial_cmp(a.1.iter()).unwrap())
            .map(|(idx, conf_values, expected_conf_values)| {
                let start_index = match expected_conf_values.iter().find_position(|&&x| 1.0 > x) {
                    Some((i, _)) => i,
                    None => expected_conf_values.len() - 1,
                };
                let lineage_confidence = utils::euclidean_distance_l1(
                    &conf_values[start_index..],
                    &expected_conf_values[start_index..],
                );
                EvaluationResult {
                    query_label: self.query_label,
                    lineage: &self.tree.lineages[idx],
                    confidence_values: conf_values,
                    local_signal: lineage_confidence,
                    global_signal: leaf_confidence,
                }
            })
            .collect_vec()
    }

    fn get_confidence(&self, node: &Node) -> f64 {
        self.confidence_prefix_sum[node.confidence_range.1]
            - self.confidence_prefix_sum[node.confidence_range.0]
    }

    fn eval_recurse(
        &mut self,
        node: &Node,
        confidence_prefix: &[f64],
        expected_confidence_prefix: &[f64],
    ) -> bool {
        let mut no_child_significant = true;
        let mut pushed_result = false;
        for c in &node.children {
            let child_conf =
                (self.get_confidence(c) * self.rounding_factor).round() / self.rounding_factor;
            if child_conf == 0.0 {
                continue;
            }
            no_child_significant = false;
            let mut conf_prefix = confidence_prefix.to_vec();
            let mut expected_conf_prefix = expected_confidence_prefix.to_vec();
            conf_prefix.push(child_conf);
            expected_conf_prefix.push(
                (c.confidence_range.1 - c.confidence_range.0) as f64 / self.tree.num_tips as f64,
            );
            let child_pushed_result = self.eval_recurse(c, &conf_prefix, &expected_conf_prefix);
            if !child_pushed_result && self.tree.is_taxon_leaf(c) {
                self.confidence_vectors.push((
                    c.confidence_range.0,
                    conf_prefix,
                    expected_conf_prefix,
                ));
                pushed_result = true;
            }
            pushed_result |= child_pushed_result;
        }
        if no_child_significant && self.tree.is_inner_taxon_node(node) {
            let mut conf_prefix = confidence_prefix.to_vec();
            let mut expected_conf_prefix = expected_confidence_prefix.to_vec();
            let mut current_node = node;
            while self.tree.is_inner_taxon_node(current_node) {
                current_node = current_node
                    .children
                    .iter()
                    .max_by(|c, d| {
                        self.get_confidence(c)
                            .partial_cmp(&self.get_confidence(d))
                            .unwrap()
                    })
                    .unwrap();
                conf_prefix.push(1.0 / self.rounding_factor);
                expected_conf_prefix.push(
                    (current_node.confidence_range.1 - current_node.confidence_range.0) as f64
                        / self.tree.num_tips as f64,
                );
            }
            self.confidence_vectors.push((
                current_node.confidence_range.0,
                conf_prefix,
                expected_conf_prefix,
            ));
            pushed_result = true;
        }
        pushed_result
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::{
        lineage::{EvaluationResult, Lineage},
        tree::Tree,
    };

    #[test]
    fn test_tree_construction() {
        let lineages = vec![
            String::from("Animalia,Chordata,Mammalia,Primates,Hominidae,Homo"),
            "Animalia,Chordata,Mammalia,Primates,Hominidae,Pan".into(),
            "Animalia,Chordata,Mammalia,Carnivora,Canidae,Canis".into(),
            "Animalia,Chordata,Mammalia,Carnivora,Felidae,Felis".into(),
            "Animalia,Chordata,Mammalia,Carnivora,Felidae,Felis".into(),
        ];
        let sequences = vec![
            [0b00].repeat(9),
            [0b00].repeat(9),
            [0b00].repeat(9),
            [0b00].repeat(9),
            [0b00].repeat(9),
        ];
        let tree = Tree::new(lineages, sequences).unwrap();
        let confidence_values = vec![0.1, 0.3, 0.4, 0.004, 0.004];
        tree.print();
        let query_label = String::from("q");
        let lineage = Lineage::new(&query_label, &tree, confidence_values);
        let result = lineage.evaluate();
        assert_eq!(
            result
                .into_iter()
                .map(
                    |EvaluationResult {
                         lineage,
                         confidence_values,
                         ..
                     }| (lineage, confidence_values)
                )
                .collect_vec(),
            vec![
                (
                    &String::from("Animalia,Chordata,Mammalia,Carnivora,Felidae,Felis"),
                    vec![0.81, 0.81, 0.81, 0.8, 0.7, 0.7,],
                ),
                (
                    &"Animalia,Chordata,Mammalia,Carnivora,Canidae,Canis".into(),
                    vec![0.81, 0.81, 0.81, 0.8, 0.1, 0.1,],
                ),
                (
                    &"Animalia,Chordata,Mammalia,Primates,Hominidae,Pan".into(),
                    vec![0.81, 0.81, 0.81, 0.01, 0.01, 0.01,],
                ),
            ]
        );
    }

    #[test]
    fn test_variable_lineage_length() {
        let lineages = vec![
            String::from("Animalia,Chordata,Mammalia,Primates,Hominidae,Homo,Homo_sapiens"),
            "Animalia,Chordata,Mammalia,Primates,Hominidae,Pan".into(),
            "Animalia,Chordata,Mammalia,Carnivora,Canidae,Canis".into(),
            "Animalia,Chordata,Mammalia,Carnivora,Doggo".into(),
            "Animalia,Chordata,Mammalia,Mouse".into(),
            "Animalia,Chordata,Mammalia,Carnivora,Felidae,Felis".into(),
            "Animalia,Chordata,Mammalia,Carnivora,Felidae,Felis".into(),
        ];
        let sequences = vec![
            [0b00].repeat(9),
            [0b00].repeat(9),
            [0b00].repeat(9),
            [0b00].repeat(9),
            [0b00].repeat(9),
            [0b00].repeat(9),
            [0b00].repeat(9),
        ];
        let tree = Tree::new(lineages, sequences).unwrap();
        let confidence_values = vec![0.05, 0.1, 0.3, 0.4, 0.1, 0.004, 0.004];
        tree.print();
        let query_label = String::from("q");
        let lineage = Lineage::new(&query_label, &tree, confidence_values);
        let result = lineage.evaluate();
        dbg!(&result);
        assert_eq!(
            result
                .into_iter()
                .map(
                    |EvaluationResult {
                         lineage,
                         confidence_values,
                         ..
                     }| (lineage, confidence_values)
                )
                .collect_vec(),
            vec![
                (
                    &String::from("Animalia,Chordata,Mammalia,Carnivora,Felidae,Felis"),
                    vec![0.96, 0.96, 0.96, 0.85, 0.7, 0.7,],
                ),
                (
                    &"Animalia,Chordata,Mammalia,Carnivora,Doggo".into(),
                    vec![0.96, 0.96, 0.96, 0.85, 0.1,],
                ),
                (
                    &"Animalia,Chordata,Mammalia,Carnivora,Canidae,Canis".into(),
                    vec![0.96, 0.96, 0.96, 0.85, 0.05, 0.05,],
                ),
                (
                    &"Animalia,Chordata,Mammalia,Mouse".into(),
                    vec![0.96, 0.96, 0.96, 0.1],
                ),
                (
                    &"Animalia,Chordata,Mammalia,Primates,Hominidae,Pan".into(),
                    vec![0.96, 0.96, 0.96, 0.01, 0.01, 0.01,],
                ),
            ]
        );
    }

    #[test]
    fn test_likelihood_edge_case() {
        let lineages = vec![
            "Animalia,Chordata,Mammalia,Carnivora,Felidae,Felis".into(),
            "Animalia,Chordata,Mammalia,Carnivora,Felidae,Felis_ferrocius".into(),
            "Animalia,Chordata,Mammalia,Carnivora,Canidae,Canis".into(),
        ];
        let sequences = vec![[0b00].repeat(9), [0b00].repeat(9), [0b00].repeat(9)];
        let tree = Tree::new(lineages, sequences).unwrap();
        let confidence_values = vec![0.004, 0.004, 0.004];
        tree.print();
        let query_label = String::from("q");
        let lineage = Lineage::new(&query_label, &tree, confidence_values);
        let result = lineage.evaluate();
        assert_eq!(
            result
                .into_iter()
                .map(
                    |EvaluationResult {
                         lineage,
                         confidence_values,
                         ..
                     }| (lineage, confidence_values)
                )
                .collect_vec(),
            vec![(
                &String::from("Animalia,Chordata,Mammalia,Carnivora,Felidae,Felis_ferrocius"),
                vec![0.01, 0.01, 0.01, 0.01, 0.01, 0.01,],
            ),]
        );
    }
}
