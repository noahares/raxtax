use anyhow::Result;
use std::{
    fs::File,
    io::{Read, Write},
    path::PathBuf,
};

use ahash::{HashMap, HashMapExt};
use indicatif::{ProgressIterator, ProgressStyle};
use itertools::Itertools;
use log::{log_enabled, Level};
use logging_timer::time;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::utils::{self, map_four_to_two_bit_repr};

#[derive(Debug)]
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
                .map(|s| s.to_string())
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
    pub fn new(query_label: &'b String, tree: &'a Tree, confidence_values: &[f64]) -> Self {
        let mut confidence_prefix_sum = vec![0.0];
        confidence_prefix_sum.extend(confidence_values.iter().scan(0.0, |sum, i| {
            *sum += i;
            Some(*sum)
        }));
        Self {
            query_label,
            tree,
            confidence_values: confidence_values.to_vec(),
            confidence_prefix_sum,
            confidence_vectors: Vec::new(),
            rounding_factor: 10_u32.pow(utils::F64_OUTPUT_ACCURACY) as f64,
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
                let start_index = expected_conf_values
                    .iter()
                    .find_position(|&&x| 1.0 - x > std::f64::EPSILON)
                    .unwrap()
                    .0;
                let secondary_conf = utils::euclidean_distance_l1(
                    &conf_values[start_index..],
                    &expected_conf_values[start_index..],
                );
                EvaluationResult {
                    query_label: self.query_label,
                    lineage: &self.tree.lineages[idx],
                    confidence_values: conf_values,
                    local_signal: secondary_conf,
                    global_signal: leaf_confidence,
                }
            })
            .collect_vec()
    }

    fn get_confidence(&self, node: &TreeNode) -> f64 {
        self.confidence_prefix_sum[node.confidence_range.1]
            - self.confidence_prefix_sum[node.confidence_range.0]
    }

    fn eval_recurse(
        &mut self,
        node: &TreeNode,
        confidence_prefix: &[f64],
        expected_confidence_prefix: &[f64],
    ) {
        let mut no_child_significant = true;
        for c in &node.children {
            let mut conf_prefix = confidence_prefix.to_vec();
            let mut expected_conf_prefix = expected_confidence_prefix.to_vec();
            let child_conf =
                (self.get_confidence(c) * self.rounding_factor).round() / self.rounding_factor;
            if child_conf == 0.0 {
                continue;
            }
            no_child_significant = false;
            conf_prefix.push(child_conf);
            expected_conf_prefix.push(
                (c.confidence_range.1 - c.confidence_range.0) as f64 / self.tree.num_tips as f64,
            );
            self.eval_recurse(c, &conf_prefix, &expected_conf_prefix);
            if self.tree.is_taxon_leaf(c) {
                self.confidence_vectors.push((
                    c.confidence_range.0,
                    conf_prefix,
                    expected_conf_prefix,
                ));
            }
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
        }
    }
}

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct Tree {
    root: TreeNode,
    pub lineages: Vec<String>,
    pub sequences: HashMap<Vec<u8>, Vec<usize>>,
    pub k_mer_map: Vec<Vec<usize>>,
    pub num_tips: usize,
}

impl Tree {
    #[time("info", "Tree::{}")]
    pub fn new(lineages: Vec<String>, sequences: Vec<Vec<u8>>) -> Result<Self> {
        let mut root = TreeNode::new(String::from("root"), 0, NodeType::Inner);
        let mut sequence_map: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
        let mut k_mer_map: Vec<Vec<usize>> = vec![Vec::new(); 2 << 15];
        let mut lineage_sequence_pairs = lineages.into_iter().zip_eq(sequences).collect_vec();
        lineage_sequence_pairs.sort_by(|(l1, _), (l2, _)| l1.cmp(l2));
        let mut confidence_idx = 0_usize;
        let _ = lineage_sequence_pairs
            .iter()
            .enumerate()
            .progress_with_style(
                ProgressStyle::with_template(
                    "[{elapsed_precise}] {bar:80.cyan/blue} {pos:>7}/{len:7}[ETA:{eta}] {msg}",
                )
                .unwrap()
                .progress_chars("##-"),
            )
            .with_message("Creating lineage tree and k-mer map...")
            .map(|(idx, (lineage, sequence))| -> Result<()> {
                let levels = lineage.split(',').collect_vec();
                let last_level_idx = levels.len() - 1;
                let mut current_node = &mut root;
                for (level, &label) in levels.iter().enumerate() {
                    let node_type = if level == last_level_idx {
                        NodeType::Taxon
                    } else {
                        NodeType::Inner
                    };
                    match &current_node.get_last_child_label() {
                        Some(name) => {
                            if name.as_str() != label {
                                current_node.add_child(TreeNode::new(
                                    label.to_owned(),
                                    confidence_idx,
                                    node_type,
                                ));
                            }
                            current_node.confidence_range.1 = confidence_idx + 1;
                        }
                        None => {
                            current_node.add_child(TreeNode::new(
                                label.to_owned(),
                                confidence_idx,
                                node_type,
                            ));
                            current_node.confidence_range.1 = confidence_idx + 1;
                        }
                    };
                    if level == last_level_idx {
                        confidence_idx += 1;
                    }
                    current_node = current_node.children.last_mut().unwrap();
                }
                current_node.add_child(TreeNode::new(
                    current_node.label.to_owned(),
                    confidence_idx - 1,
                    NodeType::Sequence,
                ));
                current_node.confidence_range.1 = confidence_idx;

                sequence_map.entry(sequence.clone()).or_default().push(idx);

                sequence.windows(8).for_each(|vals| {
                    if let Some(k_mer) = vals
                        .iter()
                        .enumerate()
                        .map(|(j, v)| map_four_to_two_bit_repr(*v).map(|c| c << (14 - j * 2)))
                        .fold_options(0_u16, |acc, c| acc | c)
                    {
                        k_mer_map[k_mer as usize].push(idx);
                    }
                });
                Ok(())
            })
            .collect::<Result<Vec<()>>>()?;
        root.confidence_range.1 = confidence_idx;
        let (sorted_lineages, _): (Vec<String>, Vec<Vec<u8>>) =
            lineage_sequence_pairs.into_iter().unzip();
        Ok(Self {
            root,
            lineages: sorted_lineages,
            sequences: sequence_map,
            k_mer_map: k_mer_map
                .into_par_iter()
                .map(|seqs| seqs.into_iter().unique().sorted().collect_vec())
                .collect(),
            num_tips: confidence_idx,
        })
    }

    pub fn print(&self) {
        self.root.print(0);
    }

    #[time("info")]
    pub fn save_to_file(&self, mut output: Box<dyn Write>) -> Result<()> {
        if log_enabled!(Level::Info) {
            println!("Writing database to file...");
        }
        bincode::serialize_into(&mut output, &self)?;
        Ok(())
    }

    pub fn load_from_file(path: &PathBuf) -> Result<Self> {
        if log_enabled!(Level::Info) {
            println!("Trying to read from database file...");
        }
        let mut file = File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;
        let decoded: Self = bincode::deserialize(&buffer)?;
        Ok(decoded)
    }

    pub fn get_shared_exact_match(&self, num_levels: usize, num_shared: usize) -> Vec<f64> {
        let mut values = vec![1.0; num_levels];
        values.push(1.0 / num_shared as f64);
        values
    }

    fn is_inner_taxon_node(&self, node: &TreeNode) -> bool {
        node.node_type == NodeType::Inner
    }

    fn is_taxon_leaf(&self, node: &TreeNode) -> bool {
        node.node_type == NodeType::Taxon
    }
}

#[derive(PartialEq, Eq, Debug, Serialize, Deserialize)]
enum NodeType {
    Inner,
    Taxon,
    Sequence,
}

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct TreeNode {
    label: String,
    confidence_range: (usize, usize),
    children: Vec<TreeNode>,
    node_type: NodeType,
}

impl TreeNode {
    fn new(label: String, confidence_idx: usize, node_type: NodeType) -> Self {
        Self {
            label,
            confidence_range: (confidence_idx, confidence_idx + 1),
            children: vec![],
            node_type,
        }
    }

    fn add_child(&mut self, child: TreeNode) {
        self.children.push(child);
    }

    fn get_last_child_label(&self) -> Option<&String> {
        match &self.children.last() {
            Some(c) => Some(&c.label),
            None => None,
        }
    }

    fn print(&self, depth: usize) {
        println!(
            "{}{} {:?}",
            "  ".repeat(depth),
            self.label,
            self.confidence_range
        );
        for child in &self.children {
            child.print(depth + 1);
        }
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::lineage::{EvaluationResult, Lineage, Tree};

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
        let confidence_values = &[0.1, 0.3, 0.4, 0.004, 0.004];
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
        let confidence_values = &[0.05, 0.1, 0.3, 0.4, 0.1, 0.004, 0.004];
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
        let confidence_values = &[0.004, 0.004, 0.004];
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
