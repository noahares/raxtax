use anyhow::{bail, Result};
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

use crate::utils;

pub struct Lineage<'a> {
    tree: &'a Tree,
    confidence_prefix_sum: Vec<f64>,
    confidence_vectors: Vec<(usize, Vec<f64>)>,
    rounding_factor: f64,
}

impl<'a> Lineage<'a> {
    pub fn new(tree: &'a Tree, confidence_values: &[f64]) -> Self {
        let mut confidence_prefix_sum = vec![0.0];
        confidence_prefix_sum.extend(confidence_values.iter().scan(0.0, |sum, i| {
            *sum += i;
            Some(*sum)
        }));
        Self {
            tree,
            confidence_prefix_sum,
            confidence_vectors: Vec::new(),
            rounding_factor: 10_u32.pow(utils::F64_OUTPUT_ACCURACY) as f64,
        }
    }

    #[time("info")]
    pub fn evaluate(mut self) -> Vec<(&'a String, Vec<f64>)> {
        self.eval_recurse(&self.tree.root, &[]);
        self.confidence_vectors
            .into_iter()
            .sorted_by(|a, b| b.1.iter().partial_cmp(a.1.iter()).unwrap())
            .map(|(idx, conf_values)| (&self.tree.lineages[idx], conf_values))
            .collect_vec()
    }

    fn get_confidence(&self, node: &TreeNode) -> f64 {
        self.confidence_prefix_sum[node.confidence_range.1]
            - self.confidence_prefix_sum[node.confidence_range.0]
    }

    fn eval_recurse(&mut self, node: &TreeNode, confidence_prefix: &[f64]) {
        let mut no_child_significant = true;
        for c in &node.children {
            let mut conf_prefix = confidence_prefix.to_vec();
            let child_conf =
                (self.get_confidence(c) * self.rounding_factor).round() / self.rounding_factor;
            if child_conf == 0.0 {
                continue;
            }
            no_child_significant = false;
            conf_prefix.push(child_conf);
            self.eval_recurse(c, &conf_prefix);
            if c.level == self.tree.num_levels - 1 {
                self.confidence_vectors
                    .push((c.confidence_range.0, conf_prefix));
            }
        }
        if no_child_significant && node.level < self.tree.num_levels - 1 {
            let mut conf_prefix = confidence_prefix.to_vec();
            let mut current_node = node;
            while current_node.level < self.tree.num_levels - 1 {
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
            }
            self.confidence_vectors
                .push((current_node.confidence_range.0, conf_prefix));
        }
    }
}

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct Tree {
    root: TreeNode,
    pub lineages: Vec<String>,
    pub sequences: HashMap<Vec<u8>, Vec<usize>>,
    pub k_mer_map: Vec<Vec<usize>>,
    pub num_levels: usize,
    pub num_tips: usize,
}

impl Tree {
    #[time("info")]
    pub fn new(lineages: Vec<String>, sequences: Vec<Vec<u8>>) -> Result<Self> {
        let mut root = TreeNode::new(String::from("root"), 0, 0);
        let mut sequence_map: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
        let mut k_mer_map: Vec<Vec<usize>> = vec![Vec::new(); 2 << 15];
        let mut lineage_sequence_pairs = lineages.into_iter().zip_eq(sequences).collect_vec();
        lineage_sequence_pairs.sort_by(|(l1, _), (l2, _)| l1.cmp(l2));
        let num_levels = lineage_sequence_pairs[0].0.split(',').count() + 1;
        let mut confidence_idx = 0_usize;
        let _ = lineage_sequence_pairs.iter().enumerate()
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
            if levels.len() + 1 != num_levels {
                    bail!("Number of taxonomical levels do not match for {lineage}, expexted {}, found {}", num_levels - 1, levels.len())
            }
            let mut current_node = &mut root;
            for (level, &label) in levels.iter().enumerate() {
                match &current_node.get_last_child_label() {
                    Some(name) => {
                        if name.as_str() != label {
                            current_node.add_child(TreeNode::new(
                                label.to_owned(),
                                level + 1,
                                confidence_idx,
                            ));
                        }
                        current_node.confidence_range.1 = confidence_idx + 1;
                    }
                    None => {
                        current_node.add_child(TreeNode::new(
                            label.to_owned(),
                            level + 1,
                            confidence_idx,
                        ));
                        current_node.confidence_range.1 = confidence_idx + 1;
                    }
                };
                if level + 2 == num_levels {
                    confidence_idx += 1;
                }
                current_node = current_node.children.last_mut().unwrap();
            }
            current_node.add_child(TreeNode::new(
                current_node.label.to_owned(),
                num_levels,
                confidence_idx - 1,
            ));
            current_node.confidence_range.1 = confidence_idx;

            sequence_map.entry(sequence.clone()).or_default().push(idx);

            let mut k_mer: u16 = 0;
            sequence[0..8]
                .iter()
                .enumerate()
                .for_each(|(j, c)| k_mer |= (*c as u16) << (14 - j * 2));
            k_mer_map[k_mer as usize].push(idx);
            sequence[8..].iter().for_each(|c| {
                k_mer = (k_mer << 2) | *c as u16;
                k_mer_map[k_mer as usize].push(idx);
            });
            Ok(())
        }).collect::<Result<Vec<()>>>()?;
        root.confidence_range.1 = confidence_idx;
        let (sorted_lineages, _): (Vec<String>, Vec<Vec<u8>>) =
            lineage_sequence_pairs.into_iter().unzip();
        dbg!(num_levels);
        Ok(Self {
            root,
            lineages: sorted_lineages,
            sequences: sequence_map,
            k_mer_map: k_mer_map
                .into_par_iter()
                .map(|seqs| seqs.into_iter().unique().sorted().collect_vec())
                .collect(),
            num_levels,
            num_tips: confidence_idx,
        })
    }

    pub fn print(&self) {
        self.root.print();
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
            println!("Reading from database file...");
        }
        let mut file = File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;
        let decoded: Self = bincode::deserialize(&buffer)?;
        Ok(decoded)
    }

    pub fn get_shared_exact_match(&self, num_shared: usize) -> Vec<f64> {
        let mut values = vec![1.0; self.num_levels - 2];
        values.push(1.0 / num_shared as f64);
        values
    }
}

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct TreeNode {
    label: String,
    level: usize,
    confidence_range: (usize, usize),
    children: Vec<TreeNode>,
}

impl TreeNode {
    fn new(label: String, level: usize, confidence_idx: usize) -> Self {
        Self {
            label,
            level,
            confidence_range: (confidence_idx, confidence_idx + 1),
            children: vec![],
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

    fn print(&self) {
        println!(
            "{}{} {:?}",
            "  ".repeat(self.level),
            self.label,
            self.confidence_range
        );
        for child in &self.children {
            child.print();
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::lineage::{Lineage, Tree};

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
        let lineage = Lineage::new(&tree, confidence_values);
        let result = lineage.evaluate();
        dbg!(result);
        assert_eq!(0, 1);
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
        let lineage = Lineage::new(&tree, confidence_values);
        let result = lineage.evaluate();
        dbg!(result);
        assert_eq!(0, 1);
    }
}
