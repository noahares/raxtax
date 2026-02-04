use anyhow::Result;
use std::{
    fs::File,
    io::{BufWriter, Read},
    path::PathBuf,
};

use ahash::HashMap;
use indicatif::{ProgressIterator, ProgressStyle};
use itertools::Itertools;
use log::{log_enabled, Level};
use logging_timer::time;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::utils::map_four_to_two_bit_repr;

#[cfg(feature = "huge_db")]
type IndexType = usize;

#[cfg(not(feature = "huge_db"))]
type IndexType = u32;

#[cfg(not(feature = "huge_db"))]
fn check_lineage_size(db_size: usize) {
    assert!(
        u32::try_from(db_size).is_ok(),
        "Too many database sequences to run with 32-bit indices!\n
            Re-compile raxtax with '--features huge_db' to enable usize indices."
    );
}

#[cfg(feature = "huge_db")]
fn check_lineage_size() {}

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct Tree {
    pub root: Node,
    pub lineages: Vec<String>,
    pub bins: Vec<String>,
    pub sequences: HashMap<Vec<u8>, Vec<IndexType>>,
    pub k_mer_map: Vec<Vec<IndexType>>,
    pub bin_idx_to_lineage_idxs: Vec<Vec<usize>>,
    pub lineage_idx_to_bin_idx: Vec<Option<usize>>,
    pub num_tips: usize,
}

impl Tree {
    #[time("debug", "Tree::{}")]
    pub fn new(labels: Vec<(String, Option<String>)>, sequences: Vec<Vec<u8>>) -> Result<Self> {
        check_lineage_size(labels.len());
        let mut root = Node::new(String::from("root"), 0, NodeType::Inner);
        let mut sequence_map: HashMap<Vec<u8>, Vec<IndexType>> =
            sequences.iter().map(|s| (s.clone(), Vec::new())).collect();
        let mut k_mer_map: Vec<Vec<IndexType>> = vec![Vec::new(); 2 << 15];
        let mut lineage_sequence_pairs = labels.into_iter().zip_eq(sequences).collect_vec();
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
            .map(|(idx, ((lineage, _), sequence))| -> Result<()> {
                let levels = lineage.split(',').collect_vec();
                let last_level_idx = levels.len() - 1;
                let mut current_node = &mut root;
                for (level, label) in levels.into_iter().enumerate() {
                    let node_type = if level == last_level_idx {
                        NodeType::Taxon
                    } else {
                        NodeType::Inner
                    };
                    match &current_node.get_last_child_label() {
                        Some(name) => {
                            if name.as_str() != label {
                                current_node.add_child(Node::new(
                                    label.to_string(),
                                    confidence_idx,
                                    node_type,
                                ));
                            }
                            current_node.confidence_range.1 = confidence_idx + 1;
                        }
                        None => {
                            current_node.add_child(Node::new(
                                label.to_string(),
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
                current_node.add_child(Node::new(
                    current_node.label.clone(),
                    confidence_idx - 1,
                    NodeType::Sequence,
                ));
                current_node.confidence_range.1 = confidence_idx;

                sequence_map
                    .get_mut(sequence)
                    .unwrap()
                    .push(idx as IndexType);

                sequence.windows(8).for_each(|vals| {
                    if let Some(k_mer) = vals
                        .iter()
                        .enumerate()
                        .map(|(j, v)| map_four_to_two_bit_repr(*v).map(|c| c << (14 - j * 2)))
                        .fold_options(0_u16, |acc, c| acc | c)
                    {
                        k_mer_map[k_mer as usize].push(idx as IndexType);
                    }
                });
                Ok(())
            })
            .collect::<Result<Vec<()>>>()?;
        root.confidence_range.1 = confidence_idx;
        let (sorted_lineages, _): (Vec<(String, Option<String>)>, Vec<Vec<u8>>) =
            lineage_sequence_pairs.into_iter().unzip();

        let mut bin_idx_to_lineage_idxs: Vec<Vec<usize>> = Vec::new();
        let bin_id_idx_pairs = sorted_lineages
            .iter()
            .enumerate()
            .flat_map(|(idx, (_, bin_id))| bin_id.clone().map(|b| (b, idx)))
            .sorted()
            .collect_vec();
        let mut lineage_idx_to_bin_idx: Vec<Option<usize>> = vec![None; confidence_idx];
        if let Some(pair) = bin_id_idx_pairs.first() {
            let mut current_bin = &pair.0;
            let mut current_idxs = Vec::new();
            let mut bin_idx = 0_usize;
            for (bin, idx) in bin_id_idx_pairs.iter() {
                if *bin != *current_bin {
                    bin_idx_to_lineage_idxs.push(current_idxs.clone());
                    current_bin = bin;
                    bin_idx += 1;
                    current_idxs.clear();
                }
                current_idxs.push(*idx);
                lineage_idx_to_bin_idx[*idx] = Some(bin_idx);
            }
            bin_idx_to_lineage_idxs.push(current_idxs);
        };
        let (bins, _): (Vec<String>, Vec<usize>) = bin_id_idx_pairs.into_iter().unzip();
        let (lineages, _): (Vec<String>, Vec<Option<String>>) = sorted_lineages.into_iter().unzip();
        Ok(Self {
            root,
            lineages,
            bins: bins.into_iter().unique().collect_vec(),
            sequences: sequence_map,
            k_mer_map: k_mer_map
                .into_par_iter()
                .map(|seqs| seqs.into_iter().unique().sorted().collect_vec())
                .collect(),
            bin_idx_to_lineage_idxs,
            lineage_idx_to_bin_idx,
            num_tips: confidence_idx,
        })
    }

    pub fn print(&self) {
        self.root.print(0);
    }

    #[time("debug")]
    pub fn save_to_file(&self, mut output: BufWriter<File>) -> Result<()> {
        if log_enabled!(Level::Info) {
            eprintln!("[INFO ] Writing database to file...");
        }
        bincode::serialize_into(&mut output, &self)?;
        Ok(())
    }

    pub fn load_from_file(path: &PathBuf) -> Result<Self> {
        if log_enabled!(Level::Info) {
            eprintln!("[INFO ] Trying to read from database file...");
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

    pub fn is_inner_taxon_node(&self, node: &Node) -> bool {
        node.node_type == NodeType::Inner
    }

    pub fn is_taxon_leaf(&self, node: &Node) -> bool {
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
pub struct Node {
    pub label: String,
    pub confidence_range: (usize, usize),
    pub children: Vec<Node>,
    node_type: NodeType,
}

impl Node {
    fn new(label: String, confidence_idx: usize, node_type: NodeType) -> Self {
        Self {
            label,
            confidence_range: (confidence_idx, confidence_idx + 1),
            children: vec![],
            node_type,
        }
    }

    fn add_child(&mut self, child: Node) {
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
