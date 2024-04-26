use anyhow::{bail, Context, Result};
use indicatif::{ProgressIterator, ProgressStyle};
use itertools::Itertools;
use log::{info, Level};
use logging_timer::{time, timer};
use rayon::prelude::*;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    path::PathBuf,
};

use crate::utils;

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct LookupTables {
    pub labels: Vec<String>,
    pub sequences: HashMap<Vec<u8>, Vec<usize>>,
    pub level_hierarchy_maps: Vec<Vec<Vec<usize>>>,
    pub k_mer_map: Vec<Vec<usize>>,
    pub sequence_species_map: Vec<usize>,
}

pub fn parse_reference_fasta_file(sequence_path: &PathBuf) -> Result<LookupTables> {
    let mut fasta_str = String::new();
    let _ = utils::get_reader(sequence_path)?.read_to_string(&mut fasta_str);
    parse_reference_fasta_str(&fasta_str)
}

#[time("info")]
pub fn parse_reference_fasta_str(fasta_str: &str) -> Result<LookupTables> {
    // Level 0: Phylum
    // Level 1: Class
    // Level 2: Order
    // Level 3: Family
    // Level 4: Genus
    // Level 5: Species
    let mut level_sets: [HashSet<String>; 6] = Default::default();
    let mut k_mer_map: Vec<Vec<usize>> = vec![Vec::new(); 2 << 15];
    let regex = Regex::new(r".*tax=p:(.*?),c:(.*?),o:(.*?),f:(.*?),g:(.*?),s:(.*?);")?;
    let (labels, sequences) = {
        let _tmr = timer!(Level::Info; "Read file and create k-mer mapping");
        let lines: Vec<String> = fasta_str
            .lines()
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty() && !l.starts_with(';'))
            .collect();
        if !lines[0].starts_with('>') {
            bail!("Not a valid FASTA file")
        }
        let mut labels: Vec<String> = Vec::new();
        let mut sequences: Vec<Vec<u8>> = Vec::new();
        let mut current_sequence = Vec::<u8>::new();

        // create label and sequence vectors
        lines
            .into_iter()
            .progress_with_style(
                ProgressStyle::with_template(
                    "[{elapsed_precise}] {bar:80.cyan/blue} {pos:>7}/{len:7}[ETA:{eta}] {msg}",
                )
                .unwrap()
                .progress_chars("##-"),
            )
            .with_message("Parsing Reference...")
            .map(|line| -> Result<()> {
                if let Some(label) = line.strip_prefix('>') {
                    let (_, pieces): (&str, [&str; 6]) = regex
                        .captures(label)
                        .map(|caps| caps.extract())
                        .context(format!(
                            "Unexpected taxonomical annotation detected in {}",
                            label
                        ))?;
                    labels.push(pieces.join("|"));
                    if !current_sequence.is_empty() {
                        sequences.push(current_sequence.clone());
                        current_sequence = Vec::new();
                    }
                } else {
                    current_sequence.extend(line.chars().map(|c| -> u8 {
                        match c.to_ascii_uppercase() {
                            'A' => 0b00,
                            'C' => 0b01,
                            'G' => 0b10,
                            'T' => 0b11,
                            _ => panic!("Unexpected character: {}", c),
                        }
                    }))
                }
                Ok(())
            })
            .collect::<Result<Vec<()>>>()?;
        sequences.push(current_sequence);
        let mut ordered_labels: Vec<String> = Vec::new();
        let mut sequence_map: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
        labels
            .iter()
            .zip_eq(sequences)
            .sorted_by_key(|(l, _)| *l)
            .enumerate()
            .progress_with_style(
                ProgressStyle::with_template(
                    "[{elapsed_precise}] {bar:80.cyan/blue} {pos:>7}/{len:7}[ETA:{eta}] {msg}",
                )
                .unwrap()
                .progress_chars("##-"),
            )
            .with_message("Creating k-mer map...")
            .for_each(|(idx, (l, s))| {
                let taxon_info = l.split('|').collect_vec();
                taxon_info[0..6]
                    .iter()
                    .zip_eq(0..level_sets.len())
                    .for_each(|(name, idx)| {
                        level_sets[idx].insert(name.to_string());
                    });
                ordered_labels.push(l.to_string());
                sequence_map.entry(s.clone()).or_default().push(idx);
                let mut k_mer: u16 = 0;
                s[0..8]
                    .iter()
                    .enumerate()
                    .for_each(|(j, c)| k_mer |= (*c as u16) << (14 - j * 2));
                k_mer_map[k_mer as usize].push(idx);
                s[8..].iter().for_each(|c| {
                    k_mer = (k_mer << 2) | *c as u16;
                    k_mer_map[k_mer as usize].push(idx);
                });
            });
        (ordered_labels, sequence_map)
    };

    // create mapping to index for each taxonomical level
    let level_name_maps = level_sets
        .into_iter()
        .map(|set| set.into_iter().sorted().collect_vec())
        .collect::<Vec<Vec<String>>>();
    info!("Unique Phyla: {}", level_name_maps[0].len());
    info!("Unique Classes: {}", level_name_maps[1].len());
    info!("Unique Orders: {}", level_name_maps[2].len());
    info!("Unique Families: {}", level_name_maps[3].len());
    info!("Unique Genus: {}", level_name_maps[4].len());
    info!("Unique Species: {}", level_name_maps[5].len());
    info!("Unique Sequences: {}", labels.len());
    // need reverse mapping for second parsing of labels and sequences to build data structure
    let level_rev_maps = level_name_maps
        .iter()
        .map(|map| {
            map.iter()
                .enumerate()
                .map(|(i, s)| (s.as_str(), i))
                .collect::<HashMap<&str, usize>>()
        })
        .collect_vec();

    // create data structure for mapping k-mers to sequences as well as between taxonomical levels
    let mut level_hierarchy_maps: Vec<Vec<HashSet<usize>>> = Vec::new();
    for level in &level_name_maps {
        level_hierarchy_maps.push(vec![HashSet::new(); level.len()]);
    }
    labels.iter().enumerate().for_each(|(i, label)| {
        let taxon_info = label.split('|').collect_vec();
        let phylum_idx = level_rev_maps[0][taxon_info[0]];
        let class_idx = level_rev_maps[1][taxon_info[1]];
        let order_idx = level_rev_maps[2][taxon_info[2]];
        let family_idx = level_rev_maps[3][taxon_info[3]];
        let genus_idx = level_rev_maps[4][taxon_info[4]];
        let species_idx = level_rev_maps[5][taxon_info[5]];

        level_hierarchy_maps[0][phylum_idx].insert(class_idx);
        level_hierarchy_maps[1][class_idx].insert(order_idx);
        level_hierarchy_maps[2][order_idx].insert(family_idx);
        level_hierarchy_maps[3][family_idx].insert(genus_idx);
        level_hierarchy_maps[4][genus_idx].insert(species_idx);
        level_hierarchy_maps[5][species_idx].insert(i);
    });

    let mut sequence_species_map = vec![0; labels.len()];
    for (species_id, sequences) in level_hierarchy_maps[5].iter().enumerate() {
        for sequence_id in sequences {
            sequence_species_map[*sequence_id] = species_id;
        }
    }

    Ok(LookupTables {
        labels,
        sequences,
        level_hierarchy_maps: level_hierarchy_maps
            .into_iter()
            .map(|level| {
                level
                    .into_iter()
                    .map(|item| item.into_iter().sorted().collect_vec())
                    .collect_vec()
            })
            .collect_vec(),
        k_mer_map: k_mer_map
            .into_par_iter()
            .map(|seqs| seqs.into_iter().unique().sorted().collect_vec())
            .collect(),
        sequence_species_map,
    })
}

pub fn parse_query_fasta_file(sequence_path: &PathBuf) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    let mut fasta_str = String::new();
    let _ = utils::get_reader(sequence_path)?.read_to_string(&mut fasta_str);
    parse_query_fasta_str(&fasta_str)
}

#[time("info")]
pub fn parse_query_fasta_str(fasta_str: &str) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    let lines: Vec<String> = fasta_str
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty() && !l.starts_with(';'))
        .collect();
    if !lines[0].starts_with('>') {
        bail!("Not a valid FASTA file")
    }
    let mut labels: Vec<String> = Vec::new();
    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut current_sequence = Vec::<u8>::new();

    // create label and sequence vectors
    for line in lines {
        if let Some(label) = line.strip_prefix('>') {
            labels.push(label.to_string());
            if !current_sequence.is_empty() {
                sequences.push(current_sequence);
                current_sequence = Vec::new();
            }
        } else {
            current_sequence.extend(line.chars().map(|c| -> u8 {
                match c {
                    'A' => 0b00,
                    'C' => 0b01,
                    'G' => 0b10,
                    'T' => 0b11,
                    _ => panic!("Unexpected character: {}", c),
                }
            }))
        }
    }
    sequences.push(current_sequence);
    Ok((labels, sequences))
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::parser::LookupTables;

    use super::{parse_query_fasta_str, parse_reference_fasta_str};

    #[test]
    fn test_str_parser() {
        let fasta_str = r">Badabing|Badabum|Phylum1|Class1|Order1|Family1|Genus1|Species1
AAACCCTTTGGGA
>Badabing|Badabum|Phylum1|Class1|Order1|Family1|Genus1|Species2
ATACGCTTTGGGA
>Badabing|Badabum|Phylum1|Class1|Order4|Family5|Genus2|Species3
ATCCGCTATGGGA
>Badabing|Badabum|Phylum1|Class2|Order2|Family3|Genus3|Species6
ATACGCTTTGCGT
>Badabing|Badabum|Phylum1|Class1|Order1|Family1|Genus1|Species2
GTGCGCTATGCGA
>Badabing|Badabum|Phylum2|Class3|Order3|Family4|Genus4|Species5
ATACGCTTTGCGT";
        let LookupTables { k_mer_map, .. } = parse_reference_fasta_str(fasta_str).unwrap();
        for (k, v) in k_mer_map.iter().enumerate() {
            if !v.is_empty() {
                println!("{k:b}:\n {v:?}");
            }
        }
        assert_eq!(
            k_mer_map[0b1010111111110_usize].iter().collect_vec(),
            &[&0_usize]
        );
        assert_eq!(
            k_mer_map[0b11000110011111_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&1, &4, &5]
        );
        assert_eq!(
            k_mer_map[0b110011100111010_usize].iter().collect_vec(),
            &[&3]
        );
    }

    #[test]
    fn test_query_parser() {
        let fasta_str = r">label1
AAACCCTTTGGGA
>label2
ATACGCTTTGGGA
>label3
ATCCGCTATGGGA";
        let (_, sequences) = parse_query_fasta_str(fasta_str).unwrap();
        assert_eq!(sequences[0], &[0, 0, 0, 1, 1, 1, 3, 3, 3, 2, 2, 2, 0]);
    }

    #[test]
    fn test_kmers() {
        let fasta_str = r">Badabing|Badabum|Phylum1|Class1|Order1|Family1|Genus1|Species1
AAACCCCGT
>Badabing|Badabum|Phylum1|Class1|Order1|Family1|Genus1|Species1
TAACCCCGG
>Badabing|Badabum|Phylum1|Class1|Order1|Family1|Genus2|Species3
TTTAAAACC
>Badabing|Badabum|Phylum1|Class1|Order1|Family1|Genus2|Species3
TTTAAAACA
>Badabing|Badabum|Phylum1|Class2|Order2|Family2|Genus3|Species4
AAACCCCGG";
        let LookupTables { k_mer_map, .. } = parse_reference_fasta_str(fasta_str).unwrap();
        // for (k, v) in k_mer_map.iter().enumerate() {
        //     if !v.is_empty() {
        //         println!("{k:b}:\n {v:?}");
        //     }
        // }
        assert_eq!(
            k_mer_map[0b101010110_usize].iter().sorted().collect_vec(),
            &[&0, &4]
        );
        assert_eq!(
            k_mer_map[0b10101011010_usize].iter().sorted().collect_vec(),
            &[&1, &4]
        );
        assert_eq!(
            k_mer_map[0b10101011011_usize].iter().sorted().collect_vec(),
            &[&0]
        );
        assert_eq!(
            k_mer_map[0b1100000101010110_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&1]
        );
        assert_eq!(
            k_mer_map[0b1111000000000101_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&2]
        );
        assert_eq!(
            k_mer_map[0b1111000000000101_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&2]
        );
        assert_eq!(
            k_mer_map[0b1111110000000001_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&2, &3]
        );
    }
}
