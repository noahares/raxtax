use anyhow::{bail, Result};
use itertools::Itertools;
use log::{debug, Level};
use logging_timer::{time, timer};
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    fs,
    path::PathBuf,
};

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct LookupTables {
    pub species_meta_vec: Vec<SpeciesMeta>,
    pub level_name_maps: Vec<Vec<String>>,
    pub level_hierarchy_maps: Vec<Vec<Vec<usize>>>,
    pub k_mer_map: Vec<Vec<usize>>,
}

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SpeciesMeta {
    pub name: String,
    // 0 -> phylum_idx,
    // 1 -> class_idx,
    // 2 -> order_idx,
    // 3 -> family_idx,
    // 4 -> genus_idx,
    pub indices: [usize; 5],
}

pub fn parse_reference_fasta_file(sequence_path: &PathBuf) -> Result<LookupTables> {
    let fasta_str = fs::read_to_string(sequence_path)?;
    parse_reference_fasta_str(&fasta_str)
}

#[time("info")]
pub fn parse_reference_fasta_str(fasta_str: &str) -> Result<LookupTables> {
    // Level 0: Phylum
    // Level 1: Class
    // Level 2: Order
    // Level 3: Family
    // Level 4: Genus
    let mut level_sets: [HashSet<String>; 5] = Default::default();
    const ARRAY_REPEAT_VALUE: std::vec::Vec<usize> = Vec::new();
    let mut k_mer_map: Vec<Vec<usize>> = vec![ARRAY_REPEAT_VALUE; 2 << 15];
    let labels = {
        let _tmr = timer!(Level::Info; "Read file and create k-mer mapping");
        let lines: Vec<String> = fasta_str
            .lines()
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty())
            .collect();
        if !lines[0].starts_with(['>', ';']) {
            bail!("Not a valid FASTA file")
        }
        let mut labels: Vec<Vec<String>> = Vec::new();
        let mut current_sequence = Vec::<u8>::new();
        let mut idx: usize = 0;

        // create label and sequence vectors
        for line in lines {
            if line.starts_with(';') {
                continue;
            }
            if let Some(label) = line.strip_prefix('>') {
                let taxon_info = label.rsplitn(7, '|').map(|s| s.to_string()).collect_vec();
                taxon_info[1..6]
                    .iter()
                    .zip_eq((0..level_sets.len()).rev())
                    .for_each(|(name, idx)| {
                        level_sets[idx].insert(name.to_string());
                    });
                labels.push(taxon_info[0..6].to_vec());
                if labels.len() > 1 {
                    let mut k_mer: u16 = 0;
                    // dbg!(&current_sequence);
                    current_sequence[0..8]
                        .iter()
                        .enumerate()
                        .for_each(|(j, c)| k_mer |= (*c as u16) << (14 - j * 2));
                    k_mer_map[k_mer as usize].push(idx);
                    // println!("{k_mer:#b}");
                    current_sequence[8..].iter().for_each(|c| {
                        k_mer = (k_mer << 2) | *c as u16;
                        // println!("{k_mer:#b}");
                        k_mer_map[k_mer as usize].push(idx);
                    });
                    current_sequence = Vec::new();
                    idx += 1;
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
        let mut k_mer: u16 = 0;
        // dbg!(&current_sequence);
        current_sequence[0..8]
            .iter()
            .enumerate()
            .for_each(|(j, c)| k_mer |= (*c as u16) << (14 - j * 2));
        k_mer_map[k_mer as usize].push(idx);
        // println!("{k_mer:#b}");
        current_sequence[8..].iter().for_each(|c| {
            k_mer = (k_mer << 2) | *c as u16;
            // println!("{k_mer:#b}");
            k_mer_map[k_mer as usize].push(idx);
        });
        labels
    };

    debug!(
        "Average number of sequences per kmer: {}",
        k_mer_map.iter().map(|e| e.len() as f64).sum::<f64>() / k_mer_map.len() as f64
    );

    // create mapping to index for each taxonomical level
    let level_name_maps = level_sets
        .into_iter()
        .map(|set| set.into_iter().sorted().collect_vec())
        .collect::<Vec<Vec<String>>>();
    // dbg!(&level_maps);
    debug!("Unique Phyla: {}", level_name_maps[0].len());
    debug!("Unique Classes: {}", level_name_maps[1].len());
    debug!("Unique Orders: {}", level_name_maps[2].len());
    debug!("Unique Families: {}", level_name_maps[3].len());
    debug!("Unique Genus: {}", level_name_maps[4].len());
    debug!("Unique Species: {}", labels.iter().unique().count());
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

    // dbg!(&level_rev_maps);
    // create data structure for mapping k-mers to sequences as well as between taxonomical levels
    let mut species_meta_vec: Vec<SpeciesMeta> = Vec::new();
    let mut level_hierarchy_maps: Vec<Vec<HashSet<usize>>> = Vec::new();
    for level in &level_name_maps {
        level_hierarchy_maps.push(vec![HashSet::new(); level.len()]);
    }
    labels.into_iter().enumerate().for_each(|(i, taxon_info)| {
        // let taxon_info = label.rsplitn(7, '|').collect_vec();
        let species_meta = SpeciesMeta {
            name: taxon_info.iter().rev().join("|"),
            indices: [
                level_rev_maps[0][taxon_info[5].as_str()],
                level_rev_maps[1][taxon_info[4].as_str()],
                level_rev_maps[2][taxon_info[3].as_str()],
                level_rev_maps[3][taxon_info[2].as_str()],
                level_rev_maps[4][taxon_info[1].as_str()],
            ],
        };
        level_hierarchy_maps[0][species_meta.indices[0]].insert(species_meta.indices[1]);
        level_hierarchy_maps[1][species_meta.indices[1]].insert(species_meta.indices[2]);
        level_hierarchy_maps[2][species_meta.indices[2]].insert(species_meta.indices[3]);
        level_hierarchy_maps[3][species_meta.indices[3]].insert(species_meta.indices[4]);
        level_hierarchy_maps[4][species_meta.indices[4]].insert(i);
        species_meta_vec.push(species_meta);
    });
    // dbg!(level_hierarchy_maps);
    Ok(LookupTables {
        species_meta_vec,
        level_name_maps,
        level_hierarchy_maps: level_hierarchy_maps
            .into_iter()
            .map(|level| {
                level
                    .into_iter()
                    .map(|item| item.into_iter().collect_vec())
                    .collect_vec()
            })
            .collect_vec(),
        k_mer_map,
    })
}

pub fn parse_query_fasta_file(sequence_path: &PathBuf) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    let fasta_str = fs::read_to_string(sequence_path)?;
    parse_query_fasta_str(&fasta_str)
}

#[time("info")]
pub fn parse_query_fasta_str(fasta_str: &str) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    let lines: Vec<String> = fasta_str
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();
    if !lines[0].starts_with(['>', ';']) {
        bail!("Not a valid FASTA file")
    }
    let mut labels: Vec<String> = Vec::new();
    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut current_sequence = Vec::<u8>::new();

    // create label and sequence vectors
    for line in lines {
        if line.starts_with(';') {
            continue;
        }
        if let Some(label) = line.strip_prefix('>') {
            labels.push(label.to_string());
            if labels.len() > 1 {
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
    use crate::parser::{LookupTables, SpeciesMeta};

    use super::{parse_reference_fasta_str, parse_query_fasta_str};

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
>Badabing|Badabum|Phylum2|Class3|Order3|Family4|Genus4|Species5
ATACGCTTTGCGT";
        let LookupTables {
            species_meta_vec,
            k_mer_map,
            ..
        } = parse_reference_fasta_str(fasta_str).unwrap();
        // for (k, v) in k_mer_map.iter().enumerate() {
        //     if !v.is_empty() {
        //         println!("{k:b}:\n {v:?}");
        //     }
        // }
        assert_eq!(k_mer_map[0b1010111111110_usize], &[0]);
        assert_eq!(k_mer_map[0b11000110011111_usize], &[1, 3, 4]);
        assert_eq!(k_mer_map[0b1001111111100110_usize], &[3, 4]);
        assert_eq!(k_mer_map[0b110011100111010_usize], &[2]);
        assert_eq!(
            species_meta_vec,
            &[
                SpeciesMeta {
                    name: "Phylum1|Class1|Order1|Family1|Genus1|Species1"
                        .to_string(),
                    indices: [0, 0, 0, 0, 0]
                },
                SpeciesMeta {
                    name: "Phylum1|Class1|Order1|Family1|Genus1|Species2"
                        .to_string(),
                    indices: [0, 0, 0, 0, 0]
                },
                SpeciesMeta {
                    name: "Phylum1|Class1|Order4|Family5|Genus2|Species3"
                        .to_string(),
                    indices: [0, 0, 3, 3, 1]
                },
                SpeciesMeta {
                    name: "Phylum1|Class2|Order2|Family3|Genus3|Species6"
                        .to_string(),
                    indices: [0, 1, 1, 1, 2]
                },
                SpeciesMeta {
                    name: "Phylum2|Class3|Order3|Family4|Genus4|Species5"
                        .to_string(),
                    indices: [1, 2, 2, 2, 3]
                }
            ]
        );
        // assert_eq!(0, 1);
    }

    #[test]
    fn test_query_parser() {
        let fasta_str = r">label1
AAACCCTTTGGGA
>label2
ATACGCTTTGGGA
>label3
ATCCGCTATGGGA";
    let (labels, sequences) = parse_query_fasta_str(fasta_str).unwrap();
    assert_eq!(sequences[0], &[0, 0, 0, 1, 1, 1, 3, 3, 3, 2, 2, 2, 0]);
    }
}
