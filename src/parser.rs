use anyhow::{bail, Context, Result};
use indicatif::{ProgressIterator, ProgressStyle};
use log::Level;
use logging_timer::{time, timer};
use regex::Regex;
use std::{io::Read, path::PathBuf};

use crate::{lineage, utils};

#[time("info")]
pub fn parse_reference_fasta_file(sequence_path: &PathBuf) -> Result<(bool, lineage::Tree)> {
    if let Ok(tree) = lineage::Tree::load_from_file(sequence_path) {
        return Ok((false, tree));
    }
    let mut fasta_str = String::new();
    let _ = utils::get_reader(sequence_path)?.read_to_string(&mut fasta_str);
    Ok((true, parse_reference_fasta_str(&fasta_str)?))
}

pub fn parse_reference_fasta_str(fasta_str: &str) -> Result<lineage::Tree> {
    if fasta_str.is_empty() {
        bail!("File is empty")
    }
    let regex = Regex::new(r"tax=([^;]+);")?;
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
                    let lineage = regex
                        .captures(label)
                        .context(format!(
                            "Unexpected taxonomical annotation detected in label {}",
                            label
                        ))?
                        .get(1)
                        .context(format!("No taxonomic string found in label {}", label))?
                        .as_str()
                        .to_owned();
                    labels.push(lineage);
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
        if labels.len() != sequences.len() {
            bail!("Number of sequences does not match number of labels")
        }
        (labels, sequences)
    };
    lineage::Tree::new(labels, sequences)
}

pub fn parse_query_fasta_file(sequence_path: &PathBuf) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    let mut fasta_str = String::new();
    let _ = utils::get_reader(sequence_path)?.read_to_string(&mut fasta_str);
    parse_query_fasta_str(&fasta_str)
}

#[time("info")]
pub fn parse_query_fasta_str(fasta_str: &str) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    if fasta_str.is_empty() {
        bail!("File is empty")
    }
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
                match c.to_ascii_uppercase() {
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
    if labels.len() != sequences.len() {
        bail!("Number of sequences does not match number of labels")
    }
    Ok((labels, sequences))
}

// #[cfg(test)]
// mod tests {
//     use itertools::Itertools;
//
//
//     use super::{parse_query_fasta_str, parse_reference_fasta_str};
//
//     #[test]
//     fn test_str_parser() {
//         let fasta_str = r">Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1;
// AAACCCTTTGGGA
// >Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species2;
// ATACGCTTTGGGA
// >Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order4,f:Family5,g:Genus2,s:Species3;
// ATCCGCTATGGGA
// >Badabing|Badabum;tax=p:Phylum1,c:Class2,o:Order2,f:Family3,g:Genus3,s:Species6;
// ATACGCTTTGCGT
// >Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species2;
// GTGCGCTATGCGA
// >Badabing|Badabum;tax=p:Phylum2,c:Class3,o:Order3,f:Family4,g:Genus4,s:Species5;
// ATACGCTTTGCGT";
//         let LookupTables { k_mer_map, .. } = parse_reference_fasta_str(fasta_str).unwrap();
//         for (k, v) in k_mer_map.iter().enumerate() {
//             if !v.is_empty() {
//                 println!("{k:b}:\n {v:?}");
//             }
//         }
//         assert_eq!(
//             k_mer_map[0b1010111111110_usize].iter().collect_vec(),
//             &[&0_usize]
//         );
//         assert_eq!(
//             k_mer_map[0b11000110011111_usize]
//                 .iter()
//                 .sorted()
//                 .collect_vec(),
//             &[&1, &4, &5]
//         );
//         assert_eq!(
//             k_mer_map[0b110011100111010_usize].iter().collect_vec(),
//             &[&3]
//         );
//     }
//
//     #[test]
//     fn test_query_parser() {
//         let fasta_str = r">label1
// AAACCCTTTGGGA
// >label2
// ATACGCTTTGGGA
// >label3
// ATCCGCTATGGGA";
//         let (_, sequences) = parse_query_fasta_str(fasta_str).unwrap();
//         assert_eq!(sequences[0], &[0, 0, 0, 1, 1, 1, 3, 3, 3, 2, 2, 2, 0]);
//     }
//
//     #[test]
//     fn test_kmers() {
//         let fasta_str = r">Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1;
// AAACCCCGT
// >Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1;
// TAACCCCGG
// >Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus2,s:Species3;
// TTTAAAACC
// >Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus2,s:Species3;
// TTTAAAACA
// >Badabing|Badabum;tax=p:Phylum1,c:Class2,o:Order2,f:Family2,g:Genus3,s:Species4;
// AAACCCCGG";
//         let LookupTables { k_mer_map, .. } = parse_reference_fasta_str(fasta_str).unwrap();
//         // for (k, v) in k_mer_map.iter().enumerate() {
//         //     if !v.is_empty() {
//         //         println!("{k:b}:\n {v:?}");
//         //     }
//         // }
//         assert_eq!(
//             k_mer_map[0b101010110_usize].iter().sorted().collect_vec(),
//             &[&0, &4]
//         );
//         assert_eq!(
//             k_mer_map[0b10101011010_usize].iter().sorted().collect_vec(),
//             &[&1, &4]
//         );
//         assert_eq!(
//             k_mer_map[0b10101011011_usize].iter().sorted().collect_vec(),
//             &[&0]
//         );
//         assert_eq!(
//             k_mer_map[0b1100000101010110_usize]
//                 .iter()
//                 .sorted()
//                 .collect_vec(),
//             &[&1]
//         );
//         assert_eq!(
//             k_mer_map[0b1111000000000101_usize]
//                 .iter()
//                 .sorted()
//                 .collect_vec(),
//             &[&2]
//         );
//         assert_eq!(
//             k_mer_map[0b1111000000000101_usize]
//                 .iter()
//                 .sorted()
//                 .collect_vec(),
//             &[&2]
//         );
//         assert_eq!(
//             k_mer_map[0b1111110000000001_usize]
//                 .iter()
//                 .sorted()
//                 .collect_vec(),
//             &[&2, &3]
//         );
//     }
// }
