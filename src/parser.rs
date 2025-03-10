use anyhow::{bail, Context, Result};
use indicatif::{ProgressIterator, ProgressStyle};
use log::Level;
use logging_timer::{time, timer};
use regex::Regex;
use std::{io::Read, path::PathBuf};

use crate::{tree::Tree, utils};

fn map_dna_char(ch: char) -> u8 {
    let a: u8 = 0b0001;
    let c: u8 = 0b0010;
    let g: u8 = 0b0100;
    let t: u8 = 0b1000;
    match ch.to_ascii_uppercase() {
        'A' => a,
        'C' => c,
        'G' => g,
        'T' => t,
        'W' => a | t,
        'S' => c | g,
        'M' => a | c,
        'K' => g | t,
        'R' => a | g,
        'Y' => c | t,
        'B' => c | g | t,
        'D' => a | g | t,
        'H' => a | c | t,
        'V' => a | c | g,
        'N' => a | c | g | t,
        _ => panic!("Unexpected character: {ch}"),
    }
}

#[time("info")]
pub fn parse_reference_fasta_file(sequence_path: &PathBuf) -> Result<(bool, Tree)> {
    if let Ok(tree) = Tree::load_from_file(sequence_path) {
        return Ok((false, tree));
    }
    let mut fasta_str = String::new();
    let _ = utils::get_reader(sequence_path)?.read_to_string(&mut fasta_str);
    Ok((true, parse_reference_fasta_str(&fasta_str)?))
}

fn parse_reference_fasta_str(fasta_str: &str) -> Result<Tree> {
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
                            "Unexpected taxonomical annotation detected in label {label}"
                        ))?
                        .get(1)
                        .context(format!("No taxonomic string found in label {label}"))?
                        .as_str()
                        .to_owned();
                    labels.push(lineage);
                    if !current_sequence.is_empty() {
                        sequences.push(current_sequence.clone());
                        current_sequence = Vec::new();
                    }
                } else {
                    current_sequence.extend(line.chars().map(|c| -> u8 { map_dna_char(c) }));
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
    Tree::new(labels, sequences)
}

#[time("info")]
pub fn parse_query_fasta_file(sequence_path: &PathBuf) -> Result<Vec<(String, Vec<u8>)>> {
    let mut fasta_str = String::new();
    let _ = utils::get_reader(sequence_path)?.read_to_string(&mut fasta_str);
    parse_query_fasta_str(&fasta_str)
}

fn parse_query_fasta_str(fasta_str: &str) -> Result<Vec<(String, Vec<u8>)>> {
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
    let mut queries: Vec<(String, Vec<u8>)> = Vec::new();
    let mut current_query: (String, Vec<u8>) = (String::new(), Vec::new());

    // create label and sequence vectors
    for line in lines {
        if let Some(label) = line.strip_prefix('>') {
            if !current_query.1.is_empty() {
                queries.push(current_query.clone());
                current_query.1 = Vec::new();
            }
            current_query.0 = label.to_string();
        } else {
            current_query
                .1
                .extend(line.chars().map(|c| -> u8 { map_dna_char(c) }));
        }
    }
    queries.push(current_query);
    Ok(queries)
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::tree::Tree;

    use super::{parse_query_fasta_str, parse_reference_fasta_str};

    #[test]
    fn test_str_parser() {
        let fasta_str = r">Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1;
AAACCCTTTGGGA
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species2;
ATACGCTTTGGGA
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order4,f:Family5,g:Genus2,s:Species3;
ATCCGCTATGGGA
>Badabing|Badabum;tax=p:Phylum1,c:Class2,o:Order2,f:Family3,g:Genus3,s:Species6;
ATACGCTTTGCGT
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species2;
GTGCGCTATGCGA
>Badabing|Badabum;tax=p:Phylum2,c:Class3,o:Order3,f:Family4,g:Genus4,s:Species5;
ATACGCTTTGCGT";
        let tree = parse_reference_fasta_str(fasta_str).unwrap();
        for (k, v) in tree.k_mer_map.iter().enumerate() {
            if !v.is_empty() {
                println!("{k:b}:\n {v:?}");
            }
        }
        assert_eq!(
            tree.k_mer_map[0b1_0101_1111_1110_usize]
                .iter()
                .collect_vec(),
            &[&0]
        );
        assert_eq!(
            tree.k_mer_map[0b11_0001_1001_1111_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&1, &4, &5]
        );
        assert_eq!(
            tree.k_mer_map[0b110_0111_0011_1010_usize]
                .iter()
                .collect_vec(),
            &[&3]
        );
        assert_eq!(tree.num_tips, 6);
        assert_eq!(
            tree.lineages,
            vec![
                String::from("p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1"),
                "p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species2".into(),
                "p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species2".into(),
                "p:Phylum1,c:Class1,o:Order4,f:Family5,g:Genus2,s:Species3".into(),
                "p:Phylum1,c:Class2,o:Order2,f:Family3,g:Genus3,s:Species6".into(),
                "p:Phylum2,c:Class3,o:Order3,f:Family4,g:Genus4,s:Species5".into(),
            ]
        );
    }

    #[test]
    fn test_query_parser() {
        let fasta_str = r">label1
AAACCCTTTGGGA";
        let (_, sequence) = &parse_query_fasta_str(fasta_str).unwrap()[0];
        assert_eq!(sequence, &[1, 1, 1, 2, 2, 2, 8, 8, 8, 4, 4, 4, 1]);

        let fasta_str2 = r">label1
ACGTWSMKRYBDHVN";
        let (_, sequence) = &parse_query_fasta_str(fasta_str2).unwrap()[0];
        assert_eq!(
            sequence,
            &[1, 2, 4, 8, 9, 6, 3, 12, 5, 10, 14, 13, 11, 7, 15]
        );
    }

    #[test]
    fn test_kmers() {
        let fasta_str = r">Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1;
AAACCCCGT
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1;
TAACCCCGG
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus2,s:Species3;
TTTAAAACC
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus2,s:Species3;
TTTAAAACA
>Badabing|Badabum;tax=p:Phylum1,c:Class2,o:Order2,f:Family2,g:Genus3,s:Species4;
AAACCCCGG";
        let Tree { k_mer_map, .. } = parse_reference_fasta_str(fasta_str).unwrap();
        for (k, v) in k_mer_map.iter().enumerate() {
            if !v.is_empty() {
                println!("{k:b}:\n {v:?}");
            }
        }
        assert_eq!(
            k_mer_map[0b1_0101_0110_usize].iter().sorted().collect_vec(),
            &[&0, &4]
        );
        assert_eq!(
            k_mer_map[0b101_0101_1010_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&1, &4]
        );
        assert_eq!(
            k_mer_map[0b101_0101_1011_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&0]
        );
        assert_eq!(
            k_mer_map[0b1100_0001_0101_0110_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&1]
        );
        assert_eq!(
            k_mer_map[0b1111_0000_0000_0101_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&2]
        );
        assert_eq!(
            k_mer_map[0b1111_0000_0000_0101_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&2]
        );
        assert_eq!(
            k_mer_map[0b1111_1100_0000_0001_usize]
                .iter()
                .sorted()
                .collect_vec(),
            &[&2, &3]
        );
    }
}
