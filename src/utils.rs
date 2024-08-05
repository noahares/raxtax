use std::{
    collections::HashSet,
    fs::File,
    io::{BufReader, Read, Write},
    path::PathBuf,
};

use anyhow::{bail, Result};
use flate2::read::GzDecoder;
use itertools::Itertools;

use crate::lineage;

pub const F64_OUTPUT_ACCURACY: u32 = 2;

pub fn map_four_to_two_bit_repr(c: u8) -> Option<u16> {
    match c {
        0b0001 => Some(0b00),
        0b0010 => Some(0b01),
        0b0100 => Some(0b10),
        0b1000 => Some(0b11),
        _ => None,
    }
}

pub fn sequence_to_kmers(sequence: &[u8]) -> Vec<u16> {
    let mut k_mers = HashSet::new();
    sequence.windows(8).for_each(|vals| {
        if let Some(k_mer) = vals
            .iter()
            .enumerate()
            .map(|(j, v)| map_four_to_two_bit_repr(*v).map(|c| c << (14 - j * 2)))
            .fold_options(0_u16, |acc, c| acc | c)
        {
            k_mers.insert(k_mer);
        }
    });
    k_mers.into_iter().unique().sorted().collect_vec()
}

pub fn get_reader(path: &PathBuf) -> Result<Box<dyn Read>> {
    let file_type = match path.extension() {
        Some(ext) => match ext.to_str() {
            Some(ext_str) => ext_str.to_ascii_lowercase(),
            None => bail!("Extension could not be parsed!"),
        },
        None => "fasta".to_string(),
    };

    let file = File::open(path)?;

    match file_type.as_str() {
        "gz" | "gzip" => {
            let reader = Box::new(GzDecoder::new(file));
            Ok(Box::new(BufReader::new(reader)))
        }
        _ => Ok(Box::new(BufReader::new(file))),
    }
}

pub fn output_results(
    results: &[Vec<lineage::EvaluationResult<'_, '_>>],
    mut output: Box<dyn Write>,
) -> Result<()> {
    let output_lines: Vec<String> = results
        .iter()
        .map(|eval_results| {
            eval_results
                .iter()
                .map(|er| er.get_output_string())
                .join("\n")
        })
        .collect_vec();
    writeln!(output, "{}", output_lines.into_iter().join("\n"))?;
    Ok(())
}

pub fn decompress_sequences(sequences: &[Vec<u8>]) -> Vec<String> {
    sequences
        .iter()
        .map(|s| {
            s.iter()
                .map(|c| match c {
                    0b0001 => 'A',
                    0b0010 => 'C',
                    0b0100 => 'G',
                    0b1000 => 'T',
                    _ => '-',
                })
                .join("")
        })
        .collect_vec()
}

pub fn output_results_tsv(
    results: &[Vec<lineage::EvaluationResult<'_, '_>>],
    sequences: Vec<String>,
    mut output: Box<dyn Write>,
) -> Result<()> {
    let output_lines: Vec<String> = results
        .iter()
        .zip_eq(sequences)
        .map(|(eval_results, sequence)| {
            eval_results
                .iter()
                .map(|er| er.get_tsv_string(&sequence))
                .join("\n")
        })
        .collect_vec();
    writeln!(output, "{}", output_lines.into_iter().join("\n"))?;
    Ok(())
}
// #[cfg(test)]
// mod tests {
//     use crate::parser::parse_reference_fasta_str;
//
//     #[test]
//     fn test_accumulation() {
//         let fasta_str = r">Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1;
// AAACCCTTTGGGA
// >Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species2;
// ATACGCTTTGGGA
// >Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order4,f:Family5,g:Genus2,s:Species3;
// ATCCGCTATGGGA
// >Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order2,f:Family3,g:Genus3,s:Species6;
// ATACGCTTTGCGT
// >Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order3,f:Family4,g:Genus4,s:Species5;
// ATACGCTTTGCGT";
//         let lookup_table = parse_reference_fasta_str(fasta_str).unwrap();
//         let hit_buffer = [1.0 / 8.0, 2.0 / 8.0, 0.0, 2.0 / 8.0, 3.0 / 8.0];
//         let results = accumulate_results(&lookup_table, &hit_buffer);
//         assert_eq!(
//             results,
//             Some(vec![
//                 (
//                     &"Phylum1|Class1|Order1|Family1|Genus1|Species2".to_string(),
//                     vec![1.0_f64, 1.0, 0.38, 0.38, 0.38, 0.25],
//                 ),
//                 (
//                     &"Phylum1|Class1|Order1|Family1|Genus1|Species1".into(),
//                     vec![1.0_f64, 1.0, 0.38, 0.38, 0.38, 0.13],
//                 ),
//                 (
//                     &"Phylum1|Class1|Order2|Family3|Genus3|Species6".into(),
//                     vec![1.0_f64, 1.0, 0.38, 0.38, 0.38, 0.38],
//                 ),
//                 (
//                     &"Phylum1|Class1|Order3|Family4|Genus4|Species5".into(),
//                     vec![1.0_f64, 1.0, 0.25, 0.25, 0.25, 0.25],
//                 ),
//             ])
//         );
//     }
// }
