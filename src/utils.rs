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
    let output_lines = results.iter().flat_map(|eval_results| {
        eval_results
            .iter()
            .map(|er| er.get_output_string())
            .collect_vec()
    });
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
    let output_lines = results
        .iter()
        .zip_eq(sequences)
        .flat_map(|(eval_results, sequence)| {
            eval_results
                .iter()
                .map(|er| er.get_tsv_string(&sequence))
                .collect_vec()
        });
    writeln!(output, "{}", output_lines.into_iter().join("\n"))?;
    Ok(())
}

pub fn euclidean_distance_l1(a: &[f64], b: &[f64]) -> f64 {
    assert!(a.len() == b.len());
    if a.is_empty() {
        return 0.0;
    };
    let a_sum = a.iter().sum::<f64>();
    let b_sum = b.iter().sum::<f64>();
    assert!(a_sum > 0.0);
    assert!(b_sum > 0.0);
    a.iter()
        .zip(b)
        .map(|(x, y)| (x / a_sum - y / b_sum).powi(2))
        .sum::<f64>()
        .sqrt()
}

pub fn euclidean_norm<I, T>(v: I) -> f64
where
    I: IntoIterator<Item = T>,
    T: std::borrow::Borrow<f64>,
{
    v.into_iter()
        .map(|x| x.borrow() * x.borrow())
        .sum::<f64>()
        .sqrt()
}

pub fn cosine_similarity(vec_a: &[f64], vec_b: &[f64]) -> f64 {
    let norm_a = euclidean_norm(vec_a.iter());
    let norm_b = euclidean_norm(vec_b.iter());
    assert!(norm_a > 0.0);
    assert!(norm_b > 0.0);
    vec_a
        .iter()
        .zip(vec_b.iter())
        .map(|(a, b)| a * b)
        .sum::<f64>()
        / (norm_a * norm_b)
}

pub fn report_error(e: anyhow::Error, message: impl std::fmt::Display) {
    let prefix = "\x1b[31m[ERROR]\x1b[0m";
    log::error!("{}: {}", message, e);
    if log::log_enabled!(log::Level::Error) {
        eprintln!("{} {}: {}", prefix, message, e);
    }
}

#[cfg(test)]
mod tests {
    use statrs::assert_almost_eq;

    use crate::utils::{cosine_similarity, euclidean_distance_l1, euclidean_norm};

    #[test]
    fn test_euclidean_norm() {
        let v = [1.0, 2.0, 3.0, 4.0];
        assert_almost_eq!(euclidean_norm(v.iter()), 30_f64.sqrt(), 1e-7);
        let w = [0.5, 0.5, 0.25, 0.2];
        assert_almost_eq!(euclidean_norm(w.iter()), 0.6025_f64.sqrt(), 1e-7);
    }

    #[test]
    fn test_euclidean_distance() {
        let v = [1.0, 0.0, 0.0];
        let w = [0.0, 1.0, 0.0];
        assert_almost_eq!(euclidean_distance_l1(&v, &w), 2_f64.sqrt(), 1e-7);
        let x = [0.5, 0.1, 0.1];
        let y = [1.0, 1.0, 0.5];
        assert_almost_eq!(euclidean_distance_l1(&x, &y), 0.4100771455544949, 1e-7);
    }

    #[test]
    fn test_cosine_similarity() {
        let v = [1.0, 0.0, 0.0];
        let w = [0.0, 1.0, 0.0];
        assert_almost_eq!(cosine_similarity(&v, &w), 0.0, 1e-7);
        let x = [0.5, 0.5];
        let y = [0.5, 0.5];
        assert_almost_eq!(cosine_similarity(&x, &y), 1.0, 1e-7);
    }
}
