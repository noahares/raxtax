use std::{
    collections::HashSet,
    fs::File,
    io::{BufReader, Read, Write},
    path::PathBuf, arch::x86_64::{_mm_loadu_si128, __m128i, _mm_cmpestrm, _SIDD_BIT_MASK, _SIDD_UWORD_OPS, _SIDD_CMP_EQUAL_ANY, _mm_extract_epi32, _popcnt32, _mm_extract_epi16},
};

use flate2::read::GzDecoder;
use indicatif::{ParallelProgressIterator, ProgressStyle};
use itertools::Itertools;
use logging_timer::time;
use rayon::prelude::*;

use crate::parser::LookupTables;
use anyhow::{bail, Result};

pub const F64_OUTPUT_ACCURACY: u32 = 2;

pub fn sequence_to_kmers(sequence: &[u8]) -> Vec<u16> {
    let mut k_mer: u16 = 0;
    let mut k_mers = HashSet::new();
    sequence[0..8]
        .iter()
        .enumerate()
        .for_each(|(j, c)| k_mer |= (*c as u16) << (14 - j * 2));
    k_mers.insert(k_mer);
    sequence[8..].iter().for_each(|c| {
        k_mer = (k_mer << 2) | *c as u16;
        k_mers.insert(k_mer);
    });
    k_mers.into_iter().sorted().collect_vec()
}

fn sse_ordered_set_intersection_count(set_a: &[u16], set_b: &[u16]) -> usize {
    let mut idx_a = 0_isize;
    let mut idx_b = 0_isize;
    let mut count = 0;
    let last_bulk_index_a = (set_a.len() >> 3) << 3;
    let last_bulk_index_b = (set_b.len() >> 3) << 3;
    let remainder_a = (set_a.len() - last_bulk_index_a) as i32;
    let remainder_b = (set_b.len() - last_bulk_index_b) as i32;

    const IMM8: i32 = _SIDD_UWORD_OPS|_SIDD_CMP_EQUAL_ANY|_SIDD_BIT_MASK;
    let mut ptr_a = set_a.as_ptr();
    let mut ptr_b = set_b.as_ptr();
    while (idx_a as usize) < last_bulk_index_a && (idx_b as usize) < last_bulk_index_b {
        unsafe {
        let v_a = _mm_loadu_si128(ptr_a as *const __m128i);
        let v_b = _mm_loadu_si128(ptr_b as *const __m128i);
        let res = _mm_cmpestrm::<IMM8>(v_a, 8, v_b, 8);
        let r = _mm_extract_epi32(res, 0);
        count += _popcnt32(r) as usize;
        let a7 = _mm_extract_epi16::<7>(v_a);
        let b7 = _mm_extract_epi16::<7>(v_b);
        idx_a += (a7 <= b7) as isize * 8;
        idx_b += (a7 >= b7) as isize * 8;
        ptr_a = ptr_a.add((a7 <= b7) as usize * 8);
        ptr_b = ptr_b.add((a7 >= b7) as usize * 8);
        }
    }
    unsafe {
        let v_a = _mm_loadu_si128(ptr_a as *const __m128i);
        let v_b = _mm_loadu_si128(ptr_b as *const __m128i);
        let res = _mm_cmpestrm::<IMM8>(v_a, remainder_a, v_b, remainder_b);
        let r = _mm_extract_epi32(res, 0);
        count += _popcnt32(r) as usize;
    }
    count
}

pub fn create_intersection_distribution_parameters(database_kmers: &[Vec<u16>], query_kmers: &[Vec<u16>]) -> Vec<f64> {
    query_kmers.par_iter()
        .progress_with_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:80.cyan/blue} {pos:>7}/{len:7}[ETA:{eta}] {msg}",
            )
            .unwrap()
            .progress_chars("##-"),
        )
        .with_message("Calculating intersection parameters...")
        .flat_map_iter(|q| database_kmers.iter().map(|d| {
            sse_ordered_set_intersection_count(q, d) as f64 / q.len() as f64
        }))
    .collect()
}

pub fn check_convergence(
    hit_buffer: &[f64],
    last_hit_buffer: &[f64],
    iteration: usize,
    min_contribution: f64,
    early_stop_mse: f64,
) -> bool {
    let relevant_values = hit_buffer
        .iter()
        .zip_eq(last_hit_buffer.iter())
        .filter(|(&v, _)| v / (iteration + 1) as f64 >= min_contribution)
        .collect_vec();
    let num_values = relevant_values.len();
    if num_values == 0 {
        return true;
    }
    let mse = relevant_values.into_iter().fold(0.0, |acc, (&c, &l)| {
        acc + ((c / (iteration + 1) as f64) - (l / iteration as f64)).powi(2)
    }) / num_values as f64;
    mse < early_stop_mse
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

#[time("debug")]
pub fn accumulate_results<'a>(
    lookup_tables: &'a LookupTables,
    hit_buffer: &[f64],
    cutoff: usize,
) -> Option<Vec<(&'a String, Vec<f64>)>> {
    let rounding_factor = 10_u32.pow(F64_OUTPUT_ACCURACY) as f64;
    let mut result_lines: Option<Vec<(&'a String, Vec<f64>)>> = None;
    let mut phylum_values = vec![0.0; lookup_tables.level_hierarchy_maps[0].len()];
    let mut class_values = vec![0.0; lookup_tables.level_hierarchy_maps[1].len()];
    let mut order_values = vec![0.0; lookup_tables.level_hierarchy_maps[2].len()];
    let mut family_values = vec![0.0; lookup_tables.level_hierarchy_maps[3].len()];
    let mut genus_values = vec![0.0; lookup_tables.level_hierarchy_maps[4].len()];
    let mut species_values = vec![0.0; lookup_tables.level_hierarchy_maps[5].len()];
    for (a, phylum) in lookup_tables.level_hierarchy_maps[0].iter().enumerate() {
        for class in phylum {
            lookup_tables.level_hierarchy_maps[1][*class]
                .iter()
                .for_each(|order| {
                    lookup_tables.level_hierarchy_maps[2][*order]
                        .iter()
                        .for_each(|family| {
                            lookup_tables.level_hierarchy_maps[3][*family]
                                .iter()
                                .for_each(|genus| {
                                    lookup_tables.level_hierarchy_maps[4][*genus]
                                        .iter()
                                        .for_each(|species| {
                                            genus_values[*genus] += hit_buffer[*species];
                                            species_values[*species] =
                                                (hit_buffer[*species] * rounding_factor).round()
                                                    / rounding_factor;
                                        });
                                    family_values[*family] += genus_values[*genus];
                                    genus_values[*genus] = (genus_values[*genus] * rounding_factor)
                                        .round()
                                        / rounding_factor;
                                });
                            order_values[*order] += family_values[*family];
                            family_values[*family] = (family_values[*family] * rounding_factor)
                                .round()
                                / rounding_factor;
                        });
                    class_values[*class] += order_values[*order];
                    order_values[*order] =
                        (order_values[*order] * rounding_factor).round() / rounding_factor;
                });
            phylum_values[a] += class_values[*class];
            class_values[*class] =
                (class_values[*class] * rounding_factor).round() / rounding_factor;
        }
        phylum_values[a] = (phylum_values[a] * rounding_factor).round() / rounding_factor;
    }

    let mut out_count = 0_usize;
    for (a, _) in phylum_values
        .iter()
        .enumerate()
        .sorted_by(|a, b| b.1.partial_cmp(a.1).unwrap())
    {
        if phylum_values[a] == 0.0 {
            break;
        }
        for (_, b) in lookup_tables.level_hierarchy_maps[0][a]
            .iter()
            .enumerate()
            .sorted_by(|a, b| class_values[*b.1].partial_cmp(&class_values[*a.1]).unwrap())
        {
            if class_values[*b] == 0.0 {
                break;
            }
            for (_, c) in lookup_tables.level_hierarchy_maps[1][*b]
                .iter()
                .enumerate()
                .sorted_by(|a, b| order_values[*b.1].partial_cmp(&order_values[*a.1]).unwrap())
            {
                if order_values[*c] == 0.0 {
                    break;
                }
                for (_, d) in lookup_tables.level_hierarchy_maps[2][*c]
                    .iter()
                    .enumerate()
                    .sorted_by(|a, b| {
                        family_values[*b.1]
                            .partial_cmp(&family_values[*a.1])
                            .unwrap()
                    })
                {
                    if family_values[*d] == 0.0 {
                        break;
                    }
                    for (_, e) in lookup_tables.level_hierarchy_maps[3][*d]
                        .iter()
                        .enumerate()
                        .sorted_by(|a, b| {
                            genus_values[*b.1].partial_cmp(&genus_values[*a.1]).unwrap()
                        })
                    {
                        if genus_values[*e] == 0.0 {
                            break;
                        }
                        for species in
                            lookup_tables.level_hierarchy_maps[4][*e]
                                .iter()
                                .sorted_by(|a, b| {
                                    species_values[**b]
                                        .partial_cmp(&species_values[**a])
                                        .unwrap()
                                })
                        {
                            if species_values[*species] == 0.0 {
                                break;
                            }
                            let label = &lookup_tables.labels
                                [lookup_tables.level_hierarchy_maps[5][*species][0]];
                            let conf_values = vec![
                                phylum_values[a],
                                class_values[*b],
                                order_values[*c],
                                family_values[*d],
                                genus_values[*e],
                                species_values[*species],
                            ];
                            result_lines
                                .get_or_insert(Vec::new())
                                .push((label, conf_values));
                            out_count += 1;
                            if out_count == cutoff {
                                return result_lines;
                            }
                        }
                    }
                }
            }
        }
    }
    result_lines
}

pub fn output_results(
    results: &[(&String, Option<Vec<(&String, Vec<f64>)>>)],
    mut output: Box<dyn Write>,
    mut confidence_output: Box<dyn Write>,
    min_confidence: f64,
) -> Result<()> {
    let (output_lines, confidence_lines): (Vec<String>, Vec<String>) = results
        .iter()
        .map(|(query_label, confidence_vec)| match confidence_vec {
            Some(c) => {
                let normal_output = c
                    .iter()
                    .map(|(label, values)| {
                        format!(
                            "{}\t{}\t{}",
                            query_label,
                            label,
                            values
                                .iter()
                                .map(|v| format!("{1:.0$}", F64_OUTPUT_ACCURACY as usize, v))
                                .join("|")
                        )
                    })
                    .join("\n");
                let confidence_output = c
                    .iter()
                    .map(|(label, values)| {
                        let (filter_taxon_info, filter_values): (Vec<&str>, Vec<String>) = values
                            .iter()
                            .zip_eq(label.split('|'))
                            .filter(|(&v, _)| v >= min_confidence)
                            .map(|(v, l)| (l, format!("{1:.0$}", F64_OUTPUT_ACCURACY as usize, v)))
                            .unzip();
                        if filter_taxon_info.is_empty() {
                            format!("{}\tNA\tNA", query_label)
                        } else {
                            format!(
                                "{}\t{}\t{}",
                                query_label,
                                filter_taxon_info.join("|"),
                                filter_values.join("|")
                            )
                        }
                    })
                    .unique()
                    .collect_vec();
                // .join("\n");
                let mut confidence_output_filtered: Vec<String> =
                    vec![confidence_output[0].clone()];
                for i in 1..confidence_output.len() {
                    if !confidence_output[i - 1]
                        .split('\t')
                        .nth(1)
                        .unwrap()
                        .starts_with(confidence_output[i].split('\t').nth(1).unwrap())
                    {
                        confidence_output_filtered.push(confidence_output[i].clone());
                    }
                }
                (normal_output, confidence_output_filtered.join("\n"))
            }
            None => (
                format!("{}\tNA\tNA", query_label),
                format!("{}\tNA\tNA", query_label),
            ),
        })
        .unzip();
    writeln!(output, "{}", output_lines.into_iter().join("\n"))?;
    writeln!(
        confidence_output,
        "{}",
        confidence_lines.into_iter().join("\n")
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::{parser::parse_reference_fasta_str, utils::sse_ordered_set_intersection_count};

    use super::accumulate_results;

    #[test]
    fn test_accumulation() {
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
        let lookup_table = parse_reference_fasta_str(fasta_str).unwrap();
        let hit_buffer = [1.0 / 8.0, 2.0 / 8.0, 0.0, 2.0 / 8.0, 3.0 / 8.0];
        let results = accumulate_results(&lookup_table, &hit_buffer, 4);
        assert_eq!(
            results,
            Some(vec![
                (
                    &"Phylum1|Class1|Order1|Family1|Genus1|Species2".to_string(),
                    vec![0.63_f64, 0.38, 0.38, 0.38, 0.38, 0.25],
                ),
                (
                    &"Phylum1|Class1|Order1|Family1|Genus1|Species1".into(),
                    vec![0.63_f64, 0.38, 0.38, 0.38, 0.38, 0.13],
                ),
                (
                    &"Phylum1|Class2|Order2|Family3|Genus3|Species6".into(),
                    vec![0.63_f64, 0.25, 0.25, 0.25, 0.25, 0.25],
                ),
                (
                    &"Phylum2|Class3|Order3|Family4|Genus4|Species5".into(),
                    vec![0.38_f64, 0.38, 0.38, 0.38, 0.38, 0.38],
                ),
            ])
        );
    }

    #[test]
    fn test_fast_intersection() {
        let set_a = vec![0_u16, 1, 4, 5, 8, 12, 16, 17, 18, 29, 33, 43, 55, 61, 75, 85, 90];
        let set_b = vec![0_u16, 1, 2, 5, 9, 12, 13, 17, 18, 31, 34, 43, 55, 61, 75, 85, 88, 90, 100];
        assert_eq!(sse_ordered_set_intersection_count(&set_a, &set_b), 12);
    }
}
