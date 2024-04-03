use crate::utils;
use crate::{io::Args, parser::LookupTables};
use indicatif::{ParallelProgressIterator, ProgressStyle};
use itertools::Itertools;
use log::{debug, Level};
use logging_timer::{time, timer};
use rand::seq::SliceRandom;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
use rayon::prelude::*;

#[time("info")]
pub fn sintax<'a, 'b>(
    query_data: &'b (Vec<String>, Vec<Vec<u8>>),
    lookup_table: &'a LookupTables,
    args: &Args,
) -> Vec<(&'b String, Option<Vec<(&'a String, Vec<f64>)>>)> {
    let threshold = f64::ceil(args.num_k_mers as f64 * args.min_hit_fraction) as u8;
    let mse_discard_treshold = 1.0 / 10_u32.pow(utils::F64_OUTPUT_ACCURACY) as f64;
    let min_iterations: usize =
        (args.min_iterations * args.num_iterations as f64).max(1.0) as usize;
    let (query_labels, query_sequences) = query_data;
    query_labels
        .par_iter()
        .zip_eq(query_sequences.par_iter())
        .enumerate()
        .progress_with_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:80.cyan/blue} {pos:>7}/{len:7}[ETA:{eta}] {msg}",
            )
            .unwrap()
            .progress_chars("##-"),
        )
        .with_message("Running Queries...")
        .map(|(i, (query_label, query_sequence))| {
            if let Some(label_idxs) = lookup_table.sequences.get(query_sequence) {
                debug!("Exact sequence match for query {query_label}!");
                return (
                    i,
                    query_label,
                    Some(
                        label_idxs
                            .iter()
                            .map(|&idx| {
                                (
                                    &lookup_table.labels[idx],
                                    vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0 / label_idxs.len() as f64],
                                )
                            })
                            .take(args.max_target_seqs)
                            .collect_vec(),
                    ),
                );
            }
            // WARN: if number of possible hits can get above 255, this breaks! <noahares>
            let mut buffer: Vec<u8> = vec![0; lookup_table.labels.len()];
            let mut hit_buffer: Vec<f64> = vec![0.0; lookup_table.level_hierarchy_maps[5].len()];
            let mut rng = Xoroshiro128PlusPlus::seed_from_u64(args.seed);
            let _tmr = timer!(Level::Debug; "Query Time");
            let k_mers = utils::sequence_to_kmers(query_sequence);
            let mut last_hit_buffer: Option<Vec<f64>> = match args.no_early_stopping {
                true => None,
                false => Some(vec![0.0; hit_buffer.len()]),
            };
            let mut num_completed_iterations = args.num_iterations;
            for j in 0..args.num_iterations {
                buffer.fill(0);
                k_mers
                    .choose_multiple(&mut rng, args.num_k_mers)
                    .for_each(|query_kmer| {
                        lookup_table.k_mer_map[*query_kmer as usize]
                            .iter()
                            .for_each(|species_id| {
                                unsafe { *buffer.get_unchecked_mut(*species_id) += 1 };
                            });
                    });
                let relevant_hits = buffer
                    .iter()
                    .enumerate()
                    .filter(|(_, h)| **h >= threshold)
                    .max_set_by_key(|(_, &value)| value);
                let num_hits = relevant_hits.len();
                debug_assert!(
                    relevant_hits.iter().tuple_windows().all(|(a, b)| a.0 < b.0)
                        && relevant_hits.last().unwrap_or(&(0, &0)).0 < hit_buffer.len()
                );
                relevant_hits.into_iter().for_each(|(idx, _)| {
                    unsafe {
                        *hit_buffer.get_unchecked_mut(
                            *lookup_table.sequence_species_map.get_unchecked(idx),
                        ) += 1.0 / num_hits as f64
                    };
                });
                if !args.no_early_stopping && j + 1 >= min_iterations {
                    let mse = last_hit_buffer
                        .unwrap()
                        .iter()
                        .zip_eq(hit_buffer.iter())
                        .filter(|(_, c)| **c / (j + 1) as f64 >= mse_discard_treshold)
                        .fold(0.0, |acc, (&a, b)| {
                            acc + (a / j as f64 - b / (j + 1) as f64).powi(2)
                        })
                        / hit_buffer
                            .iter()
                            .filter(|&&v| v / (j + 1) as f64 >= mse_discard_treshold)
                            .count() as f64;
                    if mse < args.early_stop_mse {
                        num_completed_iterations = j + 1;
                        debug!("{query_label} Stopped after {num_completed_iterations} iterations");
                        break;
                    }
                    last_hit_buffer = Some(hit_buffer.clone());
                }
            }
            if num_completed_iterations == args.num_iterations {
                debug!("{query_label} Stopped after {} iterations", args.num_iterations);
            }
            hit_buffer
                .iter_mut()
                .for_each(|v| *v /= num_completed_iterations as f64);
            (
                i,
                query_label,
                utils::accumulate_results(lookup_table, &hit_buffer, args.max_target_seqs),
            )
        })
        .collect::<Vec<(usize, &String, Option<Vec<(&String, Vec<f64>)>>)>>()
        .into_iter()
        .sorted_by_key(|(i, _, _)| *i)
        .map(|(_, q, v)| (q, v))
        .collect_vec()
}

#[cfg(test)]
mod tests {

    use crate::{
        io::Args,
        parser::{parse_query_fasta_str, parse_reference_fasta_str},
    };

    use super::sintax;

    #[test]
    fn test_sintax() {
        let args = Args {
            database_path: "".into(),
            query_file: "".into(),
            num_iterations: 100,
            num_k_mers: 32,
            min_hit_fraction: 1.0 / 3.0,
            max_target_seqs: 3,
            threads: 1,
            seed: 42,
            output: None,
            verbosity: clap_verbosity_flag::Verbosity::default(),
            min_confidence: 0.0,
            confidence_output: None,
        };
        let fasta_str = r">BOLD:AAP6467|SSBAE436-13|Canada|Arthropoda|Insecta|Diptera|Sciaridae|Claustropyga|Claustropyga_acanthostyla
TTTATCTTCTACATTATCTCACTCAGGGGCTTCAGTAGATCTATCTATTTTTTCTTTACATTTAGCAGGTATTTCATCAATTTTAGGAGCTGTAAATTTTATTTCTACTATTATTAATATACGAGCGCCAGGAATATCTTTTGATAAAATACCCTTATTTATTTGATCTGTATTAATTACAGCAATTTTATTATTATTATCATTA
>BOLD:AAP6467|GMOXC9016-15|Canada|Arthropoda|Insecta|Diptera|Sciaridae|Claustropyga|Claustropyga_acanthostyla
TTTATCTTCTACATTATCTCATTCAGGAGCTTCAGTAGATCTATCTATTTTTTCTTTACATTTAGCAGGTATTTCATCAATTTTAGGGGCTGTAAATTTTATTTCTACTATTATTAATATACGAGCGCCAGGAATATCTTTTGATAAAATACCCTTATTTATTTGATCTGTATTAATTACAGCAATTTTATTATTATTATCATTA
>BOLD:AAP6467|GMOAB1920-21|Canada|Arthropoda|Insecta|Diptera|Sciaridae|Claustropyga|Claustropyga_acanthostyla
TTTATCTTCTACATTATCTCATTCAGGAGCTTCAGTAGATCTATCTATTTTTTCTTTACATTTAGCAGGTATTTCATCAATTTTAGGGGCTGTAAATTTTATTTCTACTATTATTAATATGCGAGCGCCAGGAATATCTTTTGATAAAATACCCTTATTTATTTGATCTGTATTAATTACAGCAATTTTATTATTATTATCATTA
>BOLD:AAP6467|SSJAD1807-13|Canada|Arthropoda|Insecta|Diptera|Sciaridae|Claustropyga|Claustropyga_acanthostyla
TTTATCTTCTACATTATCTCACTCAGGAGCTTCAGTAGATCTATCTATTTTTTCTTTACATTTAGCAGGTATTTCATCAATTTTAGGGGCTGTAAACTTTATTTCTACTATTATTAATATACGAGCACCAGGAATATCTTTTGATAAAATACCCTTATTTATTTGATCTGTATTAATTACAGCAATTTTATTATTATTATCATTA
>BOLD:AAP6467|CNBFF377-15|Canada|Arthropoda|Insecta|Diptera|Sciaridae|Claustropyga|Claustropyga_acanthostyla
TTTATCTTCTACATTATCTCACTCAGGAGCTTCAGTAGACCTATCTATTTTTTCTTTACATTTAGCCGGTATTTCATCAATTTTAGGAGCTGTAAATTTTATTTCTACTATTATTAATATACGAGCACCAGGAATATCTTTTGATAAAATACCTTTATTTATTTGATCTGTATTAATTACAGCAATTTTATTATTATTATCATTA";
        let lookup_table = parse_reference_fasta_str(fasta_str).unwrap();
        let query_str = r" >ESV_1;size=200394
TCTTTCATCTACTTTATCTCATTCAGGGGCTTCAGTAGATCTTTCTATTTTTTCCCTTCATTTAGCTGGAATTTCTTCAATTTTAGGGGCTGTAAATTTCATTTCAACTATTATTAATATACGGACACCAGGGATATCTTTTGATAAAATGTCTTTATTTATTTGATCGGTATTAATCACGGCCATTCTTTTGCTTTTATCATTA
";
        let query_data = parse_query_fasta_str(query_str).unwrap();
        let result = sintax(&query_data, &lookup_table, &args);
        assert_eq!(
            vec![(
                &query_data.0[0],
                Some(vec![(&lookup_table.labels[0], vec![0.47; 6])])
            )],
            result
        );
    }
}
