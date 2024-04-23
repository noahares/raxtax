use crate::parser::LookupTables;
use crate::{prob, utils};
use indicatif::{ParallelProgressIterator, ProgressStyle};
use itertools::Itertools;
use log::{debug, log_enabled, warn, Level};
use logging_timer::{time, timer};
use rayon::prelude::*;

#[time("info")]
pub fn sintax<'a, 'b>(
    (query_labels, query_sequences): &'b (Vec<String>, Vec<Vec<u8>>),
    lookup_table: &'a LookupTables,
    num_k_mers: usize,
    max_target_seqs: usize,
    skip_exact_matches: bool,
) -> Vec<(&'b String, Option<Vec<(&'a String, Vec<f64>)>>)> {
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
            if !skip_exact_matches {
                if let Some(label_idxs) = lookup_table.sequences.get(query_sequence) {
                    debug!("Exact sequence match for query {query_label}");
                    if !label_idxs.iter().map(|&idx| lookup_table.labels[idx].split('|').take(5).join("|")).all_equal() {
                        warn!("Exact matches for {query_label} differ above the species level! Confidence values will be wrong!");
                        if log_enabled!(Level::Warn) {
                        eprintln!("\x1b[33m[WARN ]\x1b[0m Exact matches for {query_label} differ above the species level! Confidence values will be wrong!");
                        }
                    }
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
                                .take(max_target_seqs)
                                .collect_vec(),
                        ),
                    );
                }
            }
            let mut intersect_buffer: Vec<usize> = vec![0; lookup_table.labels.len()];
            let mut hit_buffer: Vec<f64> = vec![0.0; lookup_table.level_hierarchy_maps[5].len()];
            let _tmr = timer!(Level::Debug; "Query Time");
            let k_mers = utils::sequence_to_kmers(query_sequence);
            k_mers
                .iter()
                .for_each(|query_kmer| {
                    lookup_table.k_mer_map[*query_kmer as usize]
                        .iter()
                        .for_each(|species_id| {
                            unsafe { *intersect_buffer.get_unchecked_mut(*species_id) += 1 };
                        });
                });
            if skip_exact_matches {
                if let Some(exact_matches) = lookup_table.sequences.get(query_sequence) {
                    exact_matches.iter().for_each(|&id| unsafe { *intersect_buffer.get_unchecked_mut(id) = 0 });
                }
            }
            let higest_hit_probs = prob::highest_hit_prob_per_reference(k_mers.len(), num_k_mers, &intersect_buffer);
            let probs_sum: f64 = higest_hit_probs.iter().sum();
            if probs_sum > 0.0 {
                higest_hit_probs.iter().enumerate().for_each(|(idx, &v)| {
                    unsafe {
                        *hit_buffer.get_unchecked_mut(
                            *lookup_table.sequence_species_map.get_unchecked(idx),
                        ) += v / probs_sum
                    };
                });
            }
            (
                i,
                query_label,
                utils::accumulate_results(lookup_table, &hit_buffer, max_target_seqs),
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

    use crate::parser::{parse_query_fasta_str, parse_reference_fasta_str};

    use super::sintax;

    #[test]
    fn test_sintax() {
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
        let result = sintax(&query_data, &lookup_table, 32, 3, false);
        assert_eq!(
            vec![(
                &query_data.0[0],
                Some(vec![(&lookup_table.labels[0], vec![0.47; 6])])
            )],
            result
        );
    }
}
