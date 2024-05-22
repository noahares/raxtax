use std::sync::atomic::AtomicBool;

use crate::parser::LookupTables;
use crate::{prob, utils};
use indicatif::{ParallelProgressIterator, ProgressStyle};
use itertools::Itertools;
use log::{info, log_enabled, warn, Level};
use logging_timer::{time, timer};
use rayon::prelude::*;

#[time("info")]
pub fn raxtax<'a, 'b>(
    (query_labels, query_sequences): &'b (Vec<String>, Vec<Vec<u8>>),
    lookup_table: &'a LookupTables,
    skip_exact_matches: bool,
) -> Vec<(&'b String, Option<Vec<(&'a String, Vec<f64>)>>)> {
    let warnings = AtomicBool::new(false);
    let results = query_labels
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
                    info!("Exact sequence match for query {query_label}");
                    if !label_idxs.iter().map(|&idx| lookup_table.labels[idx].split('|').take(5).join("|")).all_equal() {
                        warn!("Exact matches for {query_label} differ above the species level! Confidence values will be wrong!");
                        warnings.store(true, std::sync::atomic::Ordering::Relaxed);
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
                                .collect_vec(),
                        ),
                    );
                }
            }
            let mut intersect_buffer: Vec<usize> = vec![0; lookup_table.labels.len()];
            let mut hit_buffer: Vec<f64> = vec![0.0; lookup_table.level_hierarchy_maps[5].len()];
            let _tmr = timer!(Level::Debug; "Query Time");
            let k_mers = utils::sequence_to_kmers(query_sequence);
            let num_trials = query_sequence.len() / 2;
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
            let higest_hit_probs = prob::highest_hit_prob_per_reference(k_mers.len(), num_trials, &intersect_buffer);
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
                utils::accumulate_results(lookup_table, &hit_buffer),
            )
        })
        .collect::<Vec<(usize, &String, Option<Vec<(&String, Vec<f64>)>>)>>()
        .into_iter()
        .sorted_by_key(|(i, _, _)| *i)
        .map(|(_, q, v)| (q, v))
        .collect_vec();

    if warnings.into_inner() && log_enabled!(Level::Warn) {
        eprintln!("\x1b[33m[WARN ]\x1b[0m Exact matches for some queries differ above the species level! Check the log file for more information!");
    }
    results
}
