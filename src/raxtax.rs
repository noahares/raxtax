use std::sync::atomic::AtomicBool;

use crate::lineage;
use crate::{prob, utils};
use indicatif::{ParallelProgressIterator, ProgressStyle};
use itertools::Itertools;
use log::{info, log_enabled, warn, Level};
use logging_timer::{time, timer};
use rayon::prelude::*;

#[time("info")]
pub fn raxtax<'a, 'b>(
    (query_labels, query_sequences): &'b (Vec<String>, Vec<Vec<u8>>),
    tree: &'a lineage::Tree,
    skip_exact_matches: bool,
) -> Vec<Vec<lineage::EvaluationResult<'a, 'b>>> {
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
                if let Some(label_idxs) = tree.sequences.get(query_sequence) {
                    info!("Exact sequence match for query {query_label}");
                    if !label_idxs.iter().map(|&idx| tree.lineages[idx].rsplit_once(',').unwrap().0).all_equal() {
                        warn!("Exact matches for {query_label} differ above the species level! Confidence values will be wrong!");
                        warnings.store(true, std::sync::atomic::Ordering::Relaxed);
                    }
                    return (
                        i,
                            label_idxs.iter()
                                .map(|&idx| {
                                    lineage::EvaluationResult{
                                        query_label,
                                        lineage: &tree.lineages[idx],
                                        confidence_values: tree.get_shared_exact_match(tree.lineages[idx].chars().filter(|c| *c == ',').count(), label_idxs.len()),
                                        local_signal: 1.0,
                                        global_signal: 1.0
                                    }
                                })
                                .collect_vec(),
                    );
                }
            }
            let mut intersect_buffer: Vec<usize> = vec![0; tree.num_tips];
            let _tmr = timer!(Level::Debug; "Query Time");
            let k_mers = utils::sequence_to_kmers(query_sequence);
            let num_trials = query_sequence.len() / 2;
            k_mers
                .iter()
                .for_each(|query_kmer| {
                    tree.k_mer_map[*query_kmer as usize]
                        .iter()
                        .for_each(|sequence_id| {
                            unsafe { *intersect_buffer.get_unchecked_mut(*sequence_id) += 1 };
                        });
                });
            if skip_exact_matches {
                if let Some(exact_matches) = tree.sequences.get(query_sequence) {
                    exact_matches.iter().for_each(|&id| unsafe { *intersect_buffer.get_unchecked_mut(id) = 0 });
                }
            }
            let higest_hit_probs = prob::highest_hit_prob_per_reference(k_mers.len(), num_trials, &intersect_buffer);
            (
                i,
                lineage::Lineage::new(query_label, tree, &higest_hit_probs).evaluate(),
            )
        })
        .collect::<Vec<(usize, Vec<lineage::EvaluationResult<'a, 'b>>)>>()
        .into_iter()
        .sorted_by_key(|(i, _)| *i)
        .map(|(_, v)| v)
        .collect_vec();

    if warnings.into_inner() && log_enabled!(Level::Warn) {
        eprintln!("\x1b[33m[WARN ]\x1b[0m Exact matches for some queries differ above the species level! Check the log file for more information!");
    }
    results
}
