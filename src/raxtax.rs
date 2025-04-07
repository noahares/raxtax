use std::time::Duration;

use crate::lineage;
use crate::tree::Tree;
use crate::{prob, utils};
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use log::{info, log_enabled, warn, Level};
use logging_timer::{time, timer};
use rayon::prelude::*;

#[time("info")]
pub fn raxtax<'a, 'b>(
    queries: &'b Vec<(String, Vec<u8>)>,
    tree: &'a Tree,
    skip_exact_matches: bool,
    raw_confidence: bool,
    chunk_size: usize,
) -> Vec<Vec<lineage::EvaluationResult<'a, 'b>>> {
    let warnings = std::sync::Mutex::new(false);
    let empty_vec = Vec::new();
    let pb = ProgressBar::new(queries.len() as u64)
        .with_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:80.cyan/blue} {pos:>7}/{len:7}[ETA:{eta}] {msg}",
            )
            .unwrap()
            .progress_chars("##-"),
        )
        .with_message("Running Queries...");
    pb.enable_steady_tick(Duration::from_millis(100));
    let results = queries
        .par_chunks(chunk_size)
        .flat_map(|q| {
            let mut intersect_buffer: Vec<u16> = vec![0; tree.num_tips];
            q.iter().map(|(query_label, query_sequence)| {
                pb.inc(1);
                intersect_buffer.fill(0);
                let exact_matches = tree.sequences.get(query_sequence).unwrap_or(&empty_vec);
                if !skip_exact_matches {
                    // check for inconsistencies for exact matches
                    let mut mtx = warnings.lock().unwrap();
                    for id in exact_matches {
                        info!("Exact sequence match for query {query_label}: {}", tree.lineages[*id as usize]);
                    }
                    if !exact_matches.iter().map(|&idx| tree.lineages[idx as usize].rsplit_once(',').unwrap().0).all_equal() {
                        warn!("Exact matches for {query_label} differ above the leafs of the lineage tree!");
                        *mtx = true;
                    }
                }
                let tmr = timer!(Level::Debug; "K-mer Intersections");
                let k_mers = utils::sequence_to_kmers(query_sequence);
                assert!(u16::try_from(k_mers.len()).is_ok());
                let num_trials = k_mers.len() / 2;
                for query_kmer in &k_mers {
                        tree.k_mer_map[*query_kmer as usize]
                            .iter()
                            .for_each(|sequence_id| {
                                unsafe { *intersect_buffer.get_unchecked_mut(*sequence_id as usize) += 1 };
                            });
                    }
                if skip_exact_matches {
                    // look for the next best match
                    for &id in exact_matches { unsafe { *intersect_buffer.get_unchecked_mut(id as usize) = 0 } }
                }
                drop(tmr);
                let highest_hit_probs = prob::highest_hit_prob_per_reference(k_mers.len() as u16, num_trials, &intersect_buffer);
                let eval_res = lineage::Lineage::new(query_label, tree, highest_hit_probs).evaluate();
                assert!(!eval_res.is_empty());
                if !raw_confidence && !skip_exact_matches {
                    // Special case: if there is exactly 1 exact match, confidence is set to 1.0
                    if let [idx] = exact_matches[..] {
                        return vec![lineage::EvaluationResult {
                            query_label,
                            lineage: &tree.lineages[idx as usize],
                            confidence_values: vec![1.0; tree.lineages[idx as usize].chars().filter(|c| *c == ',').count() + 1],
                            local_signal: eval_res[0].local_signal,
                            global_signal: eval_res[0].global_signal
                        }];
                    }
                }
                eval_res
            }).collect_vec()
        })
        .collect::<Vec<Vec<lineage::EvaluationResult<'a, 'b>>>>();

    if *warnings.lock().unwrap() && log_enabled!(Level::Warn) {
        eprintln!("\x1b[33m[WARN ]\x1b[0m Exact matches for some queries differ above the species level! Check the log file for more information!");
    }
    results
}
