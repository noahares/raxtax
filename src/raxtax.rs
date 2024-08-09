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
    let warnings = std::sync::Mutex::new(false);
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
            let exact_matches = tree.sequences.get(query_sequence).map_or(Vec::new(), |m| m.to_owned());
            if !skip_exact_matches {
                let mut mtx = warnings.lock().unwrap();
                for id in &exact_matches {
                    info!("Exact sequence match for query {query_label}: {}", tree.lineages[*id]);
                }
                if !exact_matches.iter().map(|&idx| tree.lineages[idx].rsplit_once(',').unwrap().0).all_equal() {
                    warn!("Exact matches for {query_label} differ above the leafs of the lineage tree!");
                    *mtx = true;
                }
            }
            let _tmr = timer!(Level::Debug; "Query Time");
            let mut intersect_buffer: Vec<u16> = vec![0; tree.num_tips];
            let k_mers = utils::sequence_to_kmers(query_sequence);
            assert!(k_mers.len() <= u16::max_value() as usize);
            let num_trials = query_sequence.len() / 2;
            assert!(num_trials <= u16::max_value() as usize);
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
                exact_matches.iter().for_each(|&id| unsafe { *intersect_buffer.get_unchecked_mut(id) = 0 });
            }
            let highest_hit_probs = prob::highest_hit_prob_per_reference(k_mers.len() as u16, num_trials as u16, &intersect_buffer);
            let eval_res = lineage::Lineage::new(query_label, tree, &highest_hit_probs).evaluate();
            if let [idx] = exact_matches[..] {
                assert!(!eval_res.is_empty());
                let mut best_hit = eval_res[0].clone();
                best_hit.confidence_values = tree.get_shared_exact_match(tree.lineages[idx].chars().filter(|c| *c == ',').count(), 1);
                (i, vec![best_hit])
            } else {
                (i, eval_res)
            }
        })
        .collect::<Vec<(usize, Vec<lineage::EvaluationResult<'a, 'b>>)>>()
        .into_iter()
        .sorted_by_key(|(i, _)| *i)
        .map(|(_, v)| v)
        .collect_vec();

    if *warnings.lock().unwrap() && log_enabled!(Level::Warn) {
        eprintln!("\x1b[33m[WARN ]\x1b[0m Exact matches for some queries differ above the species level! Check the log file for more information!");
    }
    results
}
