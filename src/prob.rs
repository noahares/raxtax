use std::collections::HashMap;

use itertools::Itertools;
use statrs::function::factorial::binomial;

// pub struct QueryOracle {
//     pmfs: HashMap<usize, Vec<f64>>,
//     cmfs: HashMap<usize, Vec<f64>>,
//     hit_probs: Vec<f64>
// }

// impl QueryOracle {
pub fn highest_hit_prob_per_reference(
    total_num_k_mers: u64,
    num_trials: u64,
    intersection_sizes: &[u32],
) -> Vec<f64> {
    let intersection_size_counts = intersection_sizes.iter().counts();
    let pmfs: HashMap<usize, Vec<f64>> = intersection_size_counts
        .keys()
        .map(|&&n_intersections| {
            (
                n_intersections as usize,
                (0..=num_trials)
                    .map(|i| pmf(total_num_k_mers, i, num_trials, n_intersections as u64))
                    .collect_vec(),
            )
        })
        .collect();
    let cmfs: HashMap<usize, Vec<f64>> = pmfs
        .iter()
        .map(|(&i, v)| {
            (
                i,
                v.iter()
                    .scan(0.0, |sum, pmf| {
                        *sum += pmf;
                        Some(*sum)
                    })
                    .collect_vec(),
            )
        })
        .collect();
    let cmf_prod_components = (0..=num_trials)
        .map(|i| {
            (0..=i)
                .map(|j| {
                    intersection_size_counts
                        .iter()
                        .map(|(&&size, &count)| {
                            let x = cmfs[&(size as usize)][j as usize];
                            if x < f64::EPSILON {
                                1.0
                            } else {
                                cmfs[&(size as usize)][j as usize].powi(count as i32)
                            }
                        })
                        .product::<f64>()
                })
                .collect_vec()
        })
        .collect_vec();
    let highest_hit_probs: HashMap<usize, f64> = pmfs
        .iter()
        .map(|(i, v)| {
            (
                *i,
                v.iter()
                    .enumerate()
                    .zip_eq(cmf_prod_components.iter())
                    .map(|((j, pdf), prod_components)| {
                        pdf * prod_components
                            .iter()
                            .map(|prod| {
                                let x = cmfs[i][j];
                                if x < f64::EPSILON {
                                    *prod
                                } else {
                                    prod / cmfs[i][j]
                                }
                            })
                            .sum::<f64>()
                    })
                    .sum(),
            )
        })
        .collect();
    intersection_sizes
        .iter()
        .map(|&n_intersections| highest_hit_probs[&(n_intersections as usize)])
        .collect_vec()
}

pub fn pmf(total_num_k_mers: u64, i: u64, num_trials: u64, num_intersections: u64) -> f64 {
    if num_intersections == total_num_k_mers && i == num_trials {
        return 1.0;
    }
    if num_intersections == 0 {
        if i == 0 {
            return 1.0;
        } else {
            return 0.0;
        }
    }
    // TODO: deal with cases where n or k are negative
    let num_possible_matches = binomial(num_intersections + i - 1, i);
    let num_impossible_matches = binomial(
        (total_num_k_mers - num_intersections) + (num_trials - i) - 1,
        num_trials - i,
    );
    let num_possible_kmer_sets = binomial(total_num_k_mers + num_trials - 1, num_trials);
    num_possible_matches * num_impossible_matches / num_possible_kmer_sets
}

pub fn cmf(pmfs: &[f64], max_matches: usize) -> f64 {
    pmfs[0..=max_matches].iter().sum()
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use statrs::assert_almost_eq;

    use crate::prob::cmf;

    use super::{highest_hit_prob_per_reference, pmf};

    #[test]
    fn test_pmf() {
        let p = pmf(200, 31, 32, 199);
        // dbg!(p);
        assert_almost_eq!(p, 0.119857, 1e-7);
    }

    #[test]
    fn test_cmf() {
        let pdfs = (0..=32).map(|i| pmf(200, i, 32, 199)).collect_vec();
        // dbg!(&pdfs);
        let c = cmf(&pdfs, 30);
        assert_almost_eq!(c, 0.5573349, 1e-7);
    }

    #[test]
    fn test_hit_prob() {
        let probs = highest_hit_prob_per_reference(200, 32, &(0..=200).collect_vec());
        dbg!(&probs);
        assert_eq!(probs.iter().sum::<f64>(), 1.0);
    }
}
