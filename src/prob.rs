use ahash::{HashMap, HashMapExt};

use itertools::Itertools;
use logging_timer::time;
use statrs::function::factorial::ln_binomial;

#[time("debug")]
pub fn highest_hit_prob_per_reference(
    total_num_k_mers: u16,
    num_trials: usize,
    intersection_sizes: &[u16],
) -> Vec<f64> {
    let intersection_size_counts = {
        let mut counts: HashMap<u16, usize> = HashMap::new();
        intersection_sizes
            .iter()
            .for_each(|item| *counts.entry(*item).or_default() += 1);
        counts
    };
    let num_possible_kmer_sets = ln_binomial(
        total_num_k_mers as u64 + num_trials as u64 - 1,
        num_trials as u64,
    );
    let highest_hit_probs = if intersection_size_counts
        .iter()
        .any(|(&i, _)| i == total_num_k_mers)
    {
        intersection_size_counts
            .iter()
            .map(|(&n_intersections, _)| {
                (
                    n_intersections,
                    only_last_pmf(
                        total_num_k_mers as u64,
                        num_trials as u64,
                        n_intersections as u64,
                        num_possible_kmer_sets,
                    ),
                )
            })
            .collect::<HashMap<u16, f64>>()
    } else {
        let pmfs: Vec<(u16, Vec<f64>)> = iterative_pmfs_ln(
            total_num_k_mers as u64,
            num_trials as u64,
            &intersection_size_counts,
            num_possible_kmer_sets,
        );
        let cmfs: Vec<Vec<f64>> = pmfs
            .iter()
            .map(|(_, v)| {
                v.iter()
                    .scan(0.0, |sum, &pmf| {
                        if pmf != f64::NEG_INFINITY {
                            *sum += pmf.exp();
                        }
                        Some(sum.ln())
                    })
                    .collect_vec()
            })
            .collect_vec();
        let cmf_prod_components = (0..=num_trials)
            .map(|i| {
                intersection_size_counts
                    .iter()
                    .zip_eq(cmfs.iter())
                    .map(|((_, &count), cmf)| {
                        let x = unsafe { *cmf.get_unchecked(i) };
                        (count as f64) * x
                    })
                    .sum::<f64>()
            })
            .collect_vec();
        pmfs.into_iter()
            .zip_eq(cmfs.into_iter())
            .map(|((i, pmf), cmf)| {
                (
                    i,
                    itertools::izip!(pmf.into_iter(), cmf.into_iter(), cmf_prod_components.iter())
                        .map(|(p, c, &prod_components)| {
                            if c == f64::NEG_INFINITY || prod_components == f64::NEG_INFINITY {
                                0.0
                            } else {
                                (p + prod_components - c).exp()
                            }
                        })
                        .sum::<f64>(),
                )
            })
            .collect::<HashMap<u16, f64>>()
    };
    let highest_hit_probs = intersection_sizes
        .iter()
        .map(|&n_intersections| highest_hit_probs[&n_intersections])
        .collect_vec();

    let probs_sum: f64 = highest_hit_probs.iter().sum();
    assert!(probs_sum > 0.0);
    highest_hit_probs
        .into_iter()
        .map(|v| v / probs_sum)
        .collect_vec()
}

fn only_last_pmf(
    total_num_k_mers: u64,
    num_trials: u64,
    num_intersections: u64,
    num_possible_kmer_sets: f64,
) -> f64 {
    if num_intersections == total_num_k_mers {
        return 1.0;
    }
    if num_intersections == 0 {
        return 0.0;
    }
    let num_possible_matches = ln_binomial(num_intersections + num_trials - 1, num_trials);
    (num_possible_matches - num_possible_kmer_sets).exp()
}

fn iterative_pmfs_ln(
    total_num_k_mers: u64,
    num_trials: u64,
    intersection_sizes: &HashMap<u16, usize>,
    num_possible_kmer_sets: f64,
) -> Vec<(u16, Vec<f64>)> {
    intersection_sizes
        .iter()
        .map(|(&num_intersections, _)| {
            if num_intersections as u64 == total_num_k_mers {
                let mut res = vec![f64::NEG_INFINITY; num_trials as usize + 1];
                res[num_trials as usize] = 0.0;
                (num_intersections, res)
            } else if num_intersections == 0 {
                let mut res = vec![f64::NEG_INFINITY; num_trials as usize + 1];
                res[0] = 0.0;
                (num_intersections, res)
            } else {
                let num_possible_matches = (1..=num_trials).scan(0.0, |sum, i| {
                    *sum += ((num_intersections as u64 + i - 1) as f64 / i as f64).ln();
                    Some(*sum)
                });
                let impossible_init = ln_binomial(
                    total_num_k_mers - num_intersections as u64 + num_trials - 1,
                    num_trials,
                );
                let num_impossible_matches = (1..num_trials)
                    .scan(impossible_init, |sum, i| {
                        *sum -= ((total_num_k_mers - num_intersections as u64 + num_trials - i)
                            as f64
                            / (num_trials - i + 1) as f64)
                            .ln();
                        Some(*sum)
                    })
                    .chain([0.0]);
                (
                    num_intersections,
                    [impossible_init - num_possible_kmer_sets]
                        .into_iter()
                        .chain(
                            num_possible_matches
                                .zip_eq(num_impossible_matches)
                                .map(|(p, i)| p + i - num_possible_kmer_sets),
                        )
                        .collect_vec(),
                )
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use ahash::HashMap;
    use itertools::Itertools;
    use statrs::{assert_almost_eq, function::factorial::ln_binomial};

    use crate::prob::iterative_pmfs_ln;

    use super::highest_hit_prob_per_reference;

    fn pmf(
        total_num_k_mers: u64,
        i: u64,
        num_trials: u64,
        num_intersections: u64,
        num_possible_kmer_sets: f64,
    ) -> f64 {
        if num_intersections == total_num_k_mers {
            if i == num_trials {
                return 1.0;
            }
            return 0.0;
        }
        if num_intersections == 0 {
            if i == 0 {
                return 1.0;
            }
            return 0.0;
        }
        let num_possible_matches = ln_binomial(num_intersections + i - 1, i);
        let num_impossible_matches = ln_binomial(
            (total_num_k_mers - num_intersections) + (num_trials - i) - 1,
            num_trials - i,
        );
        (num_possible_matches + num_impossible_matches - num_possible_kmer_sets).exp()
    }
    #[test]
    fn test_pmf() {
        let num_possible_kmer_sets = ln_binomial(200 + 32 - 1, 32);
        let p = iterative_pmfs_ln(
            200,
            32,
            &HashMap::from_iter([(50, 4)]),
            num_possible_kmer_sets,
        );
        let p2 = (0..=32)
            .map(|i| pmf(200, i, 32, 50, num_possible_kmer_sets))
            .collect_vec();
        assert_almost_eq!(p[0].1.iter().map(|p| p.exp()).sum::<f64>(), 1.0, 1e-7);
        assert_almost_eq!(p2.iter().sum::<f64>(), 1.0, 1e-7);
        // assert_eq!(0, 1);
        p[0].1
            .iter()
            .zip(p2)
            .for_each(|(&a, b)| assert_almost_eq!(a.exp(), b, 1e-7));
    }

    #[test]
    fn test_hit_prob() {
        let probs = highest_hit_prob_per_reference(400, 200, &(0..=400).collect_vec());
        dbg!(&probs);
        assert_almost_eq!(probs.iter().sum::<f64>(), 1.0, 1e-7);
        assert!(probs.windows(2).all(|w| w[0] <= w[1]));
    }
}
