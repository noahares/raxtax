use ahash::HashMap;

use itertools::Itertools;
use logging_timer::time;
use statrs::function::factorial::ln_binomial;

#[time("debug")]
pub fn highest_hit_prob_per_reference(
    total_num_k_mers: u16,
    num_trials: usize,
    intersection_sizes: &[u16],
) -> Vec<f64> {
    let intersection_size_counts = intersection_sizes.iter().counts().into_iter().collect_vec();
    let highest_hit_probs = if intersection_size_counts
        .iter()
        .any(|(&i, _)| i == total_num_k_mers)
    {
        intersection_size_counts
            .iter()
            .map(|(&n_intersections, _)| {
                (
                    n_intersections,
                    pmf(
                        total_num_k_mers as u64,
                        num_trials as u64,
                        num_trials as u64,
                        n_intersections as u64,
                    ),
                )
            })
            .collect::<HashMap<u16, f64>>()
    } else {
        let pmfs: Vec<(u16, Vec<f64>)> = intersection_size_counts
            .iter()
            .map(|(&n_intersections, _)| {
                (
                    n_intersections,
                    (0..=num_trials)
                        .map(|i| {
                            pmf(
                                total_num_k_mers as u64,
                                i as u64,
                                num_trials as u64,
                                n_intersections as u64,
                            )
                        })
                        .collect_vec(),
                )
            })
            .collect_vec();
        let cmfs: Vec<Vec<f64>> = pmfs
            .iter()
            .map(|(_, v)| {
                v.iter()
                    .scan(0.0, |sum, pmf| {
                        *sum += pmf;
                        Some(*sum)
                    })
                    .collect_vec()
            })
            .collect_vec();
        let cmf_prod_components = (0..=num_trials)
            .map(|i| {
                intersection_size_counts
                    .iter()
                    .zip_eq(cmfs.iter())
                    .map(|(&(_, count), cmf)| {
                        let x = unsafe { *cmf.get_unchecked(i) };
                        (count as f64) * x.ln()
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
                            if c == 0.0 || prod_components == f64::NEG_INFINITY {
                                0.0
                            } else {
                                (p.ln() + prod_components - c.ln()).exp()
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

pub fn pmf(total_num_k_mers: u64, i: u64, num_trials: u64, num_intersections: u64) -> f64 {
    if num_intersections == total_num_k_mers {
        if i == num_trials {
            return 1.0;
        } else {
            return 0.0;
        }
    }
    if num_intersections == 0 {
        if i == 0 {
            return 1.0;
        } else {
            return 0.0;
        }
    }
    let num_possible_matches = ln_binomial(num_intersections + i - 1, i);
    let num_impossible_matches = ln_binomial(
        (total_num_k_mers - num_intersections) + (num_trials - i) - 1,
        num_trials - i,
    );
    let num_possible_kmer_sets = ln_binomial(total_num_k_mers + num_trials - 1, num_trials);
    (num_possible_matches + num_impossible_matches - num_possible_kmer_sets).exp()
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use statrs::assert_almost_eq;

    use super::{highest_hit_prob_per_reference, pmf};

    #[test]
    fn test_pmf() {
        let p = pmf(200, 31, 32, 199);
        dbg!(p);
        assert_almost_eq!(p, 0.119857, 1e-7);
    }

    #[test]
    fn test_hit_prob() {
        let probs = highest_hit_prob_per_reference(400, 200, &(0..=400).collect_vec());
        dbg!(&probs);
        assert_almost_eq!(probs.iter().sum::<f64>(), 1.0, 1e-7);
    }
}
