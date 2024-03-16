use crate::utils;
use crate::{io::Args, parser::LookupTables};
use itertools::Itertools;
use rand::seq::SliceRandom;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

pub fn sintax(
    query_data: &(Vec<String>, Vec<Vec<u8>>),
    lookup_table: &LookupTables,
    args: &Args,
) -> Vec<String> {
    let (query_labels, query_sequences) = query_data;
    query_labels
        .par_iter()
        .zip_eq(query_sequences.par_iter())
        .enumerate()
        .map(|(i, (query_label, query_sequence))| {
            // let mut purged_runs = 0_usize;
            let mut buffer: Vec<usize> = vec![0; lookup_table.labels.len()];
            let mut hit_buffer: Vec<f64> = vec![0.0; lookup_table.labels.len()];
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(args.seed);
            let k_mers = utils::sequence_to_kmers(query_sequence);
            for _ in 0..args.num_rounds {
                buffer.fill(0);
                let selected_kmers = k_mers
                    .choose_multiple(&mut rng, args.num_k_mers)
                    .collect_vec();
                for query_kmer in selected_kmers {
                    for species_id in &lookup_table.k_mer_map[*query_kmer as usize] {
                        buffer[*species_id] += 1;
                    }
                }
                let relevant_hits = buffer
                    .iter()
                    .enumerate()
                    .filter(|(_, h)| *h >= &11_usize)
                    .max_set_by_key(|(_, &value)| value);
                let num_hits = relevant_hits.len();
                relevant_hits.into_iter().for_each(|(idx, _)| {
                    hit_buffer[idx] += 1.0 / (num_hits * args.num_rounds) as f64
                });
            }
            (
                i,
                utils::accumulate_results(lookup_table, &hit_buffer, args.num_results, query_label),
            )
        })
        .collect::<Vec<(usize, Vec<String>)>>()
        .into_iter()
        .sorted_by_key(|(i, _)| *i)
        .map(|(i, r)| r.join("\n"))
        .collect_vec()
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::parser::LookupTables;

    use crate::parser::{parse_query_fasta_str, parse_reference_fasta_str};

    #[test]
    fn test_sintax() {
        let fasta_str = r" >BOLD:AAP6467|SSBAE436-13|Canada|Arthropoda|Insecta|Diptera|Sciaridae|Claustropyga|Claustropyga_acanthostyla
TTTATCTTCTACATTATCTCACTCAGGGGCTTCAGTAGATCTATCTATTTTTTCTTTACATTTAGCAGGTATTTCATCAATTTTAGGAGCTGTAAATTTTATTTCTACTATTATTAATATACGAGCGCCAGGAATATCTTTTGATAAAATACCCTTATTTATTTGATCTGTATTAATTACAGCAATTTTATTATTATTATCATTA";
        let LookupTables { k_mer_map, .. } = parse_reference_fasta_str(fasta_str).unwrap();
        let query_str = r" >ESV_1;size=200394
TCTTTCATCTACTTTATCTCATTCAGGGGCTTCAGTAGATCTTTCTATTTTTTCCCTTCATTTAGCTGGAATTTCTTCAATTTTAGGGGCTGTAAATTTCATTTCAACTATTATTAATATACGGACACCAGGGATATCTTTTGATAAAATGTCTTTATTTATTTGATCGGTATTAATCACGGCCATTCTTTTGCTTTTATCATTA
";
        let (_, sequences) = parse_query_fasta_str(query_str).unwrap();
        // let k_mer_match = k_mer_map.keys
    }
}
