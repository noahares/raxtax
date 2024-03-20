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
pub fn sintax(
    query_data: &(Vec<String>, Vec<Vec<u8>>),
    lookup_table: &LookupTables,
    args: &Args,
) -> Vec<String> {
    let threshold = f64::ceil(args.num_k_mers as f64 * args.min_hit_fraction) as u8;
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
            if let Some(&label_idx) = lookup_table.sequences.get(query_sequence) {
                debug!("Exact sequence match for query {query_label}!");
                return (
                    i,
                    vec![format!(
                        "{}\t{}\t{}",
                        query_label,
                        lookup_table.labels[label_idx],
                        ["1.0"; 6].join("|")
                    )],
                );
            }
            // WARN: if number of possible hits can get above 255, this breaks! <noahares>
            let mut buffer: Vec<u8> = vec![0; lookup_table.labels.len()];
            let mut hit_buffer: Vec<f64> = vec![0.0; lookup_table.labels.len()];
            let mut rng = Xoroshiro128PlusPlus::seed_from_u64(args.seed);
            let _tmr = timer!(Level::Debug; "Query Time");
            let k_mers = utils::sequence_to_kmers(query_sequence);
            (0..args.num_iterations).for_each(|_| {
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
                        && relevant_hits.last().unwrap().0 < hit_buffer.len()
                );
                relevant_hits.into_iter().for_each(|(idx, _)| {
                    unsafe {
                        *hit_buffer.get_unchecked_mut(idx) +=
                            1.0 / (num_hits * args.num_iterations) as f64
                    };
                });
            });
            (
                i,
                utils::accumulate_results(
                    lookup_table,
                    &hit_buffer,
                    args.max_target_seqs,
                    query_label,
                ),
            )
        })
        .collect::<Vec<(usize, Vec<String>)>>()
        .into_iter()
        .sorted_by_key(|(i, _)| *i)
        .map(|(_, r)| r.join("\n"))
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
            database_path: None,
            database_output: None,
            num_iterations: 100,
            num_k_mers: 32,
            min_hit_fraction: 1.0 / 3.0,
            max_target_seqs: 3,
            threads: 1,
            seed: 42,
            output: None,
            verbosity: clap_verbosity_flag::Verbosity::default(),
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
        assert_eq!(vec![
    "ESV_1;size=200394 Arthropoda|Insecta|Diptera|Sciaridae|Claustropyga|Claustropyga_acanthostyla 0.4700|0.4700|0.4700|0.4700|0.4700|0.4700".to_string()], result);
    }
}