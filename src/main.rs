use clap::Parser;
use anyhow::Result;
use itertools::Itertools;
use log::Level;
use logging_timer::timer;
use make_sintax_great_again::io;
use make_sintax_great_again::parser;
use make_sintax_great_again::utils;
use rand::seq::SliceRandom;
use rand_xoshiro::Xoshiro256PlusPlus;
use rand_xoshiro::rand_core::SeedableRng;
use rayon::prelude::*;

fn main() -> Result<()> {
    let args = io::Args::parse();
    env_logger::Builder::new()
        .filter_level(args.verbosity.log_level_filter())
        .init();
    let _total_tmr = timer!(Level::Info; "Total Runtime");
    let lookup_table = if let Some(ref database_path) = args.database_path {
        let database_str =
            std::fs::read_to_string(database_path)?;
        serde_json::from_str(&database_str)
    } else {
        let lookup_table = parser::parse_reference_fasta_file(&args.sequence_file)?;
        if let Some(ref database_output_path) = args.database_output {
            let output = std::fs::File::create(database_output_path)?;
            serde_json::to_writer(output, &lookup_table)?;
        }
        Ok(lookup_table)
    }?;
    let (query_labels, query_sequences) = parser::parse_query_fasta_file(&args.query_file)?;
    let result: Vec<(usize, Vec<String>)> = query_labels.par_iter().zip_eq(query_sequences.par_iter()).enumerate().map(|(i, (query_label, query_sequence))| {
        // let mut purged_runs = 0_usize;
        let mut buffer:Vec<usize> = vec![0; lookup_table.labels.len()];
        let mut hit_buffer:Vec<f64> = vec![0.0; lookup_table.labels.len()];
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(args.seed);
        // log::info!("Processing {query_label}");
        let k_mers = utils::sequence_to_kmers(query_sequence);
        for _ in 0..args.num_rounds {
            let selected_kmers = k_mers.choose_multiple(&mut rng, args.num_k_mers).collect_vec();
            for query_kmer in selected_kmers {
                for sequence_id in &lookup_table.k_mer_map[*query_kmer as usize] {
                    buffer[*sequence_id] += 1;
                }
            }
            // let relevant_hits = buffer.iter().filter(|h| *h >= &11_usize).enumerate().max_set_by_key(|(_, &value)| value);
            let (idx, hits) = buffer.iter().enumerate().max_by_key(|(_, &value)| value).unwrap_or((0,&0));
            if *hits >= f64::ceil(args.num_k_mers as f64 * args.threshold) as usize {
                hit_buffer[idx] += 1.0 / args.num_rounds as f64;
            // } else {
            //     purged_runs += 1;
            }
            buffer.clear();
            buffer.resize(lookup_table.labels.len(), 0);
        }
        (i, utils::accumulate_results(&lookup_table, &hit_buffer, args.num_results, query_label))
    })
    .collect();
    let mut output = args.get_output()?;
    write!(output, "{}", result.into_iter().sorted_by_key(|(i, _)| *i).map(|(_, r)| r.join("\n")).join("\n"))?;

    Ok(())
}
