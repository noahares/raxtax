use anyhow::{bail, Result};
use clap::Parser;
use log::Level;
use logging_timer::timer;
use make_sintax_great_again::io;
use make_sintax_great_again::parser;
use make_sintax_great_again::sintax::sintax;
use make_sintax_great_again::utils;

fn main() -> Result<()> {
    let args = io::Args::parse();
    if args.num_k_mers >= 255 {
        bail!("Using more than 255 k-mers will break this program because of buffer sizes!\nPlease choose fewer k-mers");
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();
    env_logger::Builder::new()
        .filter_level(args.verbosity.log_level_filter())
        .init();
    let _total_tmr = timer!(Level::Info; "Total Runtime");
    let lookup_table = parser::parse_reference_fasta_file(&args.database_path)?;
    let query_data = parser::parse_query_fasta_file(&args.query_file)?;
    let output = args.get_output()?;
    let confidence_output = args.get_confidence_output()?;
    let result = sintax(&query_data, &lookup_table, &args);
    utils::output_results(&result, output, confidence_output, args.min_confidence)?;

    Ok(())
}
