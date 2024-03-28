use std::process::exit;

use clap::Parser;
use log::Level;
use log::error;
use logging_timer::timer;
use sintax_ng::io;
use sintax_ng::parser;
use sintax_ng::sintax::sintax;
use sintax_ng::utils;

fn main() {
    let args = io::Args::parse();
    env_logger::Builder::new()
        .filter_level(args.verbosity.log_level_filter())
        .init();
    if args.num_k_mers >= 255 {
        error!("Using more than 255 k-mers will break this program because of buffer sizes! Please choose fewer k-mers");
        exit(exitcode::USAGE);
    }
    let output = match args.get_output() {
        Ok(output) => output,
        Err(e) => {
            error!("{}", e);
            exit(exitcode::CANTCREAT);
        },
    };
    let confidence_output = match args.get_confidence_output() {
        Ok(output) => output,
        Err(e) => {
            error!("{}", e);
            exit(exitcode::CANTCREAT);
        },
    };
    if let Err(e) = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global() {
            error!("{}", e);
            exit(exitcode::OSERR);
        };
    let _total_tmr = timer!(Level::Info; "Total Runtime");
    let lookup_table = match parser::parse_reference_fasta_file(&args.database_path) {
        Ok(res) => res,
        Err(e) => {
            error!("Failed to parse {}: {}", &args.database_path.display(), e);
            exit(exitcode::NOINPUT);
        },
    };
    let query_data = match parser::parse_query_fasta_file(&args.query_file) {
        Ok(res) => res,
        Err(e) => {
            error!("Failed to parse {}: {}", &args.query_file.display(), e);
            exit(exitcode::NOINPUT);
        },
    };
    let result = sintax(&query_data, &lookup_table, &args);
    match utils::output_results(&result, output, confidence_output, args.min_confidence) {
        Ok(res) => res,
        Err(e) => {
            error!("Failed to write results to files: {}", e);
            exit(exitcode::IOERR);
        },
    };
}
