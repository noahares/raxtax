use std::process::exit;

use clap::Parser;
use log::error;
use log::log_enabled;
use log::Level;
use logging_timer::timer;
use sintax_ng::io;
use sintax_ng::parser;
use sintax_ng::sintax::sintax;
use sintax_ng::utils;

fn main() {
    let args = io::Args::parse();
    let log_output = match args.get_log_output() {
        Ok(output) => output,
        Err(e) => {
            if args.verbosity.log_level_filter() >= Level::Error {
                eprintln!("\x1b[31m[ERROR]\x1b[0m {}", e);
            }
            exit(exitcode::CANTCREAT);
        }
    };
    env_logger::Builder::new()
        .target(env_logger::Target::Pipe(log_output))
        .filter_level(args.verbosity.log_level_filter())
        .format_timestamp(None)
        .format_target(false)
        .init();
    let output = match args.get_output() {
        Ok(output) => output,
        Err(e) => {
            error!("{}", e);
            if log_enabled!(Level::Error) {
                eprintln!("\x1b[31m[ERROR]\x1b[0m {}", e);
            }
            exit(exitcode::CANTCREAT);
        }
    };
    let confidence_output = match args.get_confidence_output() {
        Ok(output) => output,
        Err(e) => {
            error!("{}", e);
            if log_enabled!(Level::Error) {
                eprintln!("\x1b[31m[ERROR]\x1b[0m {}", e);
            }
            exit(exitcode::CANTCREAT);
        }
    };
    if let Err(e) = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
    {
        error!("{}", e);
        if log_enabled!(Level::Error) {
            eprintln!("\x1b[31m[ERROR]\x1b[0m {}", e);
        }
        exit(exitcode::OSERR);
    };
    let _total_tmr = timer!(Level::Info; "Total Runtime");
    let lookup_table = match parser::parse_reference_fasta_file(&args.database_path) {
        Ok(res) => res,
        Err(e) => {
            error!("Failed to parse {}: {}", &args.database_path.display(), e);
            if log_enabled!(Level::Error) {
                eprintln!(
                    "\x1b[31m[ERROR]\x1b[0m Failed to parse {}: {}",
                    &args.database_path.display(),
                    e
                );
            }
            exit(exitcode::NOINPUT);
        }
    };
    let query_data = match parser::parse_query_fasta_file(&args.query_file) {
        Ok(res) => res,
        Err(e) => {
            error!("Failed to parse {}: {}", &args.query_file.display(), e);
            if log_enabled!(Level::Error) {
                eprintln!(
                    "\x1b[31m[ERROR]\x1b[0m Failed to parse {}: {}",
                    &args.query_file.display(),
                    e
                );
            }
            exit(exitcode::NOINPUT);
        }
    };
    let result = sintax(
        &query_data,
        &lookup_table,
        args.num_k_mers,
        args.max_target_seqs,
        args.skip_exact_matches,
    );
    match utils::output_results(&result, output, confidence_output, args.min_confidence) {
        Ok(res) => res,
        Err(e) => {
            error!("Failed to write results to files: {}", e);
            if log_enabled!(Level::Error) {
                eprintln!(
                    "\x1b[31m[ERROR]\x1b[0m Failed to write results to files: {}",
                    e
                );
            }
            exit(exitcode::IOERR);
        }
    };
}
