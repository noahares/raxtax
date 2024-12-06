use std::process::exit;

use clap::Parser;
use log::Level;
use logging_timer::timer;
use raxtax::io;
use raxtax::parser;
use raxtax::raxtax::raxtax;
use raxtax::utils;

fn main() {
    // Parse args, set up files and other context
    let args = io::Args::parse();
    let io::OutputWriters {
        primary: output,
        tsv: tsv_output,
        log: log_output,
    } = args.get_output().unwrap_or_else(|e| {
        if args.verbosity.log_level_filter() >= Level::Error {
            eprintln!("\x1b[31m[ERROR]\x1b[0m {e}");
        }
        exit(exitcode::CANTCREAT);
    });
    env_logger::Builder::new()
        .target(env_logger::Target::Pipe(log_output))
        .filter_level(args.verbosity.log_level_filter())
        .format_timestamp(None)
        .format_target(false)
        .init();
    if args.pin {
        if let Err(e) = utils::setup_threadpool_pinned(args.threads) {
            utils::report_error(e, "Failed to set up thread pinning! Continuing without");
            if let Err(e) = rayon::ThreadPoolBuilder::new()
                .num_threads(args.threads)
                .build_global()
            {
                utils::report_error(anyhow::Error::from(e), "Failed to set up Multithreading");
                exit(exitcode::OSERR);
            }
        };
    } else if let Err(e) = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
    {
        utils::report_error(anyhow::Error::from(e), "Failed to set up Multithreading");
        exit(exitcode::OSERR);
    };
    if let Err(e) = std::fs::metadata(&args.query_file) {
        utils::report_error(
            anyhow::Error::from(e),
            format!("Failed to parse {}", args.query_file.display()),
        );
        exit(exitcode::NOINPUT);
    }
    let _total_tmr = timer!(Level::Info; "Total Runtime");

    // Parse reference databse
    let (store_db, tree) =
        parser::parse_reference_fasta_file(&args.database_path).unwrap_or_else(|e| {
            utils::report_error(
                e,
                format!("Failed to parse {}", args.database_path.display()),
            );
            exit(exitcode::NOINPUT);
        });
    if store_db && args.make_db {
        match args.get_db_output() {
            Ok(db_output) => {
                if let Err(e) = tree.save_to_file(db_output) {
                    utils::report_error(e, "Failed to write database");
                    exit(exitcode::IOERR);
                };
            }
            Err(e) => {
                utils::report_error(e, "Could not create database");
                exit(exitcode::CANTCREAT);
            }
        }
    }

    // Parse queries
    let queries = parser::parse_query_fasta_file(&args.query_file).unwrap_or_else(|e| {
        utils::report_error(e, format!("Failed to parse {}", args.query_file.display()));
        exit(exitcode::NOINPUT);
    });

    // Compute query results and output to files
    let n_threads = rayon::current_num_threads();
    let chunk_size = if n_threads == 1 {
        queries.len()
    } else {
        ((queries.len() / (n_threads * 10)) + 1).max(100)
    };
    let result = raxtax(
        &queries,
        &tree,
        args.skip_exact_matches,
        args.raw_confidence,
        chunk_size,
    );
    if let Some(tsv_output) = tsv_output {
        let sequences: Vec<String> = utils::decompress_sequences(&queries);
        if let Err(e) = utils::output_results_tsv(&result, sequences, tsv_output) {
            utils::report_error(e, "Failed to write results to file");
            exit(exitcode::IOERR);
        };
    }
    if let Err(e) = utils::output_results(&result, output) {
        utils::report_error(e, "Failed to write results to file");
        exit(exitcode::IOERR);
    };
    exit(exitcode::OK);
}
