use std::process::exit;

use anyhow::{anyhow, Result};
use clap::Parser;
use log::Level;
use logging_timer::timer;
use raxtax::io;
use raxtax::io::FileFingerprint;
use raxtax::parser;
use raxtax::raxtax::raxtax;
use raxtax::utils;
use std::io::Write;

fn main() {
    // Parse args, set up files and other context
    let args = io::Args::parse();
    let (
        io::OutputWriters {
            primary: mut output,
            tsv: mut tsv_output,
            log: mut log_output,
            progress: mut progress_output,
        },
        mut checkpoint,
    ) = args.get_output().unwrap_or_else(|e| {
        if args.verbosity.log_level_filter() >= Level::Error {
            eprintln!("\x1b[31m[ERROR]\x1b[0m {e}");
        }
        exit(exitcode::CANTCREAT);
    });
    if let Err(e) = io::write_build_info(&mut log_output) {
            eprintln!("\x1b[31m[ERROR]\x1b[0m {e}");
    }
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
    let _total_tmr = timer!(Level::Info; "Total Runtime");

    // Parse reference databse
    let (store_db, tree) = parser::parse_reference_fasta_file(&checkpoint.db_fingerprint.path)
        .unwrap_or_else(|e| {
            utils::report_error(
                e,
                format!(
                    "Failed to parse {}",
                    checkpoint.db_fingerprint.path.display()
                ),
            );
            exit(exitcode::NOINPUT);
        });
    if store_db && !args.skip_db {
        match args.get_db_output() {
            Ok((db_output, db_path)) => {
                if let Err(e) = tree.save_to_file(db_output) {
                    utils::report_error(e, "Failed to write database");
                    exit(exitcode::IOERR);
                };
                if let Err(e) = FileFingerprint::new(&db_path).and_then(|fp| {
                    checkpoint.db_fingerprint = fp;
                    checkpoint.save()?;
                    Ok(())
                }) {
                    utils::report_error(e, "Failed to write checkpoint! Continuing without...");
                };
            }
            Err(e) => {
                utils::report_error(
                    e,
                    "Could not create database! Rerun with --skip-db to skip this step.",
                );
                exit(exitcode::CANTCREAT);
            }
        }
    } else {
        checkpoint.save().unwrap_or_else(|e| {
            utils::report_error(e, "Failed to write checkpoint! Continuing without...")
        });
    }

    if args.only_db {
        exit(exitcode::OK);
    }

    // Parse queries
    let queries = parser::parse_query_fasta_file(
        args.query_file.as_ref().unwrap(),
        &checkpoint.processed_queries,
    )
    .unwrap_or_else(|e| {
        utils::report_error(
            e,
            format!("Failed to parse {}", args.query_file.unwrap().display()),
        );
        exit(exitcode::NOINPUT);
    });

    // Compute query results and output to files
    let n_threads = rayon::current_num_threads();
    let chunk_size = if n_threads == 1 {
        queries.len()
    } else {
        ((queries.len() / (n_threads * 10)) + 1).max(100)
    };

    let (sender, receiver) = crossbeam::channel::unbounded::<(String, String, Option<String>)>();
    let writer_handle = std::thread::spawn(move || -> Result<()> {
        for (query, results, tsv_results) in receiver {
            if let Some(ref mut tsv_output) = tsv_output {
                writeln!(tsv_output, "{}", tsv_results.unwrap())?;
            }
            writeln!(output, "{}", results)?;
            writeln!(progress_output, "{}", query)?;
        }
        Ok(())
    });
    let ok = raxtax(
        &queries,
        &tree,
        args.skip_exact_matches,
        args.raw_confidence,
        chunk_size,
        &sender,
        args.tsv,
    );
    drop(sender);
    if writer_handle.join().is_err() {
        utils::report_error(
            anyhow!("IO-thread could not be joined. Check if results are complete!"),
            "",
        );
    };
    if let Err(e) = ok {
        utils::report_error(
            e,
            "Error while sending results to IO-thread!\n
            Rerun raxtax to continue from the last checkpoint.\n
            If the problem persists, please report this issue at: https://github.com/noahares/raxtax/issues",
        );
        exit(exitcode::TEMPFAIL);
    }

    if args.clean {
        checkpoint.cleanup().unwrap_or_else(|e| {
            utils::report_error(
                e,
                "Removing checkpoint files failed! Please delete them manually.",
            )
        });
    }

    exit(exitcode::OK);
}
