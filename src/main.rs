use anyhow::Result;
use clap::Parser;
use itertools::Itertools;
use log::Level;
use logging_timer::timer;
use make_sintax_great_again::io;
use make_sintax_great_again::parser;
use make_sintax_great_again::sintax::sintax;
// use make_sintax_great_again::utils;

fn main() -> Result<()> {
    let args = io::Args::parse();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();
    env_logger::Builder::new()
        .filter_level(args.verbosity.log_level_filter())
        .init();
    let _total_tmr = timer!(Level::Info; "Total Runtime");
    let lookup_table = if let Some(ref database_path) = args.database_path {
        let database_str = std::fs::read_to_string(database_path)?;
        serde_json::from_str(&database_str)
    } else {
        let lookup_table = parser::parse_reference_fasta_file(&args.sequence_file)?;
        if let Some(ref database_output_path) = args.database_output {
            let output = std::fs::File::create(database_output_path)?;
            serde_json::to_writer(output, &lookup_table)?;
        }
        Ok(lookup_table)
    }?;
    let query_data = parser::parse_query_fasta_file(&args.query_file)?;
    let mut output = args.get_output()?;
    let result = sintax(&query_data, &lookup_table, &args);
    write!(output, "{}", result.into_iter().join("\n"))?;

    Ok(())
}
