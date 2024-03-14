use anyhow::Result;
use clap::Parser;
use clap_verbosity_flag::Verbosity;
use std::{io::Write, path::PathBuf};

#[derive(Parser)]
#[command(author, version, about)]
pub struct Args {
    /// Path to the sequence file
    #[arg(short = 'i', long)]
    pub sequence_file: PathBuf,
    /// Path to the existing database
    #[arg(short = 'd', long)]
    pub database_path: Option<PathBuf>,
    /// Path to the database output
    #[arg(short = 'z', long)]
    pub database_output: Option<PathBuf>,
    /// Path to the query file
    #[arg(short, long)]
    pub query_file: PathBuf,
    /// Number of rounds per query
    #[arg(short, long, default_value_t = 100)]
    pub num_rounds: usize,
    /// Number of 8-mers
    #[arg(short = 'k', long, default_value_t = 32)]
    pub num_k_mers: usize,
    /// 8-mer hit-threshold
    #[arg(short = 't', long, default_value_t = 1.0 / 3.0)]
    pub threshold: f64,
    /// Number of output species per query
    #[arg(short = 'r', long, default_value_t = 5)]
    pub num_results: usize,
    /// Seed
    #[arg(short, long, default_value_t = 42)]
    pub seed: u64,
    /// Output path
    #[arg(short, long)]
    pub output: Option<PathBuf>,
    #[command(flatten)]
    pub verbosity: Verbosity,
}

impl Args {
    pub fn get_output(&self) -> Result<Box<dyn Write>> {
        match self.output {
            Some(ref path) => {
                Ok(std::fs::File::create(path).map(|f| Box::new(f) as Box<dyn Write>)?)
            }
            None => Ok(Box::new(std::io::stdout())),
        }
    }
}
