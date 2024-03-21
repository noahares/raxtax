use anyhow::Result;
use clap::Parser;
use clap_verbosity_flag::Verbosity;
use std::{io::Write, path::PathBuf};

#[derive(Parser)]
#[command(author, version, about)]
pub struct Args {
    /// Path to the database fasta file
    #[arg(short, long)]
    pub database_path: PathBuf,
    /// Path to the query file
    #[arg(short = 'i', long)]
    pub query_file: PathBuf,
    /// Number of rounds per query
    #[arg(short, long, default_value_t = 100)]
    pub num_iterations: usize,
    /// Number of 8-mers
    #[arg(short = 'k', long, default_value_t = 32)]
    pub num_k_mers: usize,
    /// 8-mer hit-threshold
    #[arg(short = 'f', long, default_value_t = 1.0 / 3.0)]
    pub min_hit_fraction: f64,
    /// Confidence threshold
    #[arg(short = 'c', long, default_value_t = 0.8)]
    pub min_confidence: f64,
    /// Number of output species per query
    #[arg(short = 'm', long, default_value_t = 5)]
    pub max_target_seqs: usize,
    /// Number of threads
    #[arg(short, long, default_value_t = 0)]
    pub threads: usize,
    /// Seed
    #[arg(long, default_value_t = 42)]
    pub seed: u64,
    /// Output path
    #[arg(short, long)]
    pub output: Option<PathBuf>,
    /// confidence output path
    #[arg(short = 'u', long)]
    pub confidence_output: Option<PathBuf>,
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
    pub fn get_confidence_output(&self) -> Result<Box<dyn Write>> {
        match self.confidence_output {
            Some(ref path) => {
                Ok(std::fs::File::create(path).map(|f| Box::new(f) as Box<dyn Write>)?)
            }
            None => Ok(Box::new(std::io::stdout())),
        }
    }
}
