use anyhow::{bail, Context, Result};
use clap::Parser;
use clap_verbosity_flag::{Verbosity, WarnLevel};
use std::{io::Write, path::PathBuf};

fn normalized_ratio(s: &str) -> Result<f64> {
    let ratio: f64 = s
        .parse()
        .with_context(|| format!("`{s}` isn't a valid fraction"))?;
    if (0.0..=1.0).contains(&ratio) {
        Ok(ratio)
    } else {
        bail!("Fraction is not in range {:.2}-{:.2}", 0.0, 1.0)
    }
}

#[derive(Parser)]
#[command(author, version, about)]
pub struct Args {
    /// Path to the database fasta file
    #[arg(short, long)]
    pub database_path: PathBuf,
    /// Path to the query file
    #[arg(short = 'i', long)]
    pub query_file: PathBuf,
    /// Confidence threshold
    #[arg(short = 'c', long, default_value_t = 0.8, value_parser = normalized_ratio)]
    pub min_confidence: f64,
    /// If used for mislabling analysis, you want to skip exact sequence matches
    #[arg(long)]
    pub skip_exact_matches: bool,
    /// Number of threads
    /// If 0, uses all available threads
    #[arg(short, long, default_value_t = 0, verbatim_doc_comment)]
    pub threads: usize,
    /// Seed
    #[arg(short, long, default_value_t = 42)]
    pub seed: u64,
    /// Output prefix
    #[arg(short = 'o', long)]
    pub prefix: Option<PathBuf>,
    /// Force override of existing output files
    #[arg(long)]
    pub redo: bool,
    #[command(flatten)]
    pub verbosity: Verbosity<WarnLevel>,
}

impl Args {
    pub fn get_output(&self) -> Result<(Box<dyn Write>, Box<dyn Write>, Box<dyn Write + Send>)> {
        let prefix = self.prefix.clone().unwrap_or(self.query_file.clone().with_extension("out"));
        if prefix.is_dir() && !self.redo {
            bail!("Output folder {} already exists! Please specify another folder with -o <PATH> or run with --redo to force overriding existing files!", prefix.display());
        }
        std::fs::create_dir_all(&prefix)?;
        let result_output = prefix.join("sintax_ng.out");
        let confidence_output =
            prefix.join(format!("sintax_ng.confidence{}.out", self.min_confidence));
        let log_output = prefix.join("sintax_ng.log");
        Ok((
            std::fs::File::create(result_output).map(|f| Box::new(f) as Box<dyn Write + Send>)?,
            std::fs::File::create(confidence_output)
                .map(|f| Box::new(f) as Box<dyn Write + Send>)?,
            std::fs::File::create(log_output).map(|f| Box::new(f) as Box<dyn Write + Send>)?,
        ))
    }
}
