use anyhow::{bail, Context, Result};
use clap::Parser;
use clap_verbosity_flag::{Verbosity, WarnLevel};
use log::info;
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
};

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
    /// Path to the database fasta or bin file
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
    /// Output primary result file in tsv format
    #[arg(long)]
    pub tsv: bool,
    /// Create a binary database to load instead of a fasta file for repeated execution
    #[arg(long)]
    pub make_db: bool,
    /// Number of threads
    /// If 0, uses all available threads
    #[arg(short, long, default_value_t = 0, verbatim_doc_comment)]
    pub threads: usize,
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
    pub fn get_output(
        &self,
    ) -> Result<(
        Box<dyn Write>,
        Option<Box<dyn Write>>,
        Box<dyn Write>,
        Box<dyn Write + Send>,
    )> {
        let prefix = self
            .prefix
            .clone()
            .unwrap_or(self.query_file.clone().with_extension("out"));
        if prefix.is_dir() && !self.redo {
            bail!("Output folder {} already exists! Please specify another folder with -o <PATH> or run with --redo to force overriding existing files!", prefix.display());
        }
        std::fs::create_dir_all(&prefix)?;
        let result_output = prefix.join("raxtax.out");
        let tsv_output = if self.tsv {
            Some(
                std::fs::File::create(prefix.join("raxtax.tsv"))
                    .map(|f| Box::new(f) as Box<dyn Write>)?,
            )
        } else {
            None
        };
        let confidence_output =
            prefix.join(format!("raxtax.confidence{}.out", self.min_confidence));
        let log_output = prefix.join("raxtax.log");
        Ok((
            std::fs::File::create(result_output).map(|f| Box::new(f) as Box<dyn Write>)?,
            tsv_output,
            std::fs::File::create(confidence_output).map(|f| Box::new(f) as Box<dyn Write>)?,
            std::fs::File::create(log_output).map(|f| Box::new(f) as Box<dyn Write + Send>)?,
        ))
    }

    pub fn get_db_output(&self) -> Result<Box<dyn Write>> {
        let db_path = self.database_path.with_extension("bin");
        if db_path.is_file() && !self.redo {
            bail!("Output database file {} already exists! Delete it or run with --redo to force overriding existing files!", db_path.display());
        }
        info!("Created binary database at {}", db_path.display());
        Ok(
            std::fs::File::create(db_path)
                .map(|f| Box::new(BufWriter::new(f)) as Box<dyn Write>)?,
        )
    }
}
