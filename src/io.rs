use anyhow::{bail, Result};
use clap::Parser;
use clap_verbosity_flag::{InfoLevel, Verbosity};
use log::info;
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
};

pub struct OutputWriters {
    pub primary: Box<dyn Write>,
    pub tsv: Option<Box<dyn Write>>,
    pub log: Box<dyn Write + Send>,
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
    /// If used for mislabling analysis, you want to skip exact sequence matches
    #[arg(long)]
    pub skip_exact_matches: bool,
    /// Output primary result file in tsv format
    #[arg(long)]
    pub tsv: bool,
    /// Create a binary database to load instead of a fasta file for repeated execution
    #[arg(long)]
    pub make_db: bool,
    /// Don't adjust confidence values for 1 exact match
    #[arg(long)]
    pub raw_confidence: bool,
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
    /// Use thread pinning
    #[arg(long)]
    pub pin: bool,
    #[command(flatten)]
    pub verbosity: Verbosity<InfoLevel>,
}

impl Args {
    pub fn get_output(&self) -> Result<OutputWriters> {
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
        let log_output = prefix.join("raxtax.log");
        Ok(OutputWriters {
            primary: std::fs::File::create(result_output).map(|f| Box::new(f) as Box<dyn Write>)?,
            tsv: tsv_output,
            log: std::fs::File::create(log_output).map(|f| Box::new(f) as Box<dyn Write + Send>)?,
        })
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
