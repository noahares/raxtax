use ahash::{HashSet, HashSetExt};
use anyhow::{bail, Result};
use clap::Parser;
use clap_verbosity_flag::{InfoLevel, Verbosity};
use log::{info, Level};
use serde::{Deserialize, Serialize};
use std::{
    fs::{File, OpenOptions},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use crate::utils;

pub struct OutputWriters {
    pub primary: BufWriter<File>,
    pub tsv: Option<BufWriter<File>>,
    pub log: Box<dyn Write + Send>,
}

#[derive(Serialize, Deserialize, PartialEq, Eq)]
pub struct FileFingerprint {
    pub path: PathBuf,
    size: u64,
    modified: u64,
}

impl FileFingerprint {
    pub fn new(path: &PathBuf) -> Result<FileFingerprint> {
        let metadata = std::fs::metadata(path)?;
        let size = metadata.len();
        let modified = metadata
            .modified()?
            .duration_since(std::time::UNIX_EPOCH)?
            .as_secs();
        Ok(FileFingerprint {
            path: std::path::absolute(path)?,
            size,
            modified,
        })
    }
}

#[derive(Serialize, Deserialize)]
pub struct Checkpoint {
    pub checkpoint_file: PathBuf,
    pub db_fingerprint: FileFingerprint,
    raw_confidence: bool,
    skip_exact_matches: bool,
    tsv: bool,
    pub processed_queries: HashSet<String>,
}

impl Checkpoint {
    pub fn new(own_path: &Path, args: &Args) -> Result<Checkpoint> {
        Ok(Checkpoint {
            checkpoint_file: std::path::absolute(own_path)?,
            db_fingerprint: FileFingerprint::new(&args.database_path)?,
            raw_confidence: args.raw_confidence,
            skip_exact_matches: args.skip_exact_matches,
            tsv: args.tsv,
            processed_queries: HashSet::new(),
        })
    }

    pub fn save(&self) -> Result<()> {
        let tmp_ckp_path = self.checkpoint_file.with_extension("ckp.tmp");
        let tmp_ckp_file = std::fs::File::create(&tmp_ckp_path)?;
        serde_json::to_writer_pretty(tmp_ckp_file, self)?;
        std::fs::rename(tmp_ckp_path, &self.checkpoint_file)?;
        Ok(())
    }
}

#[derive(Parser)]
#[command(author, version, about)]
pub struct Args {
    /// Path to the database fasta or bin file
    #[arg(short, long)]
    pub database_path: PathBuf,
    /// Path to the query file
    #[arg(short = 'i', long, required_unless_present = "only_db")]
    pub query_file: Option<PathBuf>,
    /// If used for mislabling analysis, you want to skip exact sequence matches
    #[arg(long)]
    pub skip_exact_matches: bool,
    /// Output primary result file in tsv format
    #[arg(long)]
    pub tsv: bool,
    /// Create binary database and exit
    #[arg(long, conflicts_with = "skip_db")]
    pub only_db: bool,
    /// Don't create the binary database for the reference sequences
    #[arg(long)]
    pub skip_db: bool,
    /// Don't adjust confidence values for 1 exact match
    #[arg(long)]
    pub raw_confidence: bool,
    /// Number of threads
    /// If 0, uses all available threads
    #[arg(short, long, default_value_t = 0, verbatim_doc_comment)]
    pub threads: usize,
    /// Output prefix
    #[arg(short = 'o', long, default_value = "raxtax")]
    pub prefix: PathBuf,
    /// Force override of existing output files
    #[arg(long)]
    pub redo: bool,
    /// Use thread pinning
    #[arg(long)]
    pub pin: bool,
    #[command(flatten)]
    pub verbosity: Verbosity<InfoLevel>,
}

fn create_file(path: PathBuf, append: bool) -> std::io::Result<BufWriter<File>> {
    if append {
        OpenOptions::new()
            .append(true)
            .create(true)
            .open(path)
            .map(BufWriter::new)
    } else {
        OpenOptions::new()
            .write(true)
            .truncate(true)
            .create(true)
            .open(path)
            .map(BufWriter::new)
    }
}

impl Args {
    pub fn get_output(&self) -> Result<(OutputWriters, Checkpoint)> {
        let prefix = self.get_prefix();
        let ckp_path = prefix.join("raxtax.ckp");
        let checkpoint = if !self.redo && ckp_path.is_file() {
            let ckp_file = std::fs::File::open(&ckp_path)?;
            match serde_json::from_reader(ckp_file) {
                Ok(ckp) => {
                    if self.checkpoint_valid(&ckp) {
                        ckp
                    } else {
                        Checkpoint::new(&ckp_path, self)?
                    }
                }
                Err(e) => {
                    utils::report_error(e.into(), "Failed to read checkpoint!");
                    Checkpoint::new(&ckp_path, self)?
                }
            }
        } else {
            Checkpoint::new(&ckp_path, self)?
        };
        if prefix.is_dir() && !ckp_path.is_file() && !self.redo {
            bail!("Output folder {} already exists! Please specify another folder with -o <PATH> or run with --redo to force overriding existing files!", prefix.display());
        }
        std::fs::create_dir_all(&prefix)?;
        let tsv_output = if self.tsv {
            Some(create_file(prefix.join("raxtax.tsv"), !self.redo)?)
        } else {
            None
        };
        let log_path = prefix.join("raxtax.log");
        let mut log_output =
            create_file(log_path, !self.redo).map(|f| Box::new(f) as Box<dyn Write + Send>)?;
        if !self.redo && ckp_path.exists() && self.verbosity.log_level_filter() >= Level::Info {
            writeln!(
                log_output,
                "[INFO ] Restarting from checkpoint {}",
                checkpoint.checkpoint_file.display()
            )?;
            eprintln!(
                "[INFO ] Restarting from checkpoint {}",
                checkpoint.checkpoint_file.display()
            );
        }
        Ok((
            OutputWriters {
                primary: create_file(prefix.join("raxtax.out"), !self.redo)?,
                tsv: tsv_output,
                log: log_output,
            },
            checkpoint,
        ))
    }

    fn get_prefix(&self) -> PathBuf {
        self.prefix.clone()
    }

    pub fn get_db_output(&self) -> Result<(BufWriter<File>, PathBuf)> {
        let prefix = self.get_prefix();
        let db_path = prefix
            .join(
                self.database_path
                    .file_name()
                    .unwrap_or_else(|| std::ffi::OsStr::new("database")),
            )
            .with_extension("bin");
        if db_path.is_file() && !self.redo {
            bail!("Output database file {} already exists! Delete it or run with --redo to force overriding existing files!", db_path.display());
        }
        info!("Created binary database at {}", db_path.display());
        Ok((
            std::fs::File::create(db_path.clone()).map(BufWriter::new)?,
            db_path,
        ))
    }

    fn checkpoint_valid(&self, checkpoint: &Checkpoint) -> bool {
        let reevaluated_fingerprint = FileFingerprint::new(&checkpoint.db_fingerprint.path);
        match reevaluated_fingerprint {
            Ok(fp) => {
                self.tsv == checkpoint.tsv
                    && self.raw_confidence == checkpoint.raw_confidence
                    && self.skip_exact_matches == checkpoint.skip_exact_matches
                    && fp == checkpoint.db_fingerprint
            }
            Err(e) => {
                utils::report_error(e, "Could not verify checkpoint, starting from scratch!");
                false
            }
        }
    }
}
