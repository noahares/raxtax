use chrono::prelude::*;
use std::process::Command;

fn main() {
    // Get git commit hash
    let git_hash = Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| "unknown".into());
    println!("cargo:rustc-env=GIT_HASH={}", git_hash);

    // --- Determine active Cargo profile ---
    let profile = std::env::var("PROFILE").unwrap_or_else(|_| "unknown".into());
    println!("cargo:rustc-env=BUILD_PROFILE={}", profile);

    // --- Capture RUSTFLAGS (from .cargo/config.toml or env) ---
    if let Ok(encoded) = std::env::var("CARGO_ENCODED_RUSTFLAGS") {
        let flags: Vec<_> = encoded.split('\x1f').collect();
        println!("cargo:rustc-env=RUSTFLAGS={}", flags.join(" "));
    } else if let Ok(rustflags) = std::env::var("RUSTFLAGS") {
        println!("cargo:rustc-env=RUSTFLAGS={}", rustflags);
    } else {
        println!("cargo:rustc-env=RUSTFLAGS=");
    }
    println!(
        "cargo:rustc-env=BUILD_DATE={}",
        Utc::now().to_rfc3339_opts(SecondsFormat::Secs, false)
    );
}
