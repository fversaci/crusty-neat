extern crate clap;
extern crate itertools;
extern crate log;
extern crate serde;
extern crate serde_json;
extern crate serde_yaml;
extern crate simplelog;
extern crate statrs;

mod utils;

use anyhow::{anyhow, Result};
use clap::Parser;
use log::*;
use rand::rngs::StdRng;
use rand::SeedableRng;
use simplelog::*;
use std::fs::File;
use std::hash::Hash;
use std::hash::{DefaultHasher, Hasher};
use utils::cli;
use utils::config::{build_config_from_args, read_config_yaml};
use utils::file_tools::check_parent;
use utils::runner::run_neat;

/// Create a random number generator from a seed string. If no seed is provided
/// we generate a random seed.
pub fn create_rng(seed: Option<&str>) -> StdRng {
    let seed = seed
        .map(|s| {
            let mut hasher = DefaultHasher::new();
            s.hash(&mut hasher);
            hasher.finish()
        })
        .unwrap_or_else(rand::random);

    StdRng::seed_from_u64(seed)
}

/// Main function for the program. This function parses the command
/// line arguments and then runs the main script for generating reads.
fn main() -> Result<()> {
    // parse the arguments from the command line
    let args = cli::Cli::parse();

    // log filter
    let level_filter = match args.log_level.to_lowercase().as_str() {
        "trace" => LevelFilter::Trace,
        "debug" => LevelFilter::Debug,
        "info" => LevelFilter::Info,
        "warn" => LevelFilter::Warn,
        "error" => LevelFilter::Error,
        "off" => LevelFilter::Off,
        _ => {
            return Err(anyhow!(
                "Unknown log level, please set to one of \
             Trace, Debug, Info, Warn, Error, or Off (case insensitive)."
            ))
        }
    };

    // Check that the parent dir exists
    let log_destination = check_parent(&args.log_dest)?;

    // Set up the logger for the run
    CombinedLogger::init(vec![
        TermLogger::new(
            level_filter,
            Config::default(),
            TerminalMode::Stdout,
            ColorChoice::Always,
        ),
        WriteLogger::new(
            level_filter,
            Config::default(),
            File::create(log_destination)?,
        ),
    ])?;

    info!("Begin processing");

    // set up the config struct based on whether there was an input
    // config. Input config overrides any other inputs.
    let config = if !args.config.is_empty() {
        info!("Using Configuration file input: {}", &args.config);
        read_config_yaml(args.config)
    } else {
        info!("Using command line arguments.");
        debug!("Command line args: {:?}", &args);
        build_config_from_args(args)
    }?;

    // Generate the RNG used for this run
    let seed = config.rng_seed.as_deref();
    let mut rng = create_rng(seed);

    if let Some(sd) = seed {
        info!("Seed string to regenerate these exact results: {}", sd);
    }

    // Run the main script
    run_neat(config, &mut rng)
}
