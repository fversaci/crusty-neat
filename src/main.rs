extern crate clap;
extern crate itertools;
extern crate log;
extern crate serde;
extern crate serde_json;
extern crate serde_yaml;
extern crate simple_rng;
extern crate simplelog;
extern crate statrs;

mod utils;

use anyhow::{anyhow, Result};
use chrono::Utc;
use clap::Parser;
use log::*;
use simple_rng::Rng;
use simplelog::*;
use std::fs::File;
use utils::cli;
use utils::config::{build_config_from_args, read_config_yaml};
use utils::file_tools::check_parent;
use utils::runner::run_neat;

fn main() -> Result<()> {
    info!("Begin processing");
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
        _ => panic!(
            "Unknown log level, please set to one of \
            Trace, Debug, Info, Warn, Error, or Off (case insensitive)."
        ),
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
        SimpleLogger::new(LevelFilter::Trace, Config::default()),
        WriteLogger::new(
            level_filter,
            Config::default(),
            File::create(log_destination)?,
        ),
    ])?;
    // set up the config struct based on whether there was an input config. Input config
    // overrides any other inputs.
    let config = if !args.config.is_empty() {
        info!("Using Configuration file input: {}", &args.config);
        read_config_yaml(args.config)
    } else {
        info!("Using command line arguments.");
        debug!("Command line args: {:?}", &args);
        build_config_from_args(args)
    }?;
    // Generate the RNG used for this run. If not we generate a random seed using the current time
    let mut seed_vec: Vec<String> = Vec::new();
    if config.rng_seed.is_some() {
        // Force read it as a string, hopefully
        let raw_seed = config
            .rng_seed
            .clone()
            .ok_or_else(|| anyhow!("cannot clone seed"))?
            .to_string();
        for seed_term in raw_seed.split_whitespace() {
            seed_vec.push(seed_term.to_string());
        }
        info!(
            "Seed string to regenerate these exact results: {}",
            raw_seed
        );
    } else {
        // since no seed was provided, we'll use the datetime stamp
        info!(
            "No rng seed provided, using timestamp (and space-separated list of words with simple \
            characters will also work as a key)"
        );
        let timestamp = Utc::now().format("%Y %m %d %H %M %S %f").to_string();
        for item in timestamp.split_whitespace() {
            seed_vec.push(item.to_string());
        }
        info!(
            "Seed string to regenerate these exact results: {}",
            timestamp
        );
    }
    let mut rng: Rng = Rng::from_seed(seed_vec);
    // run the generate reads main script
    run_neat(config, &mut rng)
}
