mod utils;

use anyhow::Result;
use clap::Parser;
use log::*;
use rand::SeedableRng;
use rand::rngs::StdRng;
use simplelog::*;
use std::fs::File;
use std::hash::Hash;
use std::hash::{DefaultHasher, Hasher};
use utils::config::{self, RunConfiguration};
use utils::runner::run_neat;

/// Create a random number generator from an optional seed string.
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

/// Parse the command line arguments and merge them with the config file.
fn main() -> Result<()> {
    // parse cli arguments
    let args = config::Args::parse();
    // read config file or start with default values
    let mut config;
    if let Some(file) = args.config_file {
        info!("Reading configuration from file: {}", file.display());
        config = RunConfiguration::from_file(&file)?;
    } else {
        info!("Using default configuration");
        config = RunConfiguration::fill();
    }
    // override values from the command line (if any)
    if let Some(arg_config) = args.config {
        info!("Overriding configuration from command line");
        config.override_with(&arg_config)?;
    }

    // start loggers
    let _loggers;
    let term_log = TermLogger::new(
        args.log_level,
        Config::default(),
        TerminalMode::Stdout,
        ColorChoice::Always,
    );
    if let Some(log_path) = args.log_dest {
        let flog = WriteLogger::new(args.log_level, Config::default(), File::create(log_path)?);
        _loggers = CombinedLogger::init(vec![term_log, flog])?;
    } else {
        _loggers = CombinedLogger::init(vec![term_log])?;
    }

    info!("Begin processing");
    debug!("RunConfiguration: {:?}", config);

    // Generate the RNG used for this run
    let seed = config.rng_seed.as_deref();
    let mut rng = create_rng(seed);

    if let Some(sd) = seed {
        info!("Seed string to regenerate these exact results: {:?}", sd);
    }

    // Run the main function
    run_neat(config, &mut rng)
}
