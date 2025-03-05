// Various file tools needed throughout the code.

use anyhow::Result;
use log::warn;
use std::fs::File;
use std::io::{BufRead, Error};
use std::path::Path;
use std::path::PathBuf;
use std::{fs, io};

pub fn read_lines(filename: &PathBuf) -> io::Result<io::Lines<io::BufReader<File>>> {
    // This creates a buffer to read lines
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn open_file(mut filename: &mut str, overwrite_file: bool) -> Result<File, Error> {
    if overwrite_file && Path::new(filename).exists() {
        File::options()
            .truncate(true)
            .write(true)
            .open(&mut filename)
    } else {
        File::options()
            .create_new(true)
            .append(true)
            .open(&mut filename)
    }
}

pub fn check_create_dir(path_to_check: &Path) -> Result<()> {
    if !path_to_check.is_dir() {
        warn!("Directory not found, creating: {:?}", path_to_check);
        fs::create_dir(path_to_check)?
    }
    Ok(())
}
