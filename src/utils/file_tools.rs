// Various file tools needed throughout the code.

use anyhow::Result;
use log::warn;
use std::fs;
use std::fs::File;
use std::io::BufRead;
use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};

pub fn read_lines(filename: &PathBuf) -> io::Result<io::Lines<io::BufReader<File>>> {
    // This creates a buffer to read lines
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn open_file(filename: &PathBuf, overwrite_file: bool) -> Result<BufWriter<File>> {
    let file = if overwrite_file && filename.exists() {
        File::options().truncate(true).write(true).open(filename)
    } else {
        File::options().create_new(true).append(true).open(filename)
    };
    Ok(BufWriter::new(file?))
}

pub fn check_create_dir(path_to_check: &Path) -> Result<()> {
    if !path_to_check.is_dir() {
        warn!("Directory not found, creating: {:?}", path_to_check);
        fs::create_dir(path_to_check)?
    }
    Ok(())
}
