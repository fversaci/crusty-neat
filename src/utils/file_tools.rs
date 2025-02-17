// Various file tools needed throughout the code.

use anyhow::{anyhow, Result};
use log::warn;
use std::fs::File;
use std::io::{BufRead, Error};
use std::path::Path;
use std::{fs, io};

pub fn read_lines(filename: &str) -> io::Result<io::Lines<io::BufReader<File>>> {
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

pub fn check_parent(filename: &str) -> Result<&Path> {
    // checks that the parent dir exists and then if so creates the Path object open
    // and ready to write
    let file_path = Path::new(filename);
    if !file_path
        .parent()
        .ok_or_else(|| anyhow!("no parent dir found"))?
        .exists()
    {
        check_create_dir(file_path)?;
    };
    Ok(file_path)
}

pub fn check_create_dir(path_to_check: &Path) -> Result<()> {
    if !path_to_check.is_dir() {
        warn!("Directory not found, creating: {:?}", path_to_check);
        fs::create_dir(path_to_check)?
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_parent() -> Result<()> {
        let filename = "test_data/H1N1.fa";
        check_parent(filename)?;
        Ok(())
    }

    #[test]
    fn test_check_parent_fail() -> Result<()> {
        let filename = "fake/test.fa";
        let er = check_parent(filename);
        assert!(er.is_err());
        Ok(())
    }
}
