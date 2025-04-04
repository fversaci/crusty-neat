use anyhow::Result;
use flate2::{Compression, write::GzEncoder};
use log::warn;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub enum FileWriter {
    Plain(BufWriter<File>),
    Gzipped(BufWriter<GzEncoder<File>>),
}

impl Write for FileWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            FileWriter::Plain(w) => w.write(buf),
            FileWriter::Gzipped(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            FileWriter::Plain(w) => w.flush(),
            FileWriter::Gzipped(w) => w.flush(),
        }
    }
}

pub fn open_file(filename: &PathBuf, overwrite_file: bool, gzipped: bool) -> Result<FileWriter> {
    let file = if overwrite_file && filename.exists() {
        File::options().truncate(true).write(true).open(filename)
    } else {
        File::options().create_new(true).append(true).open(filename)
    };
    let file = file?;
    let writer = if gzipped {
        let encoder = GzEncoder::new(file, Compression::default());
        FileWriter::Gzipped(BufWriter::new(encoder))
    } else {
        FileWriter::Plain(BufWriter::new(file))
    };
    Ok(writer)
}

pub fn check_create_dir(path_to_check: &Path) -> Result<()> {
    if !path_to_check.is_dir() {
        warn!("Directory not found, creating: {:?}", path_to_check);
        fs::create_dir(path_to_check)?
    }
    Ok(())
}
