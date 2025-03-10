// This library contains tools needed to process fasta files as input and output.

use super::file_tools::open_file;
use super::file_tools::read_lines;
use super::nucleotides::{base_to_nuc, nuc_to_base, Nuc};
use super::types::SeqByContig;
use anyhow::{anyhow, Result};
use log::info;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Reads a fasta file and returns a hashmap with the contig names as
/// keys and the sequences as values.
///
/// # Arguments
///
///   - fasta_path: the path to the fasta file
///
/// # Returns
///
///   - A tuple with the hashmap of sequences and the order of the
///     sequences in the file
///
/// Errors if the file cannot be read or if the file is not in fasta
/// format.
pub fn read_fasta(fasta_path: &PathBuf) -> Result<(SeqByContig, Vec<String>)> {
    info!("Reading fasta: {}", fasta_path.display());

    let mut fasta_map: SeqByContig = HashMap::new();
    let mut fasta_order: Vec<String> = Vec::new();
    let mut current_key = String::new();
    let mut temp_seq: Vec<Nuc> = Vec::new();

    let lines = read_lines(fasta_path)?;
    for line in lines {
        let l = line?;
        if l.starts_with('>') {
            if !current_key.is_empty() {
                fasta_map.entry(current_key).or_insert(temp_seq.clone());
            }
            current_key = l
                .strip_prefix('>')
                .ok_or_else(|| anyhow!("prefix not found"))?
                .to_string();
            fasta_order.push(current_key.clone());
            temp_seq.clear();
        } else {
            temp_seq.extend(l.chars().map(base_to_nuc).collect::<Result<Vec<_>, _>>()?);
        }
    }

    // Need to pick up the last one
    fasta_map.entry(current_key).or_insert(temp_seq);
    Ok((fasta_map, fasta_order))
}

/// Writes a hashmap of sequences to a fasta file.
///
/// # Arguments
///
///   - `fasta_output`: a hashmap with the contig names as keys and the
///     mutated sequences as values
///   - `fasta_order`: a vector with the order of the sequences in the
///     file
///   - `output_file`: the prefix for the output file name
///
/// # Returns
///
///   - Nothing
///
/// Errors if there is a problem writing the file.
pub fn write_fasta(
    fasta_output: &SeqByContig,
    fasta_order: &[String],
    overwrite_output: bool,
    output_prefix: &Path,
) -> Result<()> {
    // writing fasta output to files
    let output_fasta = output_prefix.with_extension("fasta");
    info!("Writing {}", output_fasta.display());
    let mut outfile = open_file(&output_fasta, overwrite_output)?;

    for contig in fasta_order {
        let sequence = &fasta_output[contig];
        // Write contig name
        writeln!(outfile, ">{}", contig)?;
        // write sequences[ploid] to this_fasta
        for chunk in sequence.chunks(70) {
            let line: String = chunk.iter().map(|&b| nuc_to_base(b)).collect();
            writeln!(outfile, "{}", line)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::{create_rng, utils::nucleotides::random_seq};

    use super::*;
    use std::fs;
    use tempdir::TempDir;

    #[test]
    fn test_read_fasta() -> Result<()> {
        let test_fasta: PathBuf = "test_data/H1N1.fa".into();
        let (_test_map, map_order) = read_fasta(&test_fasta)?;
        assert_eq!(map_order[0], "H1N1_HA".to_string());
        Ok(())
    }

    #[test]
    fn test_read_bad_fasta() -> Result<()> {
        let test_fasta: PathBuf = "test_data/fake.fasta".into();
        let er = read_fasta(&test_fasta);
        assert!(er.is_err());
        Ok(())
    }

    #[test]
    fn test_write_fasta() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let seq1 = random_seq(&mut rng, 100);
        let fasta_output: SeqByContig = HashMap::from([(String::from("H1N1_HA"), seq1)]);
        let fasta_pointer = fasta_output;
        let fasta_order = vec![String::from("H1N1_HA")];

        let tmp_dir = TempDir::new("crusty_neat")?;
        let output_file = tmp_dir.path().join("test");

        write_fasta(&fasta_pointer, &fasta_order, true, &output_file)?;
        let file_name = output_file.with_extension("fasta");
        let attr = fs::metadata(&file_name)?;
        assert!(attr.len() > 0);
        fs::remove_file(file_name)?;
        Ok(())
    }
}
