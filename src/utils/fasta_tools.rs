// This library contains tools needed to process fasta files as input and output.

use crate::utils::file_tools::open_file;
use crate::utils::nucleotides::{seq_to_string, string_to_seq};
use crate::utils::types::SeqByContig;
use anyhow::Result;
use log::info;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Read a fasta file and return a hashmap with the contig names as
/// keys and the sequences as values.
pub fn read_fasta(fasta_path: &PathBuf) -> Result<(SeqByContig, Vec<String>)> {
    info!("Reading fasta: {}", fasta_path.display());

    let genome: SeqByContig = SeqByContig::new();
    let file = File::open(fasta_path)?;
    let reader = BufReader::new(file);
    // Collect lines and handle errors
    let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;

    // compute windows for each contig
    let mut contig_wins: Vec<usize> = lines
        .iter()
        .enumerate()
        .filter(|(_, line)| line.starts_with('>'))
        .map(|(idx, _)| idx)
        .collect();

    // extract contig names
    let contigs: Vec<String> = contig_wins
        .iter()
        .map(|l| lines[*l].to_string())
        .map(|c| c.strip_prefix('>').unwrap().to_owned())
        .collect();
    // add eof
    contig_wins.push(lines.len());

    // closure for reading a contig
    let read_closure = |(n, contig_name): (usize, &String)| {
        let (start, end) = (1 + contig_wins[n], contig_wins[n + 1]);
        let mut seq = Vec::new();

        for line in &lines[start..end] {
            if let Ok(partial_seq) = string_to_seq(line) {
                seq.extend(partial_seq);
            }
        }

        // Insert sequence safely
        genome.entry(contig_name.clone()).or_default().extend(seq);
    };

    // Read contigs in parallel
    info!("Parsing fasta: {}", fasta_path.display());
    contigs.par_iter().enumerate().for_each(read_closure);

    Ok((genome, contigs))
}

/// Write a (mutated) genome to a fasta file.
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
        let sequence = &fasta_output.get(contig).unwrap();
        // Write contig name
        writeln!(outfile, ">{}", contig)?;
        // write sequences[ploid] to this_fasta
        for chunk in sequence.chunks(70) {
            let line: String = seq_to_string(chunk);
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
        let fasta_output: SeqByContig = SeqByContig::new();
        fasta_output.insert(String::from("H1N1_HA"), seq1);
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
