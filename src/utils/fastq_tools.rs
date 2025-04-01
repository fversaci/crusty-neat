// This library writes either single ended or paired-ended fastq files.

use crate::utils::file_tools::open_file;
use crate::utils::nucleotides::{Nuc, reverse_complement, seq_to_string};
use crate::utils::quality_model::QualityModel;
use anyhow::Result;
use rand::Rng;
use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use rayon::prelude::ParallelSlice;
use std::fs::create_dir_all;
use std::io::Write;
use std::path::Path;

fn write_sequence<W: Write>(
    outfile: &mut W,
    name_prefix: &str,
    order_index: usize,
    sequence: &[Nuc],
    quality_score_model: &QualityModel,
    rng: &mut impl Rng,
    read_number: u8,
) -> Result<()> {
    let read_length = sequence.len();
    let quality_scores = quality_score_model.generate_quality_scores(read_length, rng)?;

    // sequence id
    writeln!(
        outfile,
        "@{}{}/{}",
        name_prefix,
        order_index + 1,
        read_number
    )?;
    // sequence
    writeln!(outfile, "{}", seq_to_string(sequence))?;
    writeln!(outfile, "+")?;
    // quality scores
    writeln!(outfile, "{}", quality_scores)?;

    Ok(())
}

/// Writes the input reads as FASTQ files (in parallel).
///
/// # Arguments
///
/// * `fastq_filename` - The prefix for the output FASTQ files.
/// * `overwrite_output` - A boolean flag to enable or disable overwriting
///   existing files.
/// * `paired_ended` - A boolean flag to enable or disable paired-end mode.
/// * `reads` - A list of `Vec<Nuc>` representing DNA sequences.
/// * `quality_score_model` - A `QualityScoreModel` object to generate
///   quality scores.
/// * `rng` - A random number generator.
pub fn write_fastq<R: Rng + Send + Sync + Clone>(
    output_prefix: &Path,
    overwrite_output: bool,
    paired_ended: bool,
    reads: Vec<&[Nuc]>,
    quality_score_model: QualityModel,
    rng: &mut R,
) -> Result<()> {
    let dir1 = output_prefix.with_extension("reads_R1");
    let dir2 = output_prefix.with_extension("reads_R2");

    create_dir_all(&dir1)?;
    if paired_ended {
        create_dir_all(&dir2)?;
    }
    let max_chunk = 100000;
    reads
        .par_chunks(max_chunk)
        .enumerate()
        .try_for_each(|(chunk_index, chunk)| -> Result<()> {
            if chunk.is_empty() {
                return Ok(());
            }
            let name_prefix = format!("crusty_neat_generated_{}_", chunk_index);
            let file1_path = dir1.join(format!("{}.fastq", chunk_index));
            let mut file1 = open_file(&file1_path, overwrite_output)?;

            for (order_index, sequence) in chunk.iter().enumerate() {
                write_sequence(
                    &mut file1,
                    &name_prefix,
                    order_index,
                    sequence,
                    &quality_score_model,
                    &mut rng.clone(),
                    1,
                )?;
            }

            if paired_ended {
                let file2_path = dir2.join(format!("{}.fastq", chunk_index));
                let mut file2 = open_file(&file2_path, overwrite_output)?;
                for (order_index, sequence) in chunk.iter().enumerate() {
                    let rev_comp_sequence = reverse_complement(sequence);
                    write_sequence(
                        &mut file2,
                        &name_prefix,
                        order_index,
                        &rev_comp_sequence,
                        &quality_score_model,
                        &mut rng.clone(),
                        2,
                    )?;
                }
            }

            Ok(())
        })
}
