// This library writes either single ended or paired-ended fastq files.

use crate::utils::file_tools::open_file;
use crate::utils::nucleotides::{Nuc, reverse_complement, seq_to_string};
use crate::utils::quality_model::QualityModel;
use anyhow::Result;
use rand::Rng;
use std::io::Write;
use std::path::Path;

/// Writes FASTQ files based on the provided dataset.
///
/// # Arguments
///
/// * `fastq_filename` - The prefix for the output FASTQ files.
/// * `overwrite_output` - A boolean flag to enable or disable overwriting
///   existing files.
/// * `paired_ended` - A boolean flag to enable or disable paired-end mode.
/// * `dataset` - A list of `Vec<Nuc>` representing DNA sequences.
/// * `dataset_order` - A list of indices to order the dataset.
/// * `quality_score_model` - A `QualityScoreModel` object to generate
///   quality scores.
/// * `rng` - A random number generator.
pub fn write_fastq<R: Rng>(
    output_prefix: &Path,
    overwrite_output: bool,
    paired_ended: bool,
    dataset: Vec<&Vec<Nuc>>,
    dataset_order: Vec<usize>,
    quality_score_model: QualityModel,
    rng: &mut R,
) -> Result<()> {
    // The prefix for read names. Reads are numbered in output order.
    let name_prefix = "neat_generated_";
    let filename1 = output_prefix.with_extension("r1.fastq");
    let mut outfile1 = open_file(&filename1, overwrite_output)?;
    // Setting up pairend ended reads. For single ended reads,
    // we have outfile2 = None
    let mut outfile2 = if paired_ended {
        let filename2 = output_prefix.with_extension("r2.fastq");
        Some(open_file(&filename2, overwrite_output)?)
    } else {
        None
    };

    // Write sequences. Ordered index is used for numbering, while
    // read_index is from the shuffled index array from a previous
    // step.
    for (order_index, &read_index) in dataset_order.iter().enumerate() {
        let sequence = dataset[read_index];
        write_sequence(
            &mut outfile1,
            name_prefix,
            order_index,
            sequence,
            &quality_score_model,
            rng,
            1,
        )?;

        // Handle paired-end reads
        if let Some(ref mut outfile2) = outfile2 {
            let rev_comp_sequence = reverse_complement(sequence);
            write_sequence(
                outfile2,
                name_prefix,
                order_index,
                &rev_comp_sequence,
                &quality_score_model,
                rng,
                2,
            )?;
        }
    }

    Ok(())
}

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
