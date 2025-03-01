// This library writes either single ended or paired-ended fastq files.

use anyhow::Result;
use rand::Rng;
use std::io::Write;

use super::fasta_tools::sequence_array_to_string;
use super::file_tools::open_file;
use super::quality_scores::QualityScoreModel;

/// Complement function for DNA nucleotides.
fn complement(nucleotide: u8) -> u8 {
    // 0 = A, 1 = C, 2 = G, 3 = T,
    match nucleotide {
        0 => 3,
        1 => 2,
        2 => 1,
        3 => 0,
        _ => 4,
    }
}

/// Reverse complement function for DNA sequences.
fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    let length = sequence.len();
    let mut rev_comp = Vec::new();
    for i in (0..length).rev() {
        rev_comp.push(complement(sequence[i]))
    }
    rev_comp
}

/// Writes FASTQ files based on the provided dataset.
///
/// # Arguments
///
/// * `fastq_filename` - The prefix for the output FASTQ files.
/// * `overwrite_output` - A boolean flag to enable or disable overwriting
///   existing files.
/// * `paired_ended` - A boolean flag to enable or disable paired-end mode.
/// * `dataset` - A list of `Vec<u8>` representing DNA sequences.
/// * `dataset_order` - A list of indices to order the dataset.
/// * `quality_score_model` - A `QualityScoreModel` object to generate
///   quality scores.
/// * `rng` - A random number generator.
///
/// # Returns
///
/// Returns `()` if successful. Throws an error if there is a problem.
///
/// # Notes
///
/// Currently, only a single R1 file is written, but future versions
/// will support both R1 and R2 files.
pub fn write_fastq<R: Rng>(
    fastq_filename: &str,
    overwrite_output: bool,
    paired_ended: bool,
    dataset: Vec<&Vec<u8>>,
    dataset_order: Vec<usize>,
    quality_score_model: QualityScoreModel,
    rng: &mut R,
) -> Result<()> {
    // The prefix for read names. Reads are numbered in output order.
    let name_prefix = "neat_generated_";
    let mut filename1 = format!("{}_r1.fastq", fastq_filename);
    let mut outfile1 = open_file(&mut filename1, overwrite_output)?;
    // Setting up pairend ended reads. For single ended reads,
    // we have outfile2 = None
    let mut outfile2 = if paired_ended {
        let mut filename2 = format!("{}_r2.fastq", fastq_filename);
        Some(open_file(&mut filename2, overwrite_output)?)
    } else {
        None
    };

    // Write sequences. Ordered index is used for numbering, while
    // read_index is from the shuffled index array from a previous
    // step.
    for (order_index, &read_index) in dataset_order.iter().enumerate() {
        let sequence = dataset[read_index];
        let read_length = sequence.len();

        // Generate quality scores for read1
        let quality_scores = quality_score_model.generate_quality_scores(read_length, rng)?;

        // sequence name
        writeln!(outfile1, "@{}{}/1", name_prefix, order_index + 1)?;
        // Array as a string
        writeln!(outfile1, "{}", sequence_array_to_string(sequence))?;
        // The stupid plus sign
        writeln!(outfile1, "+")?;
        // Qual score of all F's for the whole thing.
        writeln!(outfile1, "{}", quality_scores_to_str(quality_scores))?;

        // Handle paired-end reads
        if let Some(ref mut outfile2) = outfile2 {
            let quality_scores = quality_score_model.generate_quality_scores(read_length, rng)?;

            let rev_comp_sequence = reverse_complement(sequence);
            // sequence name
            writeln!(outfile2, "@{}{}/2", name_prefix, order_index + 1)?;
            // Array as a string
            writeln!(outfile2, "{}", sequence_array_to_string(&rev_comp_sequence))?;
            // The stupid plus sign
            writeln!(outfile2, "+")?;
            // Qual score of all F's for the whole thing.
            writeln!(outfile2, "{}", quality_scores_to_str(quality_scores))?;
        }
    }

    Ok(())
}

/// Converts a vector of quality scores to a string.
fn quality_scores_to_str(array: Vec<u32>) -> String {
    let mut score_text = String::new();
    for score in array {
        score_text += &(((score + 33) as u8) as char).to_string();
    }
    score_text
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::create_rng;
    use std::fs;
    use std::path::Path;

    #[test]
    fn test_complement() {
        let nuc1 = 0;
        let nuc2 = 1;
        let nuc3 = 2;
        let nuc4 = 3;
        let nuc5 = 4;

        assert_eq!(complement(nuc1), 3);
        assert_eq!(complement(nuc2), 2);
        assert_eq!(complement(nuc3), 1);
        assert_eq!(complement(nuc4), 0);
        assert_eq!(complement(nuc5), 4);
    }

    #[test]
    fn test_reverse_complement() {
        let read: Vec<u8> = vec![0, 0, 0, 0, 1, 1, 1, 1];
        let revcomp: Vec<u8> = vec![2, 2, 2, 2, 3, 3, 3, 3];
        assert_eq!(reverse_complement(&read), revcomp);
    }

    #[test]
    fn test_write_fastq_single() -> Result<()> {
        let fastq_filename = "test_single";
        let overwrite_output = true;
        let paired_ended = false;
        let seq1 = vec![0, 0, 0, 0, 1, 1, 1, 1];
        let seq2 = vec![2, 2, 2, 2, 3, 3, 3, 3];
        let mut rng = create_rng(Some("Hello Cruel World"));
        let dataset = vec![&seq1, &seq2];
        let dataset_order = vec![1, 0];
        let quality_score_model = QualityScoreModel::new();
        write_fastq(
            fastq_filename,
            overwrite_output,
            paired_ended,
            dataset,
            dataset_order,
            quality_score_model,
            &mut rng,
        )?;
        let outfile1 = Path::new("test_single_r1.fastq");
        let outfile2 = Path::new("test_single_r2.fastq");
        assert!(outfile1.exists());
        assert!(!outfile2.exists());
        fs::remove_file(outfile1)?;
        Ok(())
    }

    #[test]
    fn test_write_fastq_paired() -> Result<()> {
        let fastq_filename = "test_paired";
        // might as well test the o_o function as well
        let overwrite_output = false;
        let paired_ended = true;
        let seq1 = vec![0, 0, 0, 0, 1, 1, 1, 1];
        let seq2 = vec![2, 2, 2, 2, 3, 3, 3, 3];
        let mut rng = create_rng(Some("Hello Cruel World"));
        let dataset = vec![&seq1, &seq2];
        let dataset_order = vec![1, 0];
        let quality_score_model = QualityScoreModel::new();
        write_fastq(
            fastq_filename,
            overwrite_output,
            paired_ended,
            dataset,
            dataset_order,
            quality_score_model,
            &mut rng,
        )?;
        let outfile1 = Path::new("test_paired_r1.fastq");
        let outfile2 = Path::new("test_paired_r2.fastq");
        assert!(outfile1.exists());
        assert!(outfile2.exists());
        fs::remove_file(outfile1)?;
        fs::remove_file(outfile2)?;
        Ok(())
    }
}
