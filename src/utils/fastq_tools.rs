// This library writes either single ended or paired-ended fastq files.

use crate::utils::file_tools::open_file;
use crate::utils::mutation;
use crate::utils::mutation::Mutation;
use crate::utils::nucleotides::{Nuc, reverse_complement, seq_to_string};
use crate::utils::quality_model::QualityModel;
use anyhow::Result;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rayon::iter::{IndexedParallelIterator, ParallelIterator};
use rayon::prelude::ParallelSlice;
use std::fs::create_dir_all;
use std::io::Write;
use std::path::Path;

fn gen_read_mutations<R: Rng>(
    seq: &[Nuc],
    quality_nums: &[usize],
    quality_model: &QualityModel,
    rng: &mut R,
) -> Result<Vec<Mutation>> {
    let mut_model = &quality_model.seq_err_model.mut_model;
    let mut mutations: Vec<Mutation> = Vec::new();
    for (pos, q) in quality_nums.iter().enumerate() {
        let p = quality_model.get_seq_err_prob(*q)?;
        if rng.random_bool(p) {
            let mut_type = mut_model.get_mut_type(rng)?;
            let mutation = mut_model.create_mutation(mut_type, seq, pos, rng);
            // if mutation is ok, add it to the mutations list
            if let Ok(mutation) = mutation {
                mutations.push(mutation);
            }
        }
    }
    mutation::sort_filter_overlap(&mut mutations);
    Ok(mutations)
}

/// apply mutations to the read
///
/// # Arguments
///
/// * `orig_seq` - The original sequence.
/// * `mut_seq` - The mutated sequence.
/// * `orig_quality_nums` - The original quality scores.
/// * `mut_quality_nums` - The mutated quality scores.
/// * `mutations` - The mutations to be applied.
pub fn apply_read_mutations(
    orig_seq: &[Nuc],
    mut_seq: &mut Vec<Nuc>,
    orig_qnums: &[usize],
    mut_qnums: &mut Vec<usize>,
    mutations: &[Mutation],
) {
    let mut prev = 0;
    for mutation in mutations {
        let (start, end) = mutation.get_skip_range();
        mut_seq.extend(&orig_seq[prev..start]);
        mut_qnums.extend(&orig_qnums[prev..start]);
        prev = end;
        match mutation {
            Mutation::Snp { alt_base, .. } => {
                mut_seq.push(*alt_base);
                mut_qnums.push(orig_qnums[start]);
            }
            Mutation::Ins { alt_bases, .. } => {
                mut_seq.extend(alt_bases);
                // repeat same quality score for the inserted bases
                mut_qnums.extend(std::iter::repeat_n(orig_qnums[start], alt_bases.len()));
            }
            Mutation::Del { .. } => {
                // do nothing
            }
        }
    }
    mut_seq.extend(&orig_seq[prev..]);
    mut_qnums.extend(&orig_qnums[prev..]);
}

macro_rules! qs_char {
    ($q:expr) => {
        ($q + 33) as u8 as char
    };
}

fn write_sequence_with_errors<W: Write, R: Rng>(
    outfile: &mut W,
    name_prefix: &str,
    order_index: usize,
    sequence: &[Nuc],
    quality_model: &QualityModel,
    rng: &mut R,
    read_number: u8,
) -> Result<()> {
    let read_length = sequence.len();
    let quality_nums = quality_model.generate_quality_scores(read_length, rng)?;
    let mutations = gen_read_mutations(sequence, &quality_nums, quality_model, rng)?;
    let mut mut_sequence = Vec::new();
    let mut mut_qnums = Vec::new();
    apply_read_mutations(
        sequence,
        &mut mut_sequence,
        &quality_nums,
        &mut mut_qnums,
        &mutations,
    );

    let quality_scores: String = mut_qnums.iter().map(|x| qs_char!(x)).collect();
    write_sequence(
        outfile,
        name_prefix,
        order_index,
        &mut_sequence,
        &quality_scores,
        read_number,
    )
}

fn write_sequence<W: Write>(
    outfile: &mut W,
    name_prefix: &str,
    order_index: usize,
    sequence: &[Nuc],
    quality_scores: &str,
    read_number: u8,
) -> Result<()> {
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
pub fn write_fastq<R: Rng + Send + Sync>(
    output_prefix: &Path,
    overwrite_output: bool,
    paired_ended: bool,
    reads: Vec<&[Nuc]>,
    quality_model: &QualityModel,
    rng: &mut R,
) -> Result<()> {
    let dir1 = output_prefix.with_extension("reads_R1");
    let dir2 = output_prefix.with_extension("reads_R2");

    create_dir_all(&dir1)?;
    if paired_ended {
        create_dir_all(&dir2)?;
    }
    let seed = rng.next_u64();
    let max_chunk = 100000;
    reads
        .par_chunks(max_chunk)
        .enumerate()
        .try_for_each(|(chunk_index, chunk)| -> Result<()> {
            if chunk.is_empty() {
                return Ok(());
            }

            // create a local rng for this chunk
            let mut local_rng = StdRng::seed_from_u64(seed + (chunk_index as u64));

            let name_prefix = format!("crusty_neat_generated_{}_", chunk_index);
            let file1_path = dir1.join(format!("{}.fastq.gz", chunk_index));
            let mut file1 = open_file(&file1_path, overwrite_output)?;

            for (order_index, sequence) in chunk.iter().enumerate() {
                write_sequence_with_errors(
                    &mut file1,
                    &name_prefix,
                    order_index,
                    sequence,
                    quality_model,
                    &mut local_rng,
                    1,
                )?;
            }

            if paired_ended {
                let file2_path = dir2.join(format!("{}.fastq.gz", chunk_index));
                let mut file2 = open_file(&file2_path, overwrite_output)?;
                for (order_index, sequence) in chunk.iter().enumerate() {
                    let rev_comp_sequence = reverse_complement(sequence);
                    write_sequence_with_errors(
                        &mut file2,
                        &name_prefix,
                        order_index,
                        &rev_comp_sequence,
                        quality_model,
                        &mut local_rng,
                        2,
                    )?;
                }
            }

            Ok(())
        })
}
