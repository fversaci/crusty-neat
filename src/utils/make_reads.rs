// This is the core functionality of NEAT. Generate reads turns a
// mutated fasta array into short reads.
// The idea of cover_dataset is we generate a set of coordinates
// That define reads and covers the dataset coverage number times, to
// give the contig the proper read coverage. Generate reads uses this
// to create a list of coordinates to take slices from the mutated
// fasta file. These will either be read-length fragments or fragment
// model length fragments.

use crate::utils::distributions::IntDistribution;
use crate::utils::nucleotides::Nuc;
use anyhow::{Result, anyhow};
use rand::Rng;
use rand_distr::Distribution;

/// Covers a sequence by selecting random reads. Starts at the
/// beginning, advances by one read length, then jumps forward by a
/// random amount (0 to half the read length) to select the next
/// position. Repeats until the read tail exceeds the sequence length,
/// then restarts by wrapping around.
///
/// # Arguments
///
/// - `span_length`: Total number of bases in the sequence
/// - `read_length`: The length of the reads for this run
/// - `fragment_distribution`: The distribution to generate fragment
///   lengths from
/// - `coverage`: The coverage depth for the reads
/// - `rng`: The random number generator for the run
///
/// # Returns
///
/// A vector of tuples (usize, usize), denoting the start and end
/// positions of the fragments of DNA that was sequenced.
fn cover_dataset<R: Rng>(
    span_length: usize,
    read_length: usize,
    fragment_distribution: IntDistribution,
    coverage: usize,
    paired_ended: bool,
    rng: &mut R,
) -> Vec<(usize, usize)> {
    // Preallocate frag_set to avoid excessive reallocations
    let mut frag_set = Vec::with_capacity(coverage * span_length / read_length);
    let target_cov_bases = coverage * span_length;
    let mut curr_cov_bases = 0;

    let mut frag_start = 0;
    // Iterate until required coverage is achieved
    while curr_cov_bases < target_cov_bases {
        let fragment_length = fragment_distribution.sample(rng) as usize;

        let frag_end = span_length.min(frag_start + fragment_length);
        frag_set.push((frag_start, frag_end));
        if paired_ended {
            // paired ended: read at most twice the read_lenght from
            // each fragment
            curr_cov_bases += (2 * read_length).min(frag_end - frag_start);
        } else {
            // single ended: read at most read_lenght from each fragment
            curr_cov_bases += read_length.min(frag_end - frag_start);
        }

        // Randomized offset to avoid uniform patterns
        let wildcard = rng.random_range(0..(read_length / 4));
        // adds to the frag_start to give it some spice
        frag_start = (frag_end + wildcard) % span_length;
    }

    frag_set
}

/// This takes a mutated sequence and produces a set of reads based on
/// the mutated sequence. For paired ended reads, this will generate a
/// set of reads from each end, by taking the reverse complement int
/// the output
///
/// # Arguments
///
/// `mutated_sequence`: a vector of Nuc's representing the mutated sequence.
/// `read_length`: the length ef the reads for this run
/// `coverage`: the average depth of coverage for this run
/// `paired_ended`: is the run paired ended or not
/// `mean`: the mean of the fragment distribution
/// `st_dev`: the standard deviation of the fragment distribution
/// `rng`: the random number generator for the run
///
/// # Returns
///
/// Vec of vectors representing the read sequences.
pub fn generate_reads<'a, R: Rng>(
    mutated_sequence: &'a [Nuc],
    read_length: &usize,
    coverage: &usize,
    paired_ended: bool,
    mean: Option<f64>,
    st_dev: Option<f64>,
    rng: &mut R,
) -> Result<Vec<&'a [Nuc]>> {
    let fragment_distribution = if paired_ended {
        IntDistribution::new_normal(
            mean.ok_or_else(|| anyhow!("mean can't be None when paired_ended is true"))?,
            st_dev.ok_or_else(|| anyhow!("std_dev can't be None when paired_ended is true"))?,
        )?
    } else {
        IntDistribution::new_constant(*read_length as i64)
    };

    // set up some defaults and storage
    let mut read_set: Vec<&[Nuc]> = Vec::new();
    // length of the mutated sequence
    let seq_len = mutated_sequence.len();
    // Generate a vector of read positions
    let read_positions: Vec<(usize, usize)> = cover_dataset(
        seq_len,
        *read_length,
        fragment_distribution,
        *coverage,
        paired_ended,
        rng,
    );
    // Generate the reads from the read positions.
    for (start, end) in read_positions {
        read_set.push(mutated_sequence[start..end].into());
    }
    // puts the reads in the heap.
    if read_set.is_empty() {
        Err(anyhow!("No reads generated"))
    } else {
        Ok(read_set)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{create_rng, utils::nucleotides::random_seq};

    #[test]
    fn test_cover_dataset() {
        let span_length = 100;
        let read_length = 10;
        let fragment_distribution = IntDistribution::new_constant(10);
        let coverage = 1;
        let paired_ended = true;
        let mut rng = create_rng(Some("Hello Cruel World"));

        let cover = cover_dataset(
            span_length,
            read_length,
            fragment_distribution,
            coverage,
            paired_ended,
            &mut rng,
        );
        assert_eq!(cover[0], (0, 10))
    }

    #[test]
    fn test_gap_function() {
        let span_length = 100_000;
        let read_length = 100;
        let fragment_distribution = IntDistribution::new_constant(300);
        let coverage = 1;
        let paired_ended = true;
        let mut rng = create_rng(Some("Hello Cruel World"));

        let cover = cover_dataset(
            span_length,
            read_length,
            fragment_distribution,
            coverage,
            paired_ended,
            &mut rng,
        );
        assert_eq!(cover[0], (0, 300))
    }

    #[test]
    fn test_generate_reads_single() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let mutated_sequence = random_seq(&mut rng, 100);
        let read_length = 10;
        let coverage = 1;
        let paired_ended = false;
        let mean = None;
        let st_dev = None;
        let mut rng = create_rng(Some("Hello Cruel World"));
        let reads = generate_reads(
            &mutated_sequence,
            &read_length,
            &coverage,
            paired_ended,
            mean,
            st_dev,
            &mut rng,
        )?;
        println!("{:?}", reads);
        assert!(reads.contains(&&mutated_sequence[0..read_length]));
        Ok(())
    }

    #[test]
    fn test_seed_rng() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let mutated_sequence = random_seq(&mut rng, 100);
        let read_length = 10;
        let coverage = 1;
        let paired_ended = false;
        let mean = None;
        let st_dev = None;
        let mut rng = create_rng(Some("Hello Cruel World"));
        let run1 = generate_reads(
            &mutated_sequence,
            &read_length,
            &coverage,
            paired_ended,
            mean,
            st_dev,
            &mut rng,
        )?;

        let mut rng = create_rng(Some("Hello Cruel World"));
        let run2 = generate_reads(
            &mutated_sequence,
            &read_length,
            &coverage,
            paired_ended,
            mean,
            st_dev,
            &mut rng,
        )?;

        assert_eq!(run1, run2);
        Ok(())
    }

    #[test]
    fn test_generate_reads_paired() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let mutated_sequence = random_seq(&mut rng, 1000);
        let read_length = 100;
        let coverage = 1;
        let paired_ended = true;
        let mean = Some(200.0);
        let st_dev = Some(1.0);
        let mut rng = create_rng(Some("Hello Cruel World"));
        let reads = generate_reads(
            &mutated_sequence,
            &read_length,
            &coverage,
            paired_ended,
            mean,
            st_dev,
            &mut rng,
        )?;
        assert!(!reads.is_empty());
        Ok(())
    }
}
