use anyhow::{anyhow, Result};
// This is the core functionality of NEAT. Generate reads turns a
// mutated fasta array into short reads.
// The idea of cover_dataset is we generate a set of coordinates
// That define reads and covers the dataset coverage number times, to
// give the contig the proper read coverage. Generate reads uses this
// to create a list of coordinates to take slices from the mutated
// fasta file. These will either be read-length fragments or fragment
// model length fragments.

use rand::seq::SliceRandom;
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::collections::{HashSet, VecDeque};

/// This function selects the positions of the reads. It starts at the
/// beginning and goes out one read length, then picks a random jump
/// between 0 and half the read length to move And picks those
/// coordinates for a second read. Once the tail of the read is past
/// the end, we start over again at 0.
///
/// # Arguments
///
/// - `span_length`: Total number of bases in the sequence
/// - `read_length`: The length of the reads for this run
/// - `fragment_pool`: a vector of sizes for the fragments. If empty,
///   it will instead be filled by the read_length (single ended
///   reads) paired_ended: true or false if the run is paired ended
///   mode or not.
/// - `coverage`: The coverage depth for the reads
///
/// # Returns
///
/// A vector of tuples (usize, usize), denoting the start and end
/// positions of the fragment of DNA that was sequenced.
fn cover_dataset<R: Rng>(
    span_length: usize,
    read_length: usize,
    mut fragment_pool: Vec<usize>,
    coverage: usize,
    rng: &mut R,
) -> Vec<(usize, usize)> {
    // Reads that will be start and end of the fragment.
    let mut read_set: Vec<(usize, usize)> = vec![];
    let mut cover_fragment_pool: VecDeque<usize>;
    if fragment_pool.is_empty() {
        // set the shuffled fragment pool just equal to an instance of read_length
        cover_fragment_pool = VecDeque::from([read_length]);
    } else {
        // shuffle the fragment pool
        fragment_pool.shuffle(rng);
        cover_fragment_pool = VecDeque::from(fragment_pool)
    }
    // Gap size to keep track of how many uncovered bases we have per
    // layer, to help decide if we need more layers
    let mut gap_size: usize = 0;
    let mut layer_count: usize = 0;
    // start this party off at zero.
    let mut start: usize = 0;
    // create coverage number of layers
    while layer_count <= coverage {
        let fragment_length = cover_fragment_pool[0];
        cover_fragment_pool.push_back(fragment_length);
        let temp_end = start + fragment_length;
        if temp_end > span_length {
            // TODO some variation on this modulo idea will work for
            // bacterial reads
            start = temp_end % span_length;
            gap_size += start;
            //
            if gap_size >= span_length {
                // if we have accumulated enough gap, then we need to
                // run the same layer again.  We'll reset gap size but
                // not increment layer_count.
                gap_size %= span_length;
                continue;
            } else {
                layer_count += 1;
                continue;
            }
        }
        read_set.push((start, temp_end));
        // insert size is the number of bases between reads in the
        // fragment for paired ended reads if these are singled ended
        // reads, then the insert size will always be -read_length
        if fragment_length > (read_length * 2) {
            // if there's any insert size on paired ended reads, we'll add
            // that to the gap to ensure adequate coverage.
            gap_size += fragment_length - (read_length * 2)
        };
        // Picks a number between zero and a quarter of a read length
        let wildcard: usize = rng.random_range(0..(read_length / 4));
        // adds to the start to give it some spice
        start += temp_end + wildcard;
        // sanity check. If we are already out of bounds, take the modulo
        if start >= span_length {
            // get us back in bounds
            start %= span_length;
            // add the gap
            gap_size += start;
        } else {
            // still in bounds, just add the gap
            gap_size += wildcard;
        }
    }
    read_set
}

/// This takes a mutated sequence and produces a set of reads based on
/// the mutated sequence. For paired ended reads, this will generate a
/// set of reads from each end, by taking the reverse complement int
/// the output
///
/// # Arguments
///
/// `mutated_sequence`: a vector of u8's representing the mutated sequence.
/// `read_length`: the length ef the reads for this run
/// `coverage`: the average depth of coverage for this run
/// `rng`: the random number generator for the run
///
/// # Returns
///
/// HashSet of vectors representing the read sequences.
pub fn generate_reads<R: Rng>(
    mutated_sequence: &[u8],
    read_length: &usize,
    coverage: &usize,
    paired_ended: bool,
    mean: Option<f64>,
    st_dev: Option<f64>,
    rng: &mut R,
) -> Result<HashSet<Vec<u8>>> {
    let mut fragment_pool: Vec<usize> = Vec::new();
    if paired_ended {
        let num_frags = (mutated_sequence.len() / read_length) * (coverage * 2);
        let fragment_distribution = Normal::new(
            mean.ok_or_else(|| anyhow!("mean can't be None when paired_ended is true"))?,
            st_dev.ok_or_else(|| anyhow!("std_dev can't be None when paired_ended is true"))?,
        )?;
        // Add fragments to the fragment pool
        fragment_pool
            .extend((0..num_frags).map(|_| fragment_distribution.sample(rng).round() as usize));
    }
    // set up some defaults and storage
    let mut read_set: HashSet<Vec<u8>> = HashSet::new();
    // length of the mutated sequence
    let seq_len = mutated_sequence.len();
    // Generate a vector of read positions
    let read_positions: Vec<(usize, usize)> =
        cover_dataset(seq_len, *read_length, fragment_pool, *coverage, rng);
    // Generate the reads from the read positions.
    for (start, end) in read_positions {
        read_set.insert(mutated_sequence[start..end].into());
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
    use crate::create_rng;

    #[test]
    fn test_cover_dataset() {
        let span_length = 100;
        let read_length = 10;
        let fragment_pool = vec![10];
        let coverage = 1;
        let mut rng = create_rng(Some("Hello Cruel World"));

        let cover = cover_dataset(span_length, read_length, fragment_pool, coverage, &mut rng);
        assert_eq!(cover[0], (0, 10))
    }

    #[test]
    fn test_gap_function() {
        let span_length = 100_000;
        let read_length = 100;
        let fragment_pool = vec![300];
        let coverage = 1;
        let mut rng = create_rng(Some("Hello Cruel World"));

        let cover = cover_dataset(span_length, read_length, fragment_pool, coverage, &mut rng);
        assert_eq!(cover[0], (0, 300))
    }

    #[test]
    fn test_generate_reads_single() -> Result<()> {
        let mutated_sequence = vec![0, 0, 1, 0, 3, 3, 3, 3, 0, 0, 0, 0, 0, 2, 2, 2, 4, 4, 4, 4];
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
        assert!(reads.contains(&(vec![0, 0, 1, 0, 3, 3, 3, 3, 0, 0])));
        Ok(())
    }

    #[test]
    fn test_seed_rng() -> Result<()> {
        let mutated_sequence = vec![0, 0, 1, 0, 3, 3, 3, 3, 0, 0, 0, 0, 0, 2, 2, 2, 4, 4, 4, 4];
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
        let mutated_sequence: Vec<u8> = std::iter::repeat(1).take(100_000).collect();
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
