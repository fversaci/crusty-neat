// This is a basic mutation with SNPs using a basic mutation model.
// mutate_fasta takes a fasta Hashmap and returns a mutated version and the locations of the
// mutations introduced
//
// mutate_sequence adds actual mutations to the fasta sequence

use super::nucleotides::NucModel;
use anyhow::{anyhow, Result};
use log::{debug, error};
use rand::Rng;
use std::collections::{HashMap, HashSet};

/// Mutates a fasta file
///
/// # Arguments
///
/// * `file_struct` - A hashmap of contig names (keys) and a vector
///   representing the original sequence.
/// * `minimum_mutations` - is a usize or None that indicates if there
///   is a requested minimum. The default is for rusty-neat to allow 0 mutations.
/// * `rng` - random number generator for the run
///
/// # Returns
///
/// A tuple with pointers to:
/// * A hashmap with keys that are contig names and a vector with the
///   mutated sequence
/// * A vector of tuples containing the location and alt of each
///   variant
///
/// This function performs a basic calculation (length x mutation rate
/// +/- a random amount) and chooses that many positions along the
/// sequence to mutate. It then builds a return string that represents
/// the altered sequence and stores all the variants.
pub fn mutate_fasta<R: Rng>(
    file_struct: &HashMap<String, Vec<u8>>,
    minimum_mutations: Option<usize>,
    rng: &mut R,
) -> Result<(
    HashMap<String, Vec<u8>>,
    HashMap<String, Vec<(usize, u8, u8)>>,
)> {
    const MUT_RATE: f64 = 0.01;
    let mut return_struct: HashMap<String, Vec<u8>> = HashMap::new();
    let mut all_variants: HashMap<String, Vec<(usize, u8, u8)>> = HashMap::new();

    for (name, sequence) in file_struct {
        let sequence_length = sequence.len();
        debug!("Sequence {} is {} bp long", name, sequence_length);
        // Calculate the number of mutation positions based on
        // mutation rate
        let mut num_positions: f64 = sequence_length as f64 * MUT_RATE;
        // Introduce random fluctuation to the mutation count
        let factor: f64 = rng.random_range(0.0..0.10);
        let sign: f64 = if rng.random_bool(0.25) { -1.0 } else { 1.0 };
        num_positions += sign * factor;
        // Ensure the number of mutations meets the minimum threshold if provided
        let num_positions = minimum_mutations.unwrap_or(0).max(num_positions as usize);
        // Apply mutations to the sequence
        let (mutated_record, contig_mutations) = mutate_sequence(sequence, num_positions, rng)?;
        // Store the mutated sequence and corresponding mutations
        return_struct
            .entry(name.clone())
            .or_insert(mutated_record.clone());
        all_variants.entry(name.clone()).or_insert(contig_mutations);
    }
    // Return mutated sequences and their mutation details
    Ok((return_struct, all_variants))
}

/// This function takes a vector of u8's and mutates a few positions
/// at random. It returns the mutated sequence and a list of tuples
/// with the position and the alts of the SNPs.
///
/// # Arguments
///
/// * `sequence` - A vector representing the original sequence.
/// * `num_positions` - The number of mutations to add to this
///   sequence.
/// * `rng` - random number generator for the run
///
/// # Returns
///
/// A tuple with pointers to:
/// * A vector representing the mutated sequence
/// * A vector of tuples containing the location, alt and ref of each
///   variant
fn mutate_sequence<R: Rng>(
    sequence: &[u8],
    num_positions: usize,
    rng: &mut R,
) -> Result<(Vec<u8>, Vec<(usize, u8, u8)>)> {
    debug!("Adding {} mutations", num_positions);
    let mut mutated_record = sequence.to_owned();
    let mut indexes_to_mutate: HashSet<usize> = HashSet::new();
    // choose num_positions distinct indexes to mutate
    while indexes_to_mutate.len() < num_positions {
        let index = rng.random_range(0..sequence.len());
        if mutated_record[index] == 4 {
            continue;
        }
        indexes_to_mutate.insert(index);
    }
    // Build the default mutation model
    // todo incorporate custom models
    let nucleotide_mutation_model = NucModel::new()?;
    // Will hold the variants added to this sequence
    let mut sequence_variants: Vec<(usize, u8, u8)> = Vec::new();
    // for each index, picks a new base
    for index in indexes_to_mutate {
        // remember the reference for later.
        let reference_base = sequence[index];
        // pick a new base and assign the position to it.
        mutated_record[index] = nucleotide_mutation_model.choose_new_nuc(reference_base, rng)?;
        // This check, only run in debugging mode, ensures that the
        // model actually mutated the base.
        if cfg!(debug_assertions) && mutated_record[index] == reference_base {
            error!("Need to check the code choosing nucleotides");
            return Err(anyhow!(
                "BUG: Mutation model failed to mutate the base. This should not happen."
            ));
        }
        // add the location and alt base for the variant
        sequence_variants.push((index, mutated_record[index], reference_base))
    }

    Ok((mutated_record, sequence_variants))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::create_rng;

    #[test]
    fn test_mutate_sequence() -> Result<()> {
        let seq1: Vec<u8> = vec![4, 4, 0, 0, 0, 1, 1, 2, 0, 3, 1, 1, 1];
        let num_positions = 2;
        let mut rng = create_rng(Some("Hello Cruel World"));
        let mutant = mutate_sequence(&seq1, num_positions, &mut rng)?;
        assert_eq!(mutant.0.len(), seq1.len());
        assert!(!mutant.1.is_empty());
        assert_eq!(mutant.0[0], 4);
        assert_eq!(mutant.0[1], 4);
        Ok(())
    }

    #[test]
    fn test_mutate_fasta() -> Result<()> {
        let seq = vec![4, 4, 0, 0, 0, 1, 1, 2, 0, 3, 1, 1, 1];
        let file_struct: HashMap<String, Vec<u8>> =
            HashMap::from([("chr1".to_string(), seq.clone())]);
        let mut rng = create_rng(Some("Hello Cruel World"));
        let mutations = mutate_fasta(&file_struct, Some(1), &mut rng)?;
        assert!(mutations.0.contains_key("chr1"));
        assert!(mutations.1.contains_key("chr1"));
        let mutation_location = mutations.1["chr1"][0].0;
        let mutation_alt = mutations.1["chr1"][0].1;
        let mutation_ref = mutations.1["chr1"][0].2;
        assert_eq!(mutation_ref, seq[mutation_location]);
        assert_ne!(mutation_alt, mutation_ref);
        Ok(())
    }

    #[test]
    fn test_mutate_fasta_no_mutations() -> Result<()> {
        let seq = vec![4, 4, 0, 0, 0, 1, 1, 2, 0, 3, 1, 1, 1];
        let file_struct: HashMap<String, Vec<u8>> =
            HashMap::from([("chr1".to_string(), seq.clone())]);
        // if a random mutation suddenly pops up in a build, it's probably the seed for this.
        let mut rng = create_rng(Some("Hello Cruel World"));
        let mutations = mutate_fasta(&file_struct, None, &mut rng)?;
        assert!(mutations.0.contains_key("chr1"));
        assert!(mutations.1.contains_key("chr1"));
        assert!(mutations.1["chr1"].is_empty());
        Ok(())
    }
}
