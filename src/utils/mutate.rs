// This is a basic mutation with SNPs using a basic mutation model.
// mutate_fasta takes a fasta Hashmap and returns a mutated version and the locations of the
// mutations introduced
//
// mutate_sequence adds actual mutations to the fasta sequence

use super::{
    nucleotides::{Nuc, NucModel},
    types::SeqByContig,
};
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
///   is a requested minimum.
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
    file_struct: &SeqByContig,
    minimum_mutations: Option<usize>,
    mutation_rate: Option<f64>,
    rng: &mut R,
) -> Result<(SeqByContig, HashMap<String, Vec<(usize, Nuc, Nuc)>>)> {
    let mut return_struct: SeqByContig = HashMap::new();
    let mut all_variants: HashMap<String, Vec<(usize, Nuc, Nuc)>> = HashMap::new();

    for (name, sequence) in file_struct {
        let sequence_length = sequence.len();
        debug!("Sequence {} is {} bp long", name, sequence_length);
        // Calculate the number of mutation positions based on the mutation rate
        let num_positions = sequence_length as f64 * mutation_rate.unwrap_or(0.0);
        // Ensure the number of mutations meets the minimum threshold, if provided
        let num_positions = minimum_mutations.unwrap_or(0).max(num_positions as usize);
        // Generate mutations to the sequence
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

/// This function takes a vector of Nuc's and mutates a few positions
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
    sequence: &[Nuc],
    num_positions: usize,
    rng: &mut R,
) -> Result<(Vec<Nuc>, Vec<(usize, Nuc, Nuc)>)> {
    debug!("Adding {} mutations", num_positions);
    let mut mutated_record = sequence.to_owned();
    let mut indexes_to_mutate: HashSet<usize> = HashSet::new();
    // choose num_positions distinct indexes to mutate
    while indexes_to_mutate.len() < num_positions {
        let index = rng.random_range(0..sequence.len());
        if mutated_record[index] == Nuc::N {
            continue;
        }
        indexes_to_mutate.insert(index);
    }
    // Build the default mutation model
    // todo incorporate custom models
    let nucleotide_mutation_model = NucModel::new()?;
    // Will hold the variants added to this sequence
    let mut sequence_variants: Vec<(usize, Nuc, Nuc)> = Vec::new();
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
    use crate::{create_rng, utils::nucleotides::random_seq};

    #[test]
    fn test_mutate_sequence() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let seq: Vec<Nuc> = random_seq(&mut rng, 100);
        let num_positions = 1;
        let mutant = mutate_sequence(&seq, num_positions, &mut rng)?;
        assert_eq!(mutant.0.len(), seq.len());
        assert!(!mutant.1.is_empty());
        assert!(seq != mutant.0);
        Ok(())
    }

    #[test]
    fn test_mutate_fasta() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let seq: Vec<Nuc> = random_seq(&mut rng, 100);
        let file_struct: SeqByContig = HashMap::from([("chr1".to_string(), seq.clone())]);
        let (mutated, mutations) = mutate_fasta(&file_struct, Some(1), None, &mut rng)?;
        assert!(mutated.contains_key("chr1"));
        assert!(mutations.contains_key("chr1"));
        let (loc, alt_nuc, ref_nuc) = mutations["chr1"][0];
        assert_eq!(ref_nuc, seq[loc]);
        assert_ne!(alt_nuc, ref_nuc);
        Ok(())
    }

    #[test]
    fn test_mutate_fasta_no_mutations() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let seq: Vec<Nuc> = random_seq(&mut rng, 100);
        let file_struct: SeqByContig = HashMap::from([("chr1".to_string(), seq.clone())]);
        let mut rng = create_rng(Some("Hello Cruel World"));
        let (mutated, mutations) = mutate_fasta(&file_struct, None, None, &mut rng)?;
        assert!(mutated.contains_key("chr1"));
        assert!(mutations.contains_key("chr1"));
        assert!(mutations["chr1"].is_empty());
        assert_eq!(seq, mutated["chr1"]);
        Ok(())
    }
}
