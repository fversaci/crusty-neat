use super::{
    mutation::Mutation,
    nucleotides::{Nuc, NucModel},
    types::SeqByContig,
};
use anyhow::Result;
use log::debug;
use rand::Rng;
use std::collections::{HashMap, HashSet};

/// Mutates a fasta hashmap of contigs and sequences.
///
/// # Arguments
///
/// * `file_struct` - A hashmap of contig names (keys) and a vector
///   representing the original sequence.
/// * `minimum_mutations` - is a usize or None that indicates if there
///   is a requested minimum.
/// * `mutation_rate` - is a f64 or None that indicates the mutation
///   rate to use.
/// * `rng` - random number generator for the run
///
/// # Returns
///
/// A tuple with pointers to:
/// * A hashmap with keys that are contig names and a vector with the
///   mutated sequence
/// * A vector of mutations for each contig
pub fn mutate_fasta<'a, R: Rng>(
    file_struct: &'a SeqByContig,
    minimum_mutations: Option<usize>,
    mutation_rate: Option<f64>,
    rng: &mut R,
) -> Result<(SeqByContig, HashMap<String, Vec<Mutation<'a>>>)> {
    let mut return_struct: SeqByContig = HashMap::new();
    let mut all_variants: HashMap<String, Vec<Mutation<'a>>> = HashMap::new();

    for (name, sequence) in file_struct {
        let sequence_length = sequence.len();
        debug!("Sequence {} is {} bp long", name, sequence_length);
        // Calculate the number of mutation positions based on the mutation rate
        let num_positions = sequence_length as f64 * mutation_rate.unwrap_or(0.0);
        // Ensure the number of mutations meets the minimum threshold, if provided
        let num_positions = minimum_mutations.unwrap_or(0).max(num_positions as usize);
        // Generate mutations for the sequence
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

/// This function takes a vector of Nuc's and mutates num_positions at
/// random. It returns the mutated sequence and a vector of mutations.
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
/// * A vector containing the mutated sequence
/// * A vector of mutations for this sequence
fn mutate_sequence<'a, R: Rng>(
    sequence: &'a [Nuc],
    num_positions: usize,
    rng: &mut R,
) -> Result<(Vec<Nuc>, Vec<Mutation<'a>>)> {
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
    let mut sequence_variants: Vec<Mutation<'a>> = Vec::new();
    // for each index, picks a new base
    for pos in indexes_to_mutate {
        // remember the reference for later.
        let reference_base = sequence[pos];
        // pick a new base and assign the position to it.
        mutated_record[pos] = nucleotide_mutation_model.choose_new_nuc(reference_base, rng)?;
        // construct the snp mutation
        let snp = Mutation::new_snp(sequence, pos, mutated_record[pos])?;
        sequence_variants.push(snp);
    }

    Ok((mutated_record, sequence_variants))
}

#[cfg(test)]
mod tests {
    use anyhow::anyhow;

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
        let mutation = mutations["chr1"][0].clone();
        match mutation {
            Mutation::Snp {
                pos,
                ref_base,
                alt_base,
            } => {
                assert_eq!(seq[pos], *ref_base);
                assert_ne!(alt_base, *ref_base);
            }
            _ => return Err(anyhow!("Expected a SNP mutation")),
        }
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
