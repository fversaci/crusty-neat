use crate::utils::mutation;

use super::{
    mutation::{Mutation, MutationType},
    mutation_model::Region,
    nucleotides::Nuc,
    ref_mutation_model::RefMutationModel,
    types::{MutByContig, SeqByContig},
};
use anyhow::{anyhow, Result};
use log::debug;
use rand::Rng;
use rand_distr::{Binomial, Distribution};
use std::collections::HashMap;
use strum::IntoEnumIterator;

/// Mutates a genome by generating mutations for each contig.
///
/// # Arguments
///
/// * `genome` - A hashmap of contig names (keys) and a vector
///   representing the original sequence.
/// * `mutation_model` - The mutation model for the run
/// * `rng` - random number generator for the run
///
/// # Returns
///
/// A tuple with pointers to:
/// * A vector of mutations for each contig
pub fn mutate_genome<R: Rng>(
    genome: &SeqByContig,
    mutation_model: &RefMutationModel,
    rng: &mut R,
) -> Result<MutByContig> {
    let mut mutations: MutByContig = HashMap::new();
    let mut_rates = mutation_model.mut_rates.clone().unwrap();
    for (contig, regions) in mut_rates {
        let sequence = &genome[&contig];
        for region in regions {
            debug!(
                "Segment of {}: from {} to {}.",
                contig, region.start, region.end
            );
            // Generate mutations for the region
            let contig_mutations = mutate_region(sequence, mutation_model, &region, rng)?;
            // Add the generated mutations to the list of mutations for this contig
            mutations
                .entry(contig.clone())
                .or_default()
                .extend(contig_mutations);
        }
    }
    Ok(mutations)
}

/// This function takes a slice of Nuc's and mutates num_mutations at
/// random. It returns the vector of mutations.
///
/// # Arguments
///
/// * `sequence` - A vector representing the original sequence.
/// * `num_mutations` - The number of mutations to add to this
///   sequence.
/// * `rng` - random number generator for the run
///
/// # Returns
///
/// A tuple with pointers to:
/// * A vector of mutations for this region
fn mutate_region<R: Rng>(
    sequence: &[Nuc],
    mutation_model: &RefMutationModel,
    region: &Region,
    rng: &mut R,
) -> Result<Vec<Mutation>> {
    // calculate the number of mutations (for each mutation) using a
    // Binomial distribution
    let region_len = region.end - region.start;

    // generate mutations
    let mut mutations: Vec<Mutation> = Vec::new();
    for m in MutationType::iter() {
        let m_rate = region.rate * mutation_model.mm.get_mut_probability(m);
        let bin = Binomial::new(region_len as u64, m_rate)?;
        let num_mut: usize = bin.sample(rng).try_into()?;
        debug!("Adding {} {} mutations", num_mut, m);
        let positions: Vec<usize> = (0..num_mut)
            .map(|_| rng.random_range(region.start..region.end))
            .collect();
        mutations.extend(
            positions
                .into_iter()
                .filter_map(|pos| mutation_model.create_mutation(m, sequence, pos, rng).ok()),
        );
    }
    // sort and filter out incompatible mutations
    debug!("Total of {} potential mutations", mutations.len());
    mutation::sort_filter_overlap(&mut mutations);
    debug!(
        "{} mutations remain, after removing overlapping ones",
        mutations.len()
    );

    Ok(mutations)
}

/// Applies the mutations to the reference genome and returns the new genome.
pub fn apply_mutations(genome: &SeqByContig, mutations: &MutByContig) -> Result<SeqByContig> {
    let mut new_genome = genome.clone();
    for (contig, mutations) in mutations {
        let sequence = new_genome.get_mut(contig).unwrap();
        for mutation in mutations {
            match mutation {
                Mutation::Snp { pos, alt_base, .. } => {
                    sequence[*pos] = *alt_base;
                }
                _ => return Err(anyhow!("Non-SNP mutations not supported yet")),
            }
        }
    }
    Ok(new_genome)
}

#[cfg(test)]
mod tests {
    use anyhow::anyhow;

    use super::*;
    use crate::{create_rng, utils::nucleotides::random_seq};

    #[test]
    fn test_mutate_genome() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let seq: Vec<Nuc> = random_seq(&mut rng, 100);
        let genome: SeqByContig = HashMap::from([("chr1".to_string(), seq.clone())]);
        let mut_model = RefMutationModel::new_all_contigs(&genome, 0.1)?;
        let mutations = mutate_genome(&genome, &mut_model, &mut rng)?;
        assert!(mutations.contains_key("chr1"));
        let mutation = mutations["chr1"][0].clone();
        match mutation {
            Mutation::Snp {
                pos,
                ref_base,
                alt_base,
            } => {
                assert_eq!(seq[pos], ref_base);
                assert_ne!(alt_base, ref_base);
            }
            _ => return Err(anyhow!("Expected a SNP mutation")),
        }
        Ok(())
    }

    #[test]
    fn test_mutate_genome_no_mutations() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let seq: Vec<Nuc> = random_seq(&mut rng, 100);
        let genome: SeqByContig = HashMap::from([("chr1".to_string(), seq.clone())]);
        let mut_model = RefMutationModel::new_all_contigs(&genome, 0.0)?;
        let mutations = mutate_genome(&genome, &mut_model, &mut rng)?;
        assert!(mutations.contains_key("chr1"));
        assert!(mutations["chr1"].is_empty());
        Ok(())
    }
}
