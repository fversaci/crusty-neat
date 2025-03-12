use crate::utils::nucleotides::{Nuc, seq_to_string};
use anyhow::{Result, anyhow};
use serde::{Deserialize, Serialize};
use strum_macros::{Display, EnumIter, EnumString};

/// A genomic mutation (SNP, insertion, deletion)
#[derive(Debug, Clone)]
pub enum Mutation {
    Snp {
        pos: usize,
        ref_base: Nuc,
        alt_base: Nuc,
    },
    Ins {
        pos: usize,
        ref_base: Nuc,
        alt_bases: Vec<Nuc>, // ref_base + inserted sequence
    },
    Del {
        pos: usize,
        ref_bases: Vec<Nuc>, // one base + deleted sequence
    },
}

/// The type of mutation
#[derive(
    Debug, Clone, Copy, Eq, PartialEq, Hash, Serialize, Deserialize, EnumIter, EnumString, Display,
)]
pub enum MutationType {
    Snp,
    Ins,
    Del,
}

impl Mutation {
    /// Get the range of the reference sequence that should not be
    /// affected by other mutations
    pub fn get_interference_range(&self) -> (usize, usize) {
        match self {
            // trinucleotide-context SNPs depend on three positions
            Mutation::Snp { pos, .. } => (pos - 1, pos + 2),
            // Keep the nucleotide before the insertion
            Mutation::Ins { pos, .. } => (*pos, *pos + 1),
            // Keep the nucleotide before the deletion and don't
            // mutate anything else within the deleted subsequence
            Mutation::Del { pos, ref_bases } => (*pos, pos + ref_bases.len()),
        }
    }
    /// Get the range of the reference sequence that can be skip after
    /// applying the mutation
    pub fn get_skip_range(&self) -> (usize, usize) {
        match self {
            // skip changed Nuc
            Mutation::Snp { pos, .. } => (*pos, pos + 1),
            // Insertion, skip nothing
            Mutation::Ins { pos, .. } => (*pos, *pos),
            // Skip deleted part
            Mutation::Del { pos, ref_bases } => (pos + 1, pos + ref_bases.len()),
        }
    }
    /// Create a new SNP mutation
    pub fn new_snp(pos: usize, ref_base: Nuc, alt_base: Nuc) -> Result<Self> {
        if ref_base == alt_base {
            Err(anyhow!(
                "Reference and alternate base are the same, no mutation."
            ))
        } else {
            Ok(Mutation::Snp {
                pos,
                ref_base,
                alt_base,
            })
        }
    }
    /// Create a new Insertion mutation
    pub fn new_ins(pos: usize, ref_base: Nuc, alt_bases: Vec<Nuc>) -> Result<Self> {
        if alt_bases.len() <= 1 {
            Err(anyhow!("Insertion requires at least one new base."))
        } else {
            Ok(Mutation::Ins {
                pos,
                ref_base,
                alt_bases,
            })
        }
    }
    /// Create a new Deletion mutation
    pub fn new_del(pos: usize, ref_seq: Vec<Nuc>) -> Result<Self> {
        if ref_seq.is_empty() {
            Err(anyhow!("Deletion length must be greater than 0."))
        } else {
            Ok(Mutation::Del {
                pos,
                ref_bases: ref_seq,
            })
        }
    }
    /// Get reference bases
    pub fn get_ref(self) -> String {
        match self {
            Mutation::Snp { ref_base, .. } => seq_to_string(&[ref_base]),
            Mutation::Ins { ref_base, .. } => seq_to_string(&[ref_base]),
            Mutation::Del { ref_bases, .. } => seq_to_string(&ref_bases),
        }
    }
    /// Get alternative bases
    pub fn get_alt(self) -> String {
        match self {
            Mutation::Snp { alt_base, .. } => seq_to_string(&[alt_base]),
            Mutation::Ins { alt_bases, .. } => seq_to_string(&alt_bases),
            Mutation::Del { ref_bases, .. } => seq_to_string(&ref_bases[0..1]),
        }
    }
    /// Get 1-based position
    pub fn get_1_pos(self) -> usize {
        match self {
            Mutation::Snp { pos, .. } => 1 + pos,
            Mutation::Ins { pos, .. } => 1 + pos,
            Mutation::Del { pos, .. } => 1 + pos,
        }
    }
}

pub fn sort_filter_overlap(mutations: &mut Vec<Mutation>) {
    mutations.sort_by_key(|m| m.get_interference_range().0);
    let mut last_end = 0;

    mutations.retain(|m| {
        let (start, end) = m.get_interference_range();
        if start >= last_end {
            last_end = end;
            true
        } else {
            false
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_snp() {
        let ref_seq = [Nuc::A, Nuc::C, Nuc::G, Nuc::T];
        let mutation = Mutation::new_snp(2, ref_seq[2], Nuc::A).unwrap();
        assert_eq!(mutation.get_interference_range().0, 1);
        let wrong_mutation = Mutation::new_snp(2, ref_seq[2], Nuc::G);
        assert!(wrong_mutation.is_err());
    }

    #[test]
    fn test_new_ins() {
        let ref_seq = [Nuc::A, Nuc::C, Nuc::G, Nuc::T];
        let mutation = Mutation::new_ins(2, ref_seq[2], vec![Nuc::A, Nuc::C]).unwrap();
        assert_eq!(mutation.get_interference_range().0, 2);
        let wrong_mutation = Mutation::new_ins(2, ref_seq[2], Vec::<Nuc>::new());
        assert!(wrong_mutation.is_err());
    }

    #[test]
    fn test_new_del() {
        let ref_seq = [Nuc::A, Nuc::C, Nuc::G, Nuc::T];
        let mutation = Mutation::new_del(1, ref_seq[1..3].to_vec()).unwrap();
        assert_eq!(mutation.get_interference_range().0, 1);
        let wrong_mutation = Mutation::new_del(1, ref_seq[1..1].to_vec());
        assert!(wrong_mutation.is_err());
    }
}
