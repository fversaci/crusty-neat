use super::nucleotides::Nuc;
use anyhow::{anyhow, Result};
use serde::{Deserialize, Serialize};
use strum_macros::{Display, EnumIter, EnumString};

/// A genomic mutation (SNP, insertion, deletion)
#[derive(Debug, Clone)]
pub enum Mutation<'a> {
    Snp {
        pos: usize,
        ref_base: &'a Nuc,
        alt_base: Nuc,
    },
    Ins {
        pos: usize,
        ref_base: &'a Nuc,
        alt_bases: Vec<Nuc>, // starting with ref_base
    },
    Del {
        pos: usize,
        ref_seq: &'a [Nuc], // deleted sequence
        len: usize,
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

impl<'a> Mutation<'a> {
    /// Get initial position of the mutation
    pub fn get_pos(&self) -> usize {
        match self {
            Mutation::Snp { pos, .. } => *pos,
            Mutation::Ins { pos, .. } => *pos,
            Mutation::Del { pos, .. } => *pos,
        }
    }

    /// Get the window of the reference sequence that should not be
    /// affected by other mutations
    pub fn get_interference_window(&self) -> (usize, usize) {
        match self {
            // trinucleotide-context SNPs depend on three positions
            Mutation::Snp { pos, .. } => (pos - 1, pos + 1),
            // Keep the nucleotide before the insertion
            Mutation::Ins { pos, .. } => (*pos, *pos),
            // Keep the nucleotide before the deletion and don't
            // mutate anything else within the deleted subsequence
            Mutation::Del { pos, len, .. } => (*pos, pos + len),
        }
    }

    /// Create a new SNP mutation
    pub fn new_snp(ref_seq: &'a [Nuc], pos: usize, alt: Nuc) -> Result<Self> {
        if ref_seq[pos] == alt {
            Err(anyhow!(
                "Reference and alternate base are the same, no mutation."
            ))
        } else {
            Ok(Mutation::Snp {
                pos,
                ref_base: &ref_seq[pos],
                alt_base: alt,
            })
        }
    }

    /// Create a new Insertion mutation
    pub fn new_ins(ref_seq: &'a [Nuc], pos: usize, alt_bases: Vec<Nuc>) -> Result<Self> {
        if alt_bases.len() <= 1 {
            Err(anyhow!("Insertion requires at least one new base."))
        } else {
            Ok(Mutation::Ins {
                pos,
                ref_base: &ref_seq[pos],
                alt_bases,
            })
        }
    }

    /// Create a new Deletion mutation
    pub fn new_del(ref_seq: &'a [Nuc], pos: usize, len: usize) -> Result<Self> {
        if len == 0 {
            Err(anyhow!("Deletion length must be greater than 0."))
        } else if pos + len > ref_seq.len() {
            Err(anyhow!("Deletion exceeds reference sequence bounds."))
        } else {
            Ok(Mutation::Del {
                pos,
                ref_seq: &ref_seq[pos..pos + len],
                len,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_snp() {
        let ref_seq = vec![Nuc::A, Nuc::C, Nuc::G, Nuc::T];
        let mutation = Mutation::new_snp(&ref_seq, 2, Nuc::A).unwrap();
        assert_eq!(mutation.get_pos(), 2);
        let wrong_mutation = Mutation::new_snp(&ref_seq, 2, Nuc::G);
        assert!(wrong_mutation.is_err());
    }

    #[test]
    fn test_new_ins() {
        let ref_seq = vec![Nuc::A, Nuc::C, Nuc::G, Nuc::T];
        let mutation = Mutation::new_ins(&ref_seq, 2, vec![Nuc::A, Nuc::C]).unwrap();
        assert_eq!(mutation.get_pos(), 2);
        let wrong_mutation = Mutation::new_ins(&ref_seq, 2, Vec::<Nuc>::new());
        assert!(wrong_mutation.is_err());
    }

    #[test]
    fn test_new_del() {
        let ref_seq = vec![Nuc::A, Nuc::C, Nuc::G, Nuc::T];
        let mutation = Mutation::new_del(&ref_seq, 1, 2).unwrap();
        assert_eq!(mutation.get_pos(), 1);
        let wrong_mutation = Mutation::new_del(&ref_seq, 1, 0);
        assert!(wrong_mutation.is_err());
        let wrong_mutation = Mutation::new_del(&ref_seq, 2, 3);
        assert!(wrong_mutation.is_err());
    }
}
