use anyhow::{Result, anyhow};
use rand::Rng;
use rand::distr::Distribution;
use rand::distr::weighted::WeightedIndex;

/// Enum for DNA nucleotides.
#[derive(Debug, Clone, Eq, PartialEq, Copy, Hash)]
pub enum Nuc {
    A,
    C,
    G,
    T,
    N,
}

impl Nuc {
    /// Convert a nucleotide to a base character.
    pub fn to_base(self) -> char {
        nuc_to_base(self)
    }

    /// Complement function for DNA nucleotides.
    pub fn complement(self) -> Nuc {
        complement(self)
    }
}

/// Convert a base character to a nucleotide.
pub fn base_to_nuc(base: char) -> Result<Nuc> {
    match base {
        'A' | 'a' => Ok(Nuc::A),
        'C' | 'c' => Ok(Nuc::C),
        'G' | 'g' => Ok(Nuc::G),
        'T' | 't' => Ok(Nuc::T),
        'N' | 'n' => Ok(Nuc::N),
        _ => Err(anyhow!("Invalid base")),
    }
}

/// Convert a nucleotide to a base character.
pub fn nuc_to_base(nuc: Nuc) -> char {
    match nuc {
        Nuc::A => 'A',
        Nuc::C => 'C',
        Nuc::G => 'G',
        Nuc::T => 'T',
        Nuc::N => 'N',
    }
}
/// Complement function for DNA nucleotides.
pub fn complement(nuc: Nuc) -> Nuc {
    match nuc {
        Nuc::A => Nuc::T,
        Nuc::C => Nuc::G,
        Nuc::G => Nuc::C,
        Nuc::T => Nuc::A,
        Nuc::N => Nuc::N,
    }
}

/// Reverse complement function for DNA nucleotides.
pub fn reverse_complement(seq: &[Nuc]) -> Vec<Nuc> {
    seq.iter().rev().map(|&n| n.complement()).collect()
}

pub fn random_nuc<R: Rng>(rng: &mut R) -> Nuc {
    let dist = WeightedIndex::new([0.25, 0.25, 0.25, 0.25]).unwrap();
    match dist.sample(rng) {
        0 => Nuc::A,
        1 => Nuc::C,
        2 => Nuc::G,
        3 => Nuc::T,
        _ => Nuc::N,
    }
}

pub fn random_seq<R: Rng>(rng: &mut R, len: usize) -> Vec<Nuc> {
    (0..len).map(|_| random_nuc(rng)).collect()
}

pub fn seq_to_string(seq: &[Nuc]) -> String {
    seq.iter().map(|&n| nuc_to_base(n)).collect()
}

pub fn string_to_seq(s: &str) -> Result<Vec<Nuc>> {
    s.chars().map(base_to_nuc).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::create_rng;

    #[test]
    fn test_base_to_nuc() {
        assert_eq!(base_to_nuc('A').unwrap(), Nuc::A);
        assert_eq!(base_to_nuc('C').unwrap(), Nuc::C);
        assert_eq!(base_to_nuc('G').unwrap(), Nuc::G);
        assert_eq!(base_to_nuc('T').unwrap(), Nuc::T);
        assert_eq!(base_to_nuc('N').unwrap(), Nuc::N);
        assert_eq!(base_to_nuc('a').unwrap(), Nuc::A);
        assert_eq!(base_to_nuc('c').unwrap(), Nuc::C);
        assert_eq!(base_to_nuc('g').unwrap(), Nuc::G);
        assert_eq!(base_to_nuc('t').unwrap(), Nuc::T);
        assert_eq!(base_to_nuc('n').unwrap(), Nuc::N);
        assert!(base_to_nuc('X').is_err());
    }

    #[test]
    fn test_nuc_to_base() {
        assert_eq!(nuc_to_base(Nuc::A), 'A');
        assert_eq!(nuc_to_base(Nuc::C), 'C');
        assert_eq!(nuc_to_base(Nuc::G), 'G');
        assert_eq!(nuc_to_base(Nuc::T), 'T');
        assert_eq!(nuc_to_base(Nuc::N), 'N');
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement(Nuc::A), Nuc::T);
        assert_eq!(complement(Nuc::C), Nuc::G);
        assert_eq!(complement(Nuc::G), Nuc::C);
        assert_eq!(complement(Nuc::T), Nuc::A);
        assert_eq!(complement(Nuc::N), Nuc::N);
    }

    #[test]
    fn test_reverse_complement() {
        let seq = vec![Nuc::A, Nuc::C, Nuc::G, Nuc::T, Nuc::N];
        let rev_comp = reverse_complement(&seq);
        assert_eq!(rev_comp, vec![Nuc::N, Nuc::A, Nuc::C, Nuc::G, Nuc::T]);
    }

    #[test]
    fn test_random_nuc() {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let nuc = random_nuc(&mut rng);
        assert!(matches!(nuc, Nuc::A | Nuc::C | Nuc::G | Nuc::T));
    }

    #[test]
    fn test_random_seq() {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let seq = random_seq(&mut rng, 100);
        assert_eq!(seq.len(), 100);
        for &nuc in &seq {
            assert!(matches!(nuc, Nuc::A | Nuc::C | Nuc::G | Nuc::T));
        }
    }

    #[test]
    fn test_seq_to_string() {
        let seq = vec![Nuc::A, Nuc::C, Nuc::G, Nuc::T, Nuc::N];
        assert_eq!(seq_to_string(&seq), "ACGTN");
    }

    #[test]
    fn test_string_to_seq() {
        let seq = vec![Nuc::A, Nuc::C, Nuc::G, Nuc::T, Nuc::N];
        assert_eq!(string_to_seq("ACGTN").unwrap(), seq);
        assert_eq!(string_to_seq("acgtn").unwrap(), seq);
    }
}
