// Throughout this program, we are standardizing the use of a u8 representation of the nucleotides
//     A = 0
//     C = 1
//     G = 2
//     T = 3
//     N = 4
// This is intended to make it easier to store them. We thought about using the u8 representation
// of the character as built into Rust, but we'd then have to figure out the translations and keep
// track of extra numbers. So this is intended to simplify everything
use anyhow::{anyhow, Result};
use rand::distr::weighted::WeightedIndex;
use rand::distr::Distribution;
use rand::Rng;

/// This defines the relationship between the 4 possible nucleotides
/// in DNA and a simple u8 numbering system. Everything that isn't a
/// recognized base is a 4.  Note that NEAT ignores soft masking.
///
/// # Examples
///
/// ```
/// assert_eq!(base_to_u8('A'), 0);
/// assert_eq!(base_to_u8('C'), 1);
/// assert_eq!(base_to_u8('G'), 2);
/// assert_eq!(base_to_u8('T'), 3);
/// assert_eq!(base_to_u8('N'), 4);
/// ```
pub fn base_to_u8(char_of_interest: char) -> u8 {
    match char_of_interest {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'G' | 'g' => 2,
        'T' | 't' => 3,
        _ => 4,
    }
}

/// Canonical conversion from base u8 representation back into the
/// character. No attempt to preserve or display any soft masking.
///
/// # Examples
///
/// ```
/// assert_eq!(u8_to_base(0), 'A');
/// assert_eq!(u8_to_base(1), 'C');
/// assert_eq!(u8_to_base(2), 'G');
/// assert_eq!(u8_to_base(3), 'T');
/// assert_eq!(u8_to_base(4), 'N');
/// ```
pub fn u8_to_base(nuc_num: u8) -> char {
    match nuc_num {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }
}

/// Simple mutation model. The letter in the model represents
/// the base we are mutating and the vector are the weights for mutating
/// to another base (in the same a, c, g, t order)
///
/// By definition, the model is a 4x4 matrix with zeros along the diagonal
/// because, e.g., A can't "mutate" to A.
#[derive(Debug, Clone)]
pub struct NucModel {
    a: WeightedIndex<u32>,
    c: WeightedIndex<u32>,
    g: WeightedIndex<u32>,
    t: WeightedIndex<u32>,
}

impl NucModel {
    pub fn new() -> Result<Self> {
        // Default mutation model based on the original from NEAT 2.0
        Ok(Self {
            a: WeightedIndex::new(vec![0, 17, 69, 14])?,
            c: WeightedIndex::new(vec![16, 0, 17, 67])?,
            g: WeightedIndex::new(vec![67, 17, 0, 16])?,
            t: WeightedIndex::new(vec![14, 69, 16, 0])?,
        })
    }

    // todo, once we have numbers we can implement this.
    pub fn from(weights: Vec<Vec<u32>>) -> Result<Self> {
        // Supply a vector of 4 vectors that define the mutation chance
        // from the given base to the other 4 bases.
        // First some safety checks. This should be a 4x4 matrix defining mutation from
        // ACGT (top -> down) to ACGT (left -> right)
        if weights.len() != 4 {
            return Err(anyhow!("Weights supplied to NucModel is wrong size"));
        }
        for weight_vec in &weights {
            if weight_vec.len() != 4 {
                return Err(anyhow!("Weights supplied to NucModel is wrong size"));
            }
        }
        Ok(Self {
            a: WeightedIndex::new(&weights[0])?,
            c: WeightedIndex::new(&weights[1])?,
            g: WeightedIndex::new(&weights[2])?,
            t: WeightedIndex::new(&weights[3])?,
        })
    }

    pub fn choose_new_nuc<R: Rng>(&self, base: u8, rng: &mut R) -> Result<u8> {
        // Pick the weights list for the base that was input
        let dist = match base {
            0 => &self.a,
            1 => &self.c,
            2 => &self.g,
            3 => &self.t,
            // anything else we return the N value of 4
            _ => {
                return Ok(4);
            }
        };
        // Now we create a distribution from the weights and sample our choices.
        let r = dist.sample(rng).try_into()?;
        Ok(r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::create_rng;

    #[test]
    fn test_nuc_model_from_weights() -> Result<()> {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        let mut rng = create_rng(Some("Hello Cruel World"));
        let test_model = NucModel::from(vec![a_weights, c_weights, g_weights, t_weights])?;
        // It actually mutates the base
        assert_ne!(test_model.choose_new_nuc(0, &mut rng).unwrap(), 0);
        assert_ne!(test_model.choose_new_nuc(1, &mut rng).unwrap(), 1);
        assert_ne!(test_model.choose_new_nuc(2, &mut rng).unwrap(), 2);
        assert_ne!(test_model.choose_new_nuc(3, &mut rng).unwrap(), 3);
        // It gives back N when you give it N
        assert_eq!(test_model.choose_new_nuc(4, &mut rng).unwrap(), 4);
        Ok(())
    }

    #[test]
    fn test_nuc_model_too_many_vecs() -> Result<()> {
        let a_weights = vec![0, 20, 1, 20];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        let u_weights = vec![20, 1, 20, 0];
        let er = NucModel::from(vec![a_weights, c_weights, g_weights, t_weights, u_weights]);
        assert!(er.is_err());
        Ok(())
    }

    #[test]
    fn test_nuc_model_too_many_bases() -> Result<()> {
        let a_weights = vec![0, 20, 1, 20, 1];
        let c_weights = vec![20, 0, 1, 1];
        let g_weights = vec![1, 1, 0, 20];
        let t_weights = vec![20, 1, 20, 0];
        let er = NucModel::from(vec![a_weights, c_weights, g_weights, t_weights]);
        assert!(er.is_err());
        Ok(())
    }
}
