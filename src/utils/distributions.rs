use anyhow::Result;
use rand::Rng;
use rand_distr::{Distribution, Normal};

/// A distribution of integers.
#[derive(Debug, Clone)]
pub enum IntDistribution {
    Normal(Normal<f64>),
    ConstantDistribution(i64),
}

impl IntDistribution {
    pub fn new_normal(mean: f64, std_dev: f64) -> Result<Self> {
        Ok(IntDistribution::Normal(Normal::new(mean, std_dev)?))
    }
    pub fn new_constant(val: i64) -> Self {
        IntDistribution::ConstantDistribution(val)
    }
}

impl Distribution<i64> for IntDistribution {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> i64 {
        match self {
            // We round the value to the nearest integer and ensure it
            // is non-negative.
            IntDistribution::Normal(n) => n.sample(rng).max(0.).round() as i64,
            // We return the constant value.
            IntDistribution::ConstantDistribution(val) => *val,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::create_rng;

    #[test]
    fn test_int_normal() {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let theor_mean = 123.;
        let dist = IntDistribution::new_normal(theor_mean, 1.).unwrap();
        let vals: Vec<i64> = (0..2000).map(|_| dist.sample(&mut rng)).collect();
        let mean = vals.iter().sum::<i64>() as f64 / vals.len() as f64;
        assert!((mean - theor_mean).abs() < 0.1);
    }

    #[test]
    fn test_int_constant() {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let dist = IntDistribution::new_constant(123);
        let val: i64 = dist.sample(&mut rng);
        assert_eq!(val, 123);
    }
}
