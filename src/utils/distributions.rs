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
