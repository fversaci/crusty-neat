use dashmap::DashMap;

use crate::utils::{mutation::Mutation, nucleotides::Nuc};
use std::collections::HashMap;

/// Maps contig names to their nucleotide sequences
pub type SeqByContig = DashMap<String, Vec<Nuc>>;

/// Maps contig names to their mutations
pub type MutByContig = HashMap<String, Vec<Mutation>>;
