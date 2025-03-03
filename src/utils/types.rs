use super::nucleotides::Nuc;
use std::collections::HashMap;

/// Maps contig names to their nucleotide sequences
pub type SeqByContig = HashMap<String, Vec<Nuc>>;
