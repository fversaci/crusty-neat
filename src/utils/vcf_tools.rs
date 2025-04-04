use crate::utils::file_tools::open_file;
use crate::utils::types::MutByContig;
use anyhow::{Result, anyhow};
use chrono::Local;
use log::info;
use rand::Rng;
use rand::seq::index::sample;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// Converts a vector of 0s and 1s representing genotype to a standard
/// VCF genotype string.
///
/// # Arguments
///
/// * `genotype` - A vector of 0s and 1s representing genotype.
///
/// # Returns
///
/// Returns a string representing the genotype in VCF format.
fn genotype_to_string(genotype: Vec<usize>) -> Result<String> {
    let mut geno_string = String::new();
    for ploid in genotype {
        geno_string += &format!("{}/", ploid)
    }
    Ok(geno_string
        .strip_suffix("/")
        .ok_or_else(|| anyhow!("suffix not found"))?
        .to_string())
}

pub struct VcfParams<'a> {
    /// Lengths of contigs in the reference genome.
    pub contig_lengths: HashMap<String, usize>,
    /// The generated mutations, indexed by contig.
    pub mutations: &'a MutByContig,
    /// The order of contigs in the reference genome.
    pub contig_order: &'a Vec<String>,
    /// The number of copies of each chromosome present in the organism.
    pub ploidy: usize,
    /// The probability of generating the same mutation in multiple
    /// copies of the chromosome.
    pub prob_mut_multiple: f64,
    /// The path to the reference file that this VCF is showing
    /// variants from.
    pub reference_path: &'a Path,
    /// The directory path and filename prefix for output files.
    pub output_prefix: &'a Path,
    /// The flag to overwrite existing output files.
    pub overwrite_output: bool,
}

/// Processes variant data and writes output files.
pub fn write_vcf<R: Rng>(p: VcfParams, rng: &mut R) -> Result<()> {
    let filename = p.output_prefix.with_extension("vcf");
    info!("Writing {}", filename.display());
    let gzipped = false;
    let mut outfile = open_file(&filename, p.overwrite_output, gzipped)?;

    let mut vcf_headers: Vec<String> = vec![
        // File format and reference
        "##fileformat=VCFv4.3".to_string(),
        format!("##fileDate={}", Local::now().format("%Y%m%d")),
        "##source=crusty-neat".to_string(),
        format!("##reference={}", p.reference_path.display()),
    ];
    // contig fields
    for (contig, length) in &p.contig_lengths {
        let line = format!("##contig=<ID={},length={}>", contig, length);
        vcf_headers.push(line);
    }
    vcf_headers.extend_from_slice(&[
        // INFO fields
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">".to_string(),
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">".to_string(),
        // FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
        // Column headers
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsim_sample".to_string(),
    ]);

    // Write all lines
    for line in &vcf_headers {
        writeln!(outfile, "{}", line)?;
    }

    // insert mutations
    for contig in p.contig_order {
        for mutation in &p.mutations[contig] {
            // Mutate more than copy of the chromosome?
            let mut genotype: Vec<usize> = vec![0; p.ploidy];
            let mut num_ploids: usize = 1;
            let multiple_mut = rng.random_bool(p.prob_mut_multiple);
            if multiple_mut && p.ploidy > 1 {
                num_ploids = rng.random_range(2..=p.ploidy);
            }
            // mark mutated ploids
            for i in sample(rng, p.ploidy, num_ploids) {
                genotype[i] = 1;
            }
            // Format the output line. Any fields without data
            // will be a simple period. Quality is set to 37
            // for all these variants.
            let line = format!(
                "{}\t{}\t.\t{}\t{}\t37\tPASS\t.\tGT\t{}",
                contig,
                mutation.clone().get_1_pos(), // 1-based
                mutation.clone().get_ref(),
                mutation.clone().get_alt(),
                genotype_to_string(genotype)?,
            );

            writeln!(&mut outfile, "{}", line)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::create_rng;
    use crate::utils::mutation::Mutation;
    use crate::utils::nucleotides::random_seq;
    use std::collections::HashMap;
    use std::path::PathBuf;
    use tempdir::TempDir;

    #[test]
    fn test_genotype_to_string() -> Result<()> {
        let genotype = vec![0, 1, 0];
        assert_eq!(String::from("0/1/0"), genotype_to_string(genotype)?);
        Ok(())
    }

    #[test]
    fn test_write_vcf() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let seq = random_seq(&mut rng, 100);
        let variants = HashMap::from([(
            "chr1".to_string(),
            vec![
                Mutation::new_snp(3, seq[3], seq[3].complement()).unwrap(),
                Mutation::new_snp(7, seq[7], seq[7].complement()).unwrap(),
            ],
        )]);
        let contig_order = vec!["chr1".to_string()];
        let ploidy = 2;
        let prob_mut_multiple = 0.1;
        let reference_path = PathBuf::from("/fake/path/to/H1N1.fa");
        let tmp_dir = TempDir::new("crusty_neat")?;
        let output_file_prefix = tmp_dir.path().join("test");
        let overwrite_output = false;
        let vcf_params = VcfParams {
            mutations: &variants,
            contig_lengths: HashMap::from([("chr1".to_string(), 100)]),
            contig_order: &contig_order,
            ploidy,
            prob_mut_multiple,
            reference_path: &reference_path,
            overwrite_output,
            output_prefix: &output_file_prefix,
        };
        write_vcf(vcf_params, &mut rng)?;
        let output_file = tmp_dir.path().join("test.vcf");
        assert!(output_file.exists());
        Ok(())
    }
}
