use crate::utils::file_tools::open_file;
use crate::utils::types::MutByContig;
use anyhow::{Result, anyhow};
use chrono::Local;
use log::info;
use rand::Rng;
use rand::seq::index::sample;
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

/// Processes variant data and writes output files.
///
/// # Arguments
///
/// * `mutations` - The generated mutations, indexed by contig
/// * `fasta_order` - A vector of contig names in the order they
///   appear in the reference FASTA.
/// * `ploidy` - The number of copies of each chromosome present in
///   the organism.
/// * `reference_path` - The path to the reference file that this VCF
///   is showing variants from.
/// * `output_file_prefix` - The directory path and filename prefix
///   for output files.
/// * `rng` - A random number generator for this run.
pub fn write_vcf<R: Rng>(
    mutations: &MutByContig,
    fasta_order: &Vec<String>,
    ploidy: usize,
    prob_mut_multiple: f64,
    reference_path: &Path,
    overwrite_output: bool,
    output_prefix: &Path,
    rng: &mut R,
) -> Result<()> {
    let filename = output_prefix.with_extension("vcf");
    info!("Writing {}", filename.display());

    let mut outfile = open_file(&filename, overwrite_output)?;

    let vcf_headers = [
        // File format and reference
        "##fileformat=VCFv4.5",
        &format!("##fileDate={}", Local::now().format("%Y%m%d")),
        "##source=crusty-neat",
        &format!("##reference={}", reference_path.display()),
        // INFO fields
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
        // FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        // Column headers
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsim_sample",
    ];

    // Write all lines
    for line in &vcf_headers {
        writeln!(outfile, "{}", line)?;
    }

    // insert mutations
    for contig in fasta_order {
        for mutation in &mutations[contig] {
            // Mutate more than copy of the chromosome?
            let mut genotype: Vec<usize> = vec![0; ploidy];
            let mut num_ploids: usize = 1;
            let multiple_mut = rng.random_bool(prob_mut_multiple);
            if multiple_mut && ploidy > 1 {
                num_ploids = rng.random_range(2..=ploidy);
            }
            // mark mutated ploids
            for i in sample(rng, ploidy, num_ploids) {
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
        let fasta_order = vec!["chr1".to_string()];
        let ploidy = 2;
        let prob_mut_multiple = 0.1;
        let reference_path = PathBuf::from("/fake/path/to/H1N1.fa");
        let tmp_dir = TempDir::new("crusty_neat")?;
        let output_file_prefix = tmp_dir.path().join("test");
        let overwrite_output = false;
        write_vcf(
            &variants,
            &fasta_order,
            ploidy,
            prob_mut_multiple,
            &reference_path,
            overwrite_output,
            &output_file_prefix,
            &mut rng,
        )?;
        let output_file = tmp_dir.path().join("test.vcf");
        assert!(output_file.exists());
        Ok(())
    }
}
