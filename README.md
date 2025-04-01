# Crusty-Neat Project

A diverging fork of the
[rusty_neat](https://github.com/ncsa/rusty-neat) project, designed to
simulate genomic mutations and generate sequencing reads based on a
reference genome.

In this fork, most of the original code has been rewritten to optimize
performance and shared-memory parallelism, adopting libraries such as
[DashMap](https://github.com/xacrimon/dashmap),
[Rayon](https://github.com/rayon-rs/rayon), for improved performance.

🚧 **Work in Progress**: This project is actively being developed, and
changes are ongoing.

# How to run

Clone the repository and run the software with:

```bash
cargo run -r -- -c config/crusty_test.yml --log-level debug
```

In this example:
- The `-c` flag specifies the path to the configuration file
  (`config/crusty_test.yml`).
- The `--log-level debug` flag sets the logging level to debug.

## Configuration file:

The configuration file defines the parameters for the software. The
values can be modified in the file, or you can override them directly
from the command line. Below is an example of a configuration file:

```yaml
reference: test_data/ecoli.fa
read_len: 150
coverage: 5

def_mutation_rate: 0.001
mutation_model: models/mut_model.yml
rng_seed: null

ploidy: 2
paired_ended: true
fragment_mean: 300.0
fragment_st_dev: 30.0

produce_fastq: true
produce_bam: false
produce_vcf: true
produce_fasta: true

overwrite_output: true
output_dir: /tmp/test_output
output_prefix: ecoli
```

## Overriding Configuration Values

Values defined in the configuration file can be overridden by passing
the corresponding flags directly from the command line. For example:

```bash
cargo run -r -- -c config/crusty_test.yml --log-level debug --read_len 100 --coverage 10
```

For a full list of command line parameters just run
```bash
cargo run -- -h
```

# License

This project is licensed under the **BSD-3-Clause** license.

## Authors

The original rusty_neat has been developed by the NEAT Contributors.

Crusty-Neat is being developed by:
  * Francesco Versaci, CRS4 <francesco.versaci@crs4.it>
