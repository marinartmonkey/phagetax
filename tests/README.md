# PhageTax Test Data

This folder contains test datasets and expected results to verify the correct functioning of the PhageTax pipeline.

## Structure

tests/
├── entrada_genomas/ # Real bacteriophage genome sequences for testing
├── expected_outputs/ # Optional: Example outputs for result comparison
└── README.md # This file


## Description

- **entrada_genomas/**: Contains real bacteriophage genome sequences in FASTA format. These genomes are provided for testing and validation purposes only.
- **expected_outputs/**: Optional folder with expected output files for manual comparison, such as:
  - Binary PHROG matrix
  - Genome summary
  - Final combined matrix
  - Taxonomic predictions

## Provenance of Example Genomes

The genomes included in this folder are real bacteriophage sequences obtained from public databases (e.g., NCBI RefSeq, GenBank) or scientific publications. Their use in this repository is strictly limited to testing and educational purposes.

If you plan to distribute this repository publicly, please ensure that all genomes included comply with the corresponding licenses or citation requirements.

## How to use

You can copy the test genomes into the main working directory and run PhageTax normally:

```bash
rm -rf entrada_genomas pharokka_out
cp -r tests/entrada_genomas entrada_genomas
mkdir pharokka_out

python phagetax.py annotate_batch
python phagetax.py generate_matrix
python phagetax.py genome_summary
python phagetax.py final_matrix
python phagetax.py predict
python phagetax.py phrog_report
