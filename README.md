<p align="center">
  <h1 align="center">🧬 ORFEX</h1>
  <p align="center">
    <strong>ORF Extraction & Annotation Pipeline</strong>
  </p>
  <p align="center">
    <em>Automated Open Reading Frame prediction and functional annotation for genomic sequences</em>
  </p>
  <p align="center">
    <img src="https://img.shields.io/badge/Language-Bash%20%7C%20Python3-blue?style=flat-square" alt="Language">
    <img src="https://img.shields.io/badge/License-MIT-green?style=flat-square" alt="License">
    <img src="https://img.shields.io/badge/Version-1.0-orange?style=flat-square" alt="Version">
  </p>
</p>

---

## 📖 Overview

**ORFEX** is an end-to-end pipeline for **predicting open reading frames (ORFs)** from genomic/metagenomic FASTA sequences and **annotating** them against a protein reference database. It chains together Prodigal for gene prediction, DIAMOND/BLASTP for homology search, and a custom Python annotation engine producing clean, organized output in a single command.

## ✨ Features

- **Automated ORF Prediction** - Uses [Prodigal](https://github.com/hyattp/Prodigal) with support for both single-genome and metagenomic modes
- **Flexible Homology Search** - Auto-detects and uses [DIAMOND](https://github.com/bbuchfink/diamond) (fast) or falls back to NCBI BLAST+ (comprehensive)
- **Rich Annotation Output** - Generates annotated protein FASTA, gene FASTA, and a detailed TSV annotation table
- **Organized Directory Structure** - All outputs neatly arranged in numbered subdirectories
- **Full Logging & Reproducibility** - Saves run parameters, software versions, and a complete pipeline log
- **Summary Report** - Auto-generates a human-readable summary with annotation statistics

## 🛠️ Prerequisites

| Tool                                                                                             | Purpose           | Install                              |
| ------------------------------------------------------------------------------------------------ | ----------------- | ------------------------------------ |
| [Prodigal](https://github.com/hyattp/Prodigal)                                                   | ORF prediction    | `conda install -c bioconda prodigal` |
| [DIAMOND](https://github.com/bbuchfink/diamond) **or** [BLAST+](https://blast.ncbi.nlm.nih.gov/) | Homology search   | `conda install -c bioconda diamond`  |
| Python 3.6+                                                                                      | Annotation script | Usually pre-installed                |

## 🚀 Quick Start

```bash
# Basic usage
./run_orf_analysis.sh -i genome.fasta -o results/ -d uniprot_sprot.fasta

# Metagenomic mode with custom e-value and threads
./run_orf_analysis.sh -i metagenome.fasta -o meta_results/ -d uniprot_sprot.fasta -m meta -e 1e-10 -t 8
```

## 📋 Usage

```
./run_orf_analysis.sh -i INPUT.fasta -o OUTPUT_DIR -d DATABASE.fasta [OPTIONS]
```

### Required Arguments

| Flag | Description                                                 |
| ---- | ----------------------------------------------------------- |
| `-i` | Input FASTA file (genome or contigs)                        |
| `-o` | Output directory (created automatically)                    |
| `-d` | Protein database in FASTA format (e.g., UniProt/Swiss-Prot) |

### Optional Arguments

| Flag | Description                           | Default |
| ---- | ------------------------------------- | ------- |
| `-m` | Prodigal mode: `single` or `meta`     | `meta`  |
| `-e` | E-value threshold for homology search | `0.001` |
| `-t` | Number of CPU threads                 | `20`    |
| `-h` | Show help message                     | —       |

## 🔬 How ORFEX Works

```
Input FASTA
    │
    ▼
┌─────────────────────────┐
│  Step 1: Prodigal        │  ORF prediction (single/meta mode)
│  Gene Calling            │  → proteins.faa, genes.fna, predictions.gff
└──────────┬──────────────┘
           │
           ▼
┌─────────────────────────┐
│  Step 2: DIAMOND/BLASTP  │  Homology search against protein database
│  Similarity Search       │  → blast_results.txt
└──────────┬──────────────┘
           │
           ▼
┌─────────────────────────┐
│  Step 3: Python Annotator│  Merge predictions + BLAST hits
│  annotate_orfs_new.py    │  → annotated FASTAs + TSV table
└──────────┬──────────────┘
           │
           ▼
┌─────────────────────────┐
│  Step 4: Summary Report  │  Statistics & organized output
└─────────────────────────┘
```

## 📂 Output Structure

```
OUTPUT_DIR/
├── input.fasta                          # Copy of input file
├── 01_prodigal/
│   ├── predictions.gff                  # ORF predictions (GFF3)
│   ├── proteins.faa                     # Predicted protein sequences
│   └── genes.fna                        # Predicted gene sequences
├── 02_blast/
│   └── blast_results.txt                # Tabular homology search results
├── 03_annotated/
│   ├── annotated_proteins.faa           # Proteins with functional annotations
│   ├── annotated_genes.fna              # Genes with functional annotations
│   └── annotation_table.tsv             # Comprehensive annotation table
├── run_info.txt                         # Pipeline parameters & software versions
├── pipeline.log                         # Complete execution log
└── SUMMARY.txt                          # Human-readable results summary
```

### Annotation Table Columns (`annotation_table.tsv`)

| Column              | Description                         |
| ------------------- | ----------------------------------- |
| `ORF_ID`            | Unique ORF identifier from Prodigal |
| `Contig`            | Source contig/scaffold name         |
| `Start` / `End`     | Genomic coordinates                 |
| `Strand`            | `+` (forward) or `-` (reverse)      |
| `Length_bp`         | Nucleotide length                   |
| `Protein_Length_aa` | Amino acid length                   |
| `Best_Hit`          | Top database match accession        |
| `Identity_%`        | Percent identity to best hit        |
| `E-value`           | Statistical significance            |
| `Bit_Score`         | Alignment quality score             |
| `Annotation`        | Functional description              |

## 📄 Pipeline Files

| File                   | Description                                                                                           |
| ---------------------- | ----------------------------------------------------------------------------------------------------- |
| `run_orf_analysis.sh`  | Main ORFEX pipeline wrapper script                                                                    |
| `annotate_orfs_new.py` | Python annotation engine - parses BLAST results, annotates FASTA headers, and generates the TSV table |

## 💡 Example

```bash
# Annotate extracted sequences against UniProt/Swiss-Prot
./run_orf_analysis.sh \
    -i Clado_extracted.fasta \
    -o Clado_ORF_results \
    -d uniprot_sprot.fasta \
    -m meta \
    -e 1e-5 \
    -t 20
```

## 📝 Notes

- DIAMOND is preferred over BLASTP for large datasets due to significantly faster runtime
- Use `-m single` for complete, single-organism genomes; use `-m meta` (default) for contigs, draft genomes, or metagenomic assemblies
- The pipeline automatically cleans up temporary database files after the BLAST step
- All intermediate files are preserved for inspection and reproducibility

---

<p align="center">
  <strong>ORFEX</strong> -- <em>From raw sequence to functional annotation, in one command.</em>
</p>
