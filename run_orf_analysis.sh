#!/bin/bash

################################################################################
# ORF Pipeline with Organized Output
# 
# This wrapper runs the ORF pipeline and saves EVERYTHING in one folder
#
# Usage: ./run_orf_analysis.sh -i INPUT.fasta -o OUTPUT_DIR -d DATABASE.fasta
################################################################################

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Help message
show_help() {
    cat << EOF
${BLUE}═══════════════════════════════════════════════════════════════════════════
ORF ANALYSIS - ORGANIZED OUTPUT VERSION
═══════════════════════════════════════════════════════════════════════════${NC}

${GREEN}USAGE:${NC}
  $0 -i INPUT.fasta -o OUTPUT_DIR -d DATABASE.fasta [OPTIONS]

${GREEN}REQUIRED ARGUMENTS:${NC}
  -i    Input FASTA file (genome/contigs)
  -o    Output directory (will be created if doesn't exist)
  -d    Protein database for annotation (FASTA format)

${GREEN}OPTIONAL ARGUMENTS:${NC}
  -m    Prodigal mode: 'single' or 'meta' (default: meta)
  -e    E-value threshold (default: 1e-5)
  -t    Number of threads (default: 4)
  -h    Show this help message

${GREEN}EXAMPLE:${NC}
  $0 -i genome.fasta -o my_analysis -d uniprot.fasta

${GREEN}OUTPUT STRUCTURE:${NC}
  OUTPUT_DIR/
  ├── input.fasta                  (copy of your input)
  ├── 01_prodigal/
  │   ├── predictions.gff          (GFF format)
  │   ├── proteins.faa             (protein sequences)
  │   └── genes.fna                (gene sequences)
  ├── 02_blast/
  │   └── blast_results.txt        (raw BLAST output)
  ├── 03_annotated/
  │   ├── annotated_proteins.faa   (proteins with annotations)
  │   ├── annotated_genes.fna      (genes with annotations)
  │   └── annotation_table.tsv     (detailed table)
  ├── run_info.txt                 (pipeline parameters)
  └── pipeline.log                 (full log)

${BLUE}═══════════════════════════════════════════════════════════════════════════${NC}
EOF
}

# Default values
PRODIGAL_MODE="meta"
EVALUE="0.001"
THREADS=20

# Parse arguments
while getopts "i:o:d:m:e:t:h" opt; do
    case $opt in
        i) INPUT_FASTA="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        d) DATABASE="$OPTARG" ;;
        m) PRODIGAL_MODE="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

# Check required arguments
if [ -z "$INPUT_FASTA" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$DATABASE" ]; then
    echo -e "${RED}Error: Missing required arguments${NC}"
    echo ""
    show_help
    exit 1
fi

# Check if input files exist
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input file '$INPUT_FASTA' not found"
    exit 1
fi

if [ ! -f "$DATABASE" ]; then
    echo "Error: Database file '$DATABASE' not found"
    exit 1
fi

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Check for required scripts
if [ ! -f "${SCRIPT_DIR}/annotate_orfs_new.py" ]; then
    echo "Error: annotate_orfs_new.py not found in ${SCRIPT_DIR}"
    exit 1
fi

# Create output directory structure
echo -e "${BLUE}Creating output directory structure...${NC}"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/01_prodigal"
mkdir -p "$OUTPUT_DIR/02_blast"
mkdir -p "$OUTPUT_DIR/03_annotated"

# Copy input file to output directory
cp "$INPUT_FASTA" "$OUTPUT_DIR/input.fasta"

# Start logging
LOG_FILE="$OUTPUT_DIR/pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo ""
echo "═══════════════════════════════════════════════════════════════════════════"
echo "ORF FINDING AND ANNOTATION PIPELINE"
echo "═══════════════════════════════════════════════════════════════════════════"
echo ""
echo "Started: $(date)"
echo ""
echo "Configuration:"
echo "  Input FASTA:     $INPUT_FASTA"
echo "  Output directory: $OUTPUT_DIR"
echo "  Database:        $DATABASE"
echo "  Prodigal mode:   $PRODIGAL_MODE"
echo "  E-value cutoff:  $EVALUE"
echo "  Threads:         $THREADS"
echo ""
echo "═══════════════════════════════════════════════════════════════════════════"
echo ""

# Save run info
cat > "$OUTPUT_DIR/run_info.txt" << EOF
ORF Analysis Run Information
============================

Date: $(date)
Pipeline version: 1.0

Input Parameters:
-----------------
Input file: $INPUT_FASTA
Database: $DATABASE
Prodigal mode: $PRODIGAL_MODE
E-value threshold: $EVALUE
Threads: $THREADS

Output Directory Structure:
---------------------------
01_prodigal/     - Prodigal ORF predictions
02_blast/        - BLAST homology search results
03_annotated/    - Final annotated sequences and table

Software Versions:
------------------
$(prodigal -v 2>&1 | head -1 || echo "Prodigal: not found")
$(diamond version 2>&1 | head -1 || echo "Diamond: not found")
$(blastp -version 2>&1 | head -1 || echo "BLAST: not found")
Python: $(python3 --version 2>&1)
EOF

# Step 1: Run Prodigal
echo "[1/4] Running Prodigal for ORF prediction..."
prodigal -i "$INPUT_FASTA" \
         -o "$OUTPUT_DIR/01_prodigal/predictions.gff" \
         -a "$OUTPUT_DIR/01_prodigal/proteins.faa" \
         -d "$OUTPUT_DIR/01_prodigal/genes.fna" \
         -f gff \
         -p "$PRODIGAL_MODE" 2>&1 | grep -v "^$" || true

NUM_ORFS=$(grep -c ">" "$OUTPUT_DIR/01_prodigal/proteins.faa" || echo "0")
echo "      ✓ Found $NUM_ORFS ORFs"
echo ""

# Step 2: Run BLAST/Diamond
echo "[2/4] Running homology search..."

if command -v diamond &> /dev/null; then
    echo "      Using Diamond..."
    
    diamond makedb --in "$DATABASE" -d "$OUTPUT_DIR/02_blast/db" --quiet
    
    diamond blastp -d "$OUTPUT_DIR/02_blast/db" \
                   -q "$OUTPUT_DIR/01_prodigal/proteins.faa" \
                   -o "$OUTPUT_DIR/02_blast/blast_results.txt" \
                   --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
                   --max-target-seqs 1 \
                   --evalue "$EVALUE" \
                   --threads "$THREADS" \
                   --quiet
    
    rm -f "$OUTPUT_DIR/02_blast/db.dmnd"
    
elif command -v blastp &> /dev/null; then
    echo "      Using BLASTP..."
    
    makeblastdb -in "$DATABASE" -dbtype prot -out "$OUTPUT_DIR/02_blast/db" > /dev/null 2>&1
    
    blastp -query "$OUTPUT_DIR/01_prodigal/proteins.faa" \
           -db "$OUTPUT_DIR/02_blast/db" \
           -out "$OUTPUT_DIR/02_blast/blast_results.txt" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
           -max_target_seqs 1 \
           -evalue "$EVALUE" \
           -num_threads "$THREADS"
    
    rm -f "$OUTPUT_DIR/02_blast/db".*
fi

NUM_HITS=$(wc -l < "$OUTPUT_DIR/02_blast/blast_results.txt" || echo "0")
echo "      ✓ Found $NUM_HITS hits"
echo ""

# Step 3: Annotate sequences
echo "[3/4] Creating annotated sequences..."

python3 "${SCRIPT_DIR}/annotate_orfs_new.py" \
    --blast "$OUTPUT_DIR/02_blast/blast_results.txt" \
    --proteins "$OUTPUT_DIR/01_prodigal/proteins.faa" \
    --genes "$OUTPUT_DIR/01_prodigal/genes.fna" \
    --output-prefix "$OUTPUT_DIR/03_annotated/annotated"

# Rename files to more descriptive names
mv "$OUTPUT_DIR/03_annotated/annotated_annotated.faa" "$OUTPUT_DIR/03_annotated/annotated_proteins.faa" 2>/dev/null || true
mv "$OUTPUT_DIR/03_annotated/annotated_annotated.fna" "$OUTPUT_DIR/03_annotated/annotated_genes.fna" 2>/dev/null || true
mv "$OUTPUT_DIR/03_annotated/annotated_annotations.tsv" "$OUTPUT_DIR/03_annotated/annotation_table.tsv" 2>/dev/null || true

echo ""
echo "[4/4] Creating summary report..."

# Count annotations
TOTAL_ORFS=$(grep -c ">" "$OUTPUT_DIR/01_prodigal/proteins.faa" || echo "0")
ANNOTATED=$(grep -v "hypothetical protein" "$OUTPUT_DIR/03_annotated/annotated_proteins.faa" | grep -c "^>" || echo "0")
HYPOTHETICAL=$((TOTAL_ORFS - ANNOTATED))

# Create summary file
cat > "$OUTPUT_DIR/SUMMARY.txt" << EOF
═══════════════════════════════════════════════════════════════════════════
ORF ANALYSIS SUMMARY
═══════════════════════════════════════════════════════════════════════════

Analysis completed: $(date)

INPUT:
------
Genome file: $INPUT_FASTA
Database: $DATABASE

RESULTS:
--------
Total ORFs predicted:        $TOTAL_ORFS
Annotated ORFs:             $ANNOTATED ($(awk "BEGIN {printf \"%.1f\", ($ANNOTATED/$TOTAL_ORFS)*100}")%)
Hypothetical proteins:       $HYPOTHETICAL ($(awk "BEGIN {printf \"%.1f\", ($HYPOTHETICAL/$TOTAL_ORFS)*100}")%)

OUTPUT FILES:
-------------
Main results (use these):
  03_annotated/annotated_proteins.faa  - Protein sequences with annotations
  03_annotated/annotated_genes.fna     - Gene sequences with annotations
  03_annotated/annotation_table.tsv    - Detailed annotation table

Intermediate files:
  01_prodigal/predictions.gff          - Prodigal GFF predictions
  01_prodigal/proteins.faa             - Raw protein sequences
  01_prodigal/genes.fna                - Raw gene sequences
  02_blast/blast_results.txt           - Raw BLAST results

Additional files:
  input.fasta                          - Copy of input file
  run_info.txt                         - Run parameters
  pipeline.log                         - Complete log
  SUMMARY.txt                          - This file

═══════════════════════════════════════════════════════════════════════════
EOF

# Display summary
cat "$OUTPUT_DIR/SUMMARY.txt"

echo ""
echo "═══════════════════════════════════════════════════════════════════════════"
echo "ALL FILES SAVED TO: $OUTPUT_DIR"
echo "═══════════════════════════════════════════════════════════════════════════"
echo ""
echo "Main results in: $OUTPUT_DIR/03_annotated/"
echo ""
