#!/usr/bin/env python3
"""
ORF Annotation Script - Combines Prodigal output with BLAST annotations
Keeps original headers and provides both nucleotide and protein annotated files
"""

import sys
import argparse

def parse_blast(blast_file):
    """Parse BLAST output - key by query ID"""
    annotations = {}
    
    with open(blast_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 13:
                continue
            
            query_id = fields[0]
            subject_id = fields[1]
            pident = fields[2]
            evalue = fields[10]
            bitscore = fields[11]
            description = fields[12] if len(fields) > 12 else "Unknown"
            
            # Take only the best hit (first occurrence)
            if query_id not in annotations:
                annotations[query_id] = {
                    'subject_id': subject_id,
                    'pident': pident,
                    'evalue': evalue,
                    'bitscore': bitscore,
                    'description': description
                }
    
    return annotations

def annotate_fasta(input_fasta, annotations, output_fasta):
    """Write FASTA with annotations in headers"""
    annotated = 0
    total = 0
    
    with open(input_fasta, 'r') as f_in, open(output_fasta, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                total += 1
                # Get the ID (everything before first space)
                header = line[1:].strip()
                seq_id = header.split()[0]
                
                # Look for annotation
                if seq_id in annotations:
                    ann = annotations[seq_id]
                    description = ann['description']
                    f_out.write(f">{seq_id} {description}\n")
                    annotated += 1
                else:
                    # Keep original header
                    f_out.write(f">{header}\n")
            else:
                f_out.write(line)
    
    return annotated, total

def parse_prodigal_header(header):
    """Extract ORF coordinates, strand, and contig from a Prodigal FASTA header.
    
    Prodigal header format:
      CONTIG|SEQ_START-SEQ_END|...|Name_N # ORF_START # ORF_END # STRAND # ;gc_cont=X
    
    Returns (orf_start, orf_end, strand, contig)
    """
    import re
    
    orf_start = 'NA'
    orf_end = 'NA'
    strand = 'NA'
    contig = 'NA'
    
    # Extract contig name (everything before first |)
    contig_match = re.match(r'^([^\|]+)', header)
    if contig_match:
        contig = contig_match.group(1)
    
    # Extract ORF start, end, strand from Prodigal's # fields
    # Format: # START # END # STRAND #
    prodigal_match = re.search(r'#\s*(\d+)\s*#\s*(\d+)\s*#\s*(-?1)\s*#', header)
    if prodigal_match:
        orf_start = prodigal_match.group(1)
        orf_end = prodigal_match.group(2)
        strand_val = int(prodigal_match.group(3))
        strand = '+' if strand_val == 1 else '-'
    
    return orf_start, orf_end, strand, contig


def create_annotation_table(proteins_fasta, genes_fasta, annotations, output_tsv):
    """Create annotation table with all details including ORF coordinates"""
    
    # Read protein sequences to get lengths and extract ORF coordinates from headers
    seq_info = {}
    current_id = None
    current_header = None
    current_seq = []
    
    with open(proteins_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    orf_start, orf_end, strand, contig = parse_prodigal_header(current_header)
                    seq_info[current_id] = {
                        'protein_length': len(''.join(current_seq)),
                        'start': orf_start,
                        'end': orf_end,
                        'strand': strand,
                        'contig': contig
                    }
                
                # Start new sequence
                current_header = line[1:].strip()
                current_id = current_header.split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        # Save last sequence
        if current_id:
            orf_start, orf_end, strand, contig = parse_prodigal_header(current_header)
            seq_info[current_id] = {
                'protein_length': len(''.join(current_seq)),
                'start': orf_start,
                'end': orf_end,
                'strand': strand,
                'contig': contig
            }
    
    # Write table
    with open(output_tsv, 'w') as out:
        out.write("ORF_ID\tContig\tStart\tEnd\tStrand\tLength_bp\tProtein_Length_aa\tBest_Hit\tIdentity_%\tE-value\tBit_Score\tAnnotation\n")
        
        for seq_id in sorted(seq_info.keys()):
            info = seq_info[seq_id]
            
            # Calculate nucleotide length
            if info['start'] != 'NA' and info['end'] != 'NA':
                length_bp = abs(int(info['end']) - int(info['start'])) + 1
            else:
                length_bp = info['protein_length'] * 3  # Approximate
            
            if seq_id in annotations:
                ann = annotations[seq_id]
                best_hit = ann['subject_id']
                pident = ann['pident']
                evalue = ann['evalue']
                bitscore = ann['bitscore']
                description = ann['description']
            else:
                best_hit = 'No hit'
                pident = 'NA'
                evalue = 'NA'
                bitscore = 'NA'
                description = 'hypothetical protein'
            
            out.write(f"{seq_id}\t{info['contig']}\t{info['start']}\t{info['end']}\t{info['strand']}\t{length_bp}\t{info['protein_length']}\t{best_hit}\t{pident}\t{evalue}\t{bitscore}\t{description}\n")

def main():
    parser = argparse.ArgumentParser(
        description='Annotate Prodigal ORFs with BLAST results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --blast results.txt --proteins orfs.faa --genes orfs.fna --output-prefix annotated
  
Output files:
  PREFIX_annotated.faa - Protein sequences with annotations
  PREFIX_annotated.fna - Gene sequences with annotations
  PREFIX_annotations.tsv - Tab-delimited table with all details
        """
    )
    parser.add_argument('--blast', required=True, help='BLAST output file (tabular format)')
    parser.add_argument('--proteins', required=True, help='Prodigal protein FASTA file')
    parser.add_argument('--genes', required=True, help='Prodigal gene nucleotide FASTA file')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("ORF ANNOTATION PIPELINE")
    print("=" * 60)
    
    print("\n[1/4] Parsing BLAST results...")
    annotations = parse_blast(args.blast)
    print(f"      ✓ Found {len(annotations)} BLAST annotations")
    
    if len(annotations) > 0:
        print("\n      Sample BLAST hits:")
        for i, (qid, ann) in enumerate(list(annotations.items())[:3]):
            print(f"        {i+1}. {qid[:60]}...")
            print(f"           → {ann['description'][:70]}...")
    
    print("\n[2/4] Annotating protein sequences...")
    prot_annotated, prot_total = annotate_fasta(
        args.proteins, 
        annotations, 
        f"{args.output_prefix}_annotated.faa"
    )
    print(f"      ✓ Annotated {prot_annotated} out of {prot_total} proteins")
    
    print("\n[3/4] Annotating nucleotide sequences...")
    gene_annotated, gene_total = annotate_fasta(
        args.genes, 
        annotations, 
        f"{args.output_prefix}_annotated.fna"
    )
    print(f"      ✓ Annotated {gene_annotated} out of {gene_total} genes")
    
    print("\n[4/4] Creating annotation table...")
    create_annotation_table(
        args.proteins,
        args.genes,
        annotations,
        f"{args.output_prefix}_annotations.tsv"
    )
    print(f"      ✓ Table created with all details")
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Total ORFs:     {prot_total}")
    print(f"  Annotated:      {prot_annotated} ({prot_annotated/prot_total*100:.1f}%)")
    print(f"  Hypothetical:   {prot_total - prot_annotated} ({(prot_total-prot_annotated)/prot_total*100:.1f}%)")
    print("\n" + "=" * 60)
    print("OUTPUT FILES")
    print("=" * 60)
    print(f"  {args.output_prefix}_annotated.faa     (protein sequences)")
    print(f"  {args.output_prefix}_annotated.fna     (nucleotide sequences)")
    print(f"  {args.output_prefix}_annotations.tsv   (annotation table)")
    print("=" * 60)
    print("\nDone!")

if __name__ == '__main__':
    main()
