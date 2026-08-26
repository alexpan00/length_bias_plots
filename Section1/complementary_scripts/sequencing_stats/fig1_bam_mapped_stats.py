#!/usr/bin/env python3
import sys
import os
import re
import subprocess

cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')

def parse_cigar(cigar_str):
    """
    Parses a CIGAR string to calculate:
    1. aligned_read_len: number of read bases that align to reference (excludes soft-clipping)
    2. full_read_len: total sequence length represented in BAM (includes soft-clipping)
    3. aligned_ref_span: exonic span on the reference genome (excludes introns/N)
    """
    operations = cigar_pattern.findall(cigar_str)
    aligned_read_len = 0
    full_read_len = 0
    aligned_ref_span = 0
    
    for length_str, op in operations:
        length = int(length_str)
        
        # Consumes query (read) bases
        if op in ('M', 'I', 'S', '=', 'X'):
            full_read_len += length
            if op != 'S':
                aligned_read_len += length
                
        # Consumes reference bases
        if op in ('M', 'D', 'N', '=', 'X'):
            if op != 'N': # Exclude introns/reference gaps
                aligned_ref_span += length
                
    return aligned_read_len, full_read_len, aligned_ref_span

def calculate_n50(lengths):
    """Calculates N50 from a list of lengths."""
    if not lengths:
        return 0
    lengths = sorted(lengths, reverse=True)
    total_bp = sum(lengths)
    cumsum = 0
    n50_threshold = total_bp / 2
    for l in lengths:
        cumsum += l
        if cumsum >= n50_threshold:
            return l
    return 0

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 bam_mapped_stats.py <bam_file>", file=sys.stderr)
        sys.exit(1)
        
    bam_path = sys.argv[1]
    if not os.path.exists(bam_path):
        print(f"Error: file {bam_path} does not exist", file=sys.stderr)
        sys.exit(1)
        
    # Stream mapped reads only (-F 4 filters out unmapped reads)
    cmd = ["samtools", "view", "-F", "4", bam_path]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=100000)
    
    mapped_read_lengths = []
    full_read_lengths = []
    ref_spans = []
    
    total_mapped_reads = 0
    total_mapped_bases = 0
    total_full_bases = 0
    
    for line in process.stdout:
        fields = line.strip().split('\t')
        if len(fields) < 11:
            continue
            
        cigar_str = fields[5]
        if cigar_str == '*':
            continue
            
        aligned_read_len, full_read_len, aligned_ref_span = parse_cigar(cigar_str)
        
        mapped_read_lengths.append(aligned_read_len)
        full_read_lengths.append(full_read_len)
        ref_spans.append(aligned_ref_span)
        
        total_mapped_reads += 1
        total_mapped_bases += aligned_read_len
        total_full_bases += full_read_len
        
    process.stdout.close()
    process.wait()
    
    if total_mapped_reads == 0:
        print(f"Sample: {os.path.basename(bam_path)}")
        print("No mapped reads found.")
        return
        
    # Calculate statistics
    avg_mapped_len = total_mapped_bases / total_mapped_reads
    avg_full_len = total_full_bases / total_mapped_reads
    
    n50_mapped = calculate_n50(mapped_read_lengths)
    n50_full = calculate_n50(full_read_lengths)
    
    pct_bases_mapped = (total_mapped_bases / total_full_bases * 100) if total_full_bases > 0 else 0
    
    print(f"Sample\t{os.path.basename(bam_path)}")
    print(f"Mapped_Reads\t{total_mapped_reads}")
    print(f"Average_Full_Read_Length_Mapped_Only\t{avg_full_len:.2f}")
    print(f"Average_Mapped_Read_Length\t{avg_mapped_len:.2f}")
    print(f"N50_Full_Read_Length\t{n50_full}")
    print(f"N50_Mapped_Read_Length\t{n50_mapped}")
    print(f"Percent_Bases_Mapped\t{pct_bases_mapped:.2f}")

if __name__ == "__main__":
    main()
