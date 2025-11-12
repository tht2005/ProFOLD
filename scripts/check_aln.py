#!/usr/bin/env python3
import sys

from Bio import SeqIO

def validate_aln(aln_path):
    seqs = [str(record.seq).upper() for record in SeqIO.parse(aln_path, "fasta")]
    AMINO = "ACDEFGHIKLMNPQRSTVWY-XBZUOJ"

    if len(seqs) < 2:
        print("Error: less than 2 sequences in the ALN file.")
        return False

    lengths = [len(seq) for seq in seqs]
    if len(set(lengths)) > 1:
        print("Error: sequences have different lengths!")
        return False

    for i, seq in enumerate(seqs, 1):
        for aa in seq:
            if aa not in AMINO:
                print(f"Invalid character '{aa}' on {i}-th sequence")
                return False

    unique_seqs = set(seqs)
    if len(unique_seqs) < 2:
        print(f"Warning: too few unique sequences: {len(unique_seqs)}")

    print("ALN file looks OK")
    return True

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <aln_file>")
        sys.exit(1)

    aln_path = sys.argv[1]
    validate_aln(aln_path)

if __name__ == "__main__":
    main()

