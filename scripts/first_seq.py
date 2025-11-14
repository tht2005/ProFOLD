#!/usr/bin/env python3
import sys
from Bio import SeqIO

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <aln_file>")
        sys.exit(1)

    aln_path = sys.argv[1]
    aln = next(SeqIO.parse(aln_path, "fasta"))
    seq = "".join([c for c in str(aln.seq) if c.isupper()])
    print(seq)

if __name__ == "__main__":
    main()

