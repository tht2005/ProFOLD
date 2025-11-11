from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO

def run_mafft(input_fasta, output_fasta="aligned.fasta"):
    """Run MAFFT alignment"""
    mafft_cline = MafftCommandline(input=input_fasta)
    stdout, _ = mafft_cline()
    with open(output_fasta, "w") as f:
        f.write(stdout)
    alignment = AlignIO.read(output_fasta, "fasta")
    return alignment

def write_aln_file(fasta_alignment_file, aln_file="aligned.aln"):
    """
    Convert a FASTA MSA to plain .aln format suitable for ProFOLD.
    
    Each line contains the aligned sequence (no headers), gaps as '-' are preserved.
    """
    aln_content = AlignIO.read(fasta_alignment_file, "fasta")
    
    with open(aln_file, "w") as f:
        for record in aln_content:
            f.write(str(record.seq) + "\n")
    
    return aln_file, aln_content
