import subprocess
import os
from Bio import SeqIO

def run_hhblits(log_dir, query_file, db_prefix, output_prefix, n_iter=3, e_value=1e-3, maxseq=500):
    """
    Run HHblits with given query and database prefix.
    """
    a3m_file = f"{output_prefix}.a3m"
    hhr_file = f"{output_prefix}.hhr"

    cmd = [
        "hhblits",
        "-i", query_file,
        "-d", db_prefix,
        "-oa3m", a3m_file,
        "-o", hhr_file,
        "-n", str(n_iter),
        "-e", str(e_value),
        "-maxseq", str(maxseq),
    ]

    print("Running HHblits:", " ".join(cmd))
    log_file = os.path.join(log_dir, "hhblits.log")
    with open(log_file, "a") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
        ret_code = proc.wait()
    if ret_code != 0:
        raise RuntimeError(f"hhblits failed with exit code {ret_code}. See {log_file} for details.")
    return a3m_file

def a3m_to_aln(log_dir, query_fasta, a3m_file, aln_file):
    """
    Convert an A3M file to a ProFOLD-compatible ALN using MAFFT.
    
    Parameters:
        log_dir (str): directory for logs
        query_fasta (str): query sequence FASTA (gapless)
        a3m_file (str): A3M file from HHblits
        aln_file (str): output ALN file (FASTA)
    """
    log_file = os.path.join(log_dir, "mafft.log")

    cmd = [
        "mafft",
        "--auto",
        "--addfragments", a3m_file,
        "--keeplength", query_fasta
    ]

    print("Running MAFFT:", " ".join(cmd))
    with open(aln_file, "w") as out_f, open(log_file, "a") as log:
        proc = subprocess.Popen(cmd, stdout=out_f, stderr=log, text=True)
        ret_code = proc.wait()

    if ret_code != 0:
        raise RuntimeError(f"MAFFT failed with exit code {ret_code}. See {log_file} for details.")

    print(f"Converted {a3m_file} + {query_fasta} -> {aln_file} using MAFFT")

