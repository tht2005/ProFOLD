import os
import subprocess

def prepare_blast_db(fasta_file):
    """
    Prepare a BLAST protein database from a FASTA file.
    Database files will be created in the same folder as the FASTA.
    Returns the absolute path to the database prefix.
    """
    # Decompress if needed
    if fasta_file.endswith(".gz"):
        import gzip, shutil
        decompressed = fasta_file[:-3]
        with gzip.open(fasta_file, "rt") as f_in, open(decompressed, "w") as f_out:
            shutil.copyfileobj(f_in, f_out)
        fasta_file = decompressed

    # Determine directory and base name
    fasta_dir = os.path.dirname(os.path.abspath(fasta_file))
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    db_prefix = os.path.join(fasta_dir, base_name)

    # Only make database if it doesn't exist
    if not (os.path.exists(db_prefix + ".pin") or os.path.exists(db_prefix + ".psq")):
        cmd = ["makeblastdb", "-in", fasta_file, "-dbtype", "prot", "-out", db_prefix]
        subprocess.run(cmd, check=True)

    return db_prefix  # absolute path to database prefix
