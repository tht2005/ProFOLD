from .msa import *
from .profold import *

import os
from Bio import SeqIO

def run_pipeline(root_dir, work_dir, log_dir, query_file, db_prefix, top_hits=None):

    output_prefix = os.path.join(work_dir, "query")

    a3m_file = run_hhblits(log_dir, query_file, db_prefix, output_prefix)
    records = list(SeqIO.parse(a3m_file, "fasta"))

    if not top_hits:
        top_hits = len(records)
    top_hits = min(top_hits, len(records) - 1)

    top_records = records[1:1+top_hits]
    print(f'Keep {len(top_records)} from {len(records)} sequences.')

    top_a3m_file = a3m_file.replace(".a3m", f"_top{top_hits}.a3m")
    SeqIO.write(top_records, top_a3m_file, "fasta")

    aln_file = f"{output_prefix}.aln"
    a3m_to_aln(log_dir, query_file, top_a3m_file, aln_file)

    profold.run_profold(root_dir, aln_file)
