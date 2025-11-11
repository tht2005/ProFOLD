from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

def blast_search(query_fasta, db_name, out_xml="results.xml", top_hits=5, threads=2):
    """Perform BLASTP search against a local database"""
    blastp = NcbiblastpCommandline(query=query_fasta,
                                   db=db_name,
                                   evalue=1e-5,
                                   outfmt=5,
                                   out=out_xml,
                                   num_threads=threads)
    blastp()

    # Parse top hit IDs
    hit_ids = []
    with open(out_xml) as handle:
        for record in NCBIXML.parse(handle):
            for align in record.alignments[:top_hits]:
                hit_ids.append(align.hit_def.split()[0])
    return hit_ids
