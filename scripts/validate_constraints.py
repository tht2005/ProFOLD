#!/usr/bin/env python3
import re
import math
import sys
from Bio import SeqIO

CST = "/dev/shm/wfs2k56e/minimize.cst"
FA = "work_dir/query.fasta"   # adjust path if needed

# load query length
try:
    rec = next(SeqIO.parse(FA, "fasta"))
    qlen = len(rec.seq)
except Exception as e:
    print("ERROR reading fasta", FA, e)
    sys.exit(2)

print("Query length:", qlen)
print("Checking", CST)

int_re = re.compile(r'\d+')
float_re = re.compile(r'[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?')

def is_number(s):
    try:
        x = float(s)
    except:
        return False
    return not (math.isnan(x) or math.isinf(x))

bad = []
counts = {"AtomPair":0, "Angle":0, "Dihedral":0, "Other":0}

with open(CST) as fh:
    for i,ln in enumerate(fh, start=1):
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        m = re.match(r'^(AtomPair|Angle|Dihedral)\b', ln)
        if not m:
            counts["Other"] += 1
            continue
        typ = m.group(1)
        counts[typ] += 1
        cols = ln.split()
        # Basic parsing rules (adapt if your format differs)
        # AtomPair example: AtomPair CA 10 CB 30 HARMONIC 3.1 0.5
        if typ == "AtomPair":
            # expect at least: AtomPair [atom1] [res1] [atom2] [res2] <fun> <params...>
            if len(cols) < 6:
                bad.append((i, typ, "too few tokens", ln))
                continue
            # tokens 2 and 4 should be integers (residue indices)
            try:
                r1 = int(cols[2])
                r2 = int(cols[4])
            except:
                # sometimes token positions are different (AtomPair CA 10 CB 30)
                # try find integers in line
                ints = int_re.findall(ln)
                if len(ints) < 2:
                    bad.append((i, typ, "no residue indices", ln))
                    continue
                r1 = int(ints[0]); r2 = int(ints[1])
            if r1 < 1 or r2 < 1 or r1 > qlen or r2 > qlen:
                bad.append((i, typ, f"res idx out of range: {r1},{r2}", ln))
            # check numeric params exist
            floats = float_re.findall(ln)
            # some AtomPairs can have many numbers; ensure at least one numeric param
            if not floats:
                bad.append((i, typ, "no numeric params", ln))

        elif typ == "Angle":
            # Angle atom1 res1 atom2 res2 atom3 res3 FUNC params...
            # expect at least 8 tokens normally; find three integers
            ints = int_re.findall(ln)
            if len(ints) < 3:
                bad.append((i, typ, "angle: fewer than 3 residue indices", ln))
                continue
            r = list(map(int, ints[:3]))
            if any(x < 1 or x > qlen for x in r):
                bad.append((i, typ, f"angle res idx out of range: {r}", ln))
            # check numeric params
            floats = float_re.findall(ln)
            if not floats:
                bad.append((i, typ, "no numeric params", ln))

        elif typ == "Dihedral":
            # Dihedral expects 4 residue indices
            ints = int_re.findall(ln)
            if len(ints) < 4:
                bad.append((i, typ, "dihedral: fewer than 4 residue indices", ln))
                continue
            r = list(map(int, ints[:4]))
            if any(x < 1 or x > qlen for x in r):
                bad.append((i, typ, f"dihedral res idx out of range: {r}", ln))
            floats = float_re.findall(ln)
            if not floats:
                bad.append((i, typ, "no numeric params", ln))

# report
print("Counts:", counts)
if not bad:
    print("No immediate syntax/index issues found in", CST)
else:
    print(f"Found {len(bad)} problematic constraint lines (showing up to 30):")
    for rec in bad[:30]:
        print(f"line {rec[0]}: {rec[1]} -> {rec[2]}")
        print("   ", rec[3][:200])
    # optionally write to file
    with open("constraints_issues.txt","w") as w:
        for rec in bad:
            w.write(f"{rec[0]}\t{rec[1]}\t{rec[2]}\t{rec[3]}\n")
    print("Wrote constraints_issues.txt")

