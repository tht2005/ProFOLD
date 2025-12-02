#!/bin/bash

if [ $# -le 2 ]; then
    echo "Usage $0 <MSA> <output_dir>"
    exit 1
fi

BINROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
aln=$1
outdir=$(readlink -f "$2")
_filename=$(basename -- "$aln")
target="${_filename%.*}"
fasta=$outdir/$target.fasta

n_workers=${3:-8}
n_structs=${4:-20}
n_iter=${5:-100}

mkdir -p "$outdir"

echo ">$target" > "$fasta"
# head -1 $aln >> $fasta
python3 "$BINROOT/scripts/first_seq.py" "$aln" >> "$fasta"

echo "Predict distance--------------------------------------------------------"
"$BINROOT/distance_prediction/run_inference.py" \
    -m "$BINROOT/distance_prediction/model" \
    -i "$aln" \
    -o "$outdir/$target.npz"
feat=$outdir/$target.npz
if [ ! -e "$feat" ]; then
    echo "Predict distance failed... Stop"
    exit 1
fi

echo "Generate structure by gradient descent----------------------------------"
"$BINROOT/folding/run_builder.py" \
    -i "$fasta" \
    -f "$feat" \
    -o "$outdir" \
    --n_workers $n_workers \
    --n_structs $n_structs \
    --n_iter $n_iter

echo "Full-atom relax---------------------------------------------------------"
"$BINROOT/folding/run_relax.py" \
    -s "$fasta" \
    -f "$feat" \
    -i "$outdir/final" \
    -o "$outdir/relax" \
    --n_workers $n_workers

echo "Ranking decoy-----------------------------------------------------------"
find "$outdir/relax" -name '*.pdb' | while read -r LINE; do
    echo "$LINE" "$(grep "^pose" "$LINE" | awk '{print $NF}')"
done | sort -k 2 -n > "$outdir/rank.txt"
