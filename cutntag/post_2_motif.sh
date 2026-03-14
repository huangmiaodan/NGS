#!/usr/bin/env bash
set -euo pipefail

# Prompt for project path if not provided
if [ $# -lt 1 ]; then
    read -p "Enter the project directory path: " dir
else
    dir=$1
fi
cd "${dir}"

# Prompt for motif database path if not provided
if [ $# -lt 2 ]; then
    read -p "Enter the motif database path: " motif_db
else
    motif_db=$2
fi

# Prompt for target motif path if not provided
if [ $# -lt 3 ]; then
    read -p "Enter the target motif path: " tf_meme
else
    tf_meme=$3
fi

jobs=4
mkdir -p memechip fimo

# Function to process a single BAM file
process_sample() {
    local i=$1
    local sample=$(basename "${i%.bam}")

    echo "=== Processing sample: ${sample} ==="
    fasta_file="peaks/${sample}_peaks.fasta"

    meme_output="memechip/${sample}"
    if [ ! -d "${meme_output}" ]; then
        meme-chip -seed 0 -oc "${meme_output}" -time 3000 -dna -order 2 -minw 3 -maxw 8 -db ${motif_db} -meme-mod zoops -meme-nmotifs 3 -meme-searchsize 100000 -streme-pvt 0.01 -streme-totallength 4000000 -streme-align center -centrimo-score 5.0 -centrimo-ethresh 10.0 ${fasta_file}
        
        echo "✅ MEME-ChIP completed for ${sample}"
    else
        echo "  MEME-ChIP output exists for ${sample}, skipping."
    fi

    # Run FIMO motif scanning
    fimo_output="fimo/${sample}"
    if [ ! -d "${fimo_output}" ]; then
        fimo \
            --verbosity 1 \
            --oc "${fimo_output}" \
            "${tf_meme}" \
            "${fasta_file}"
        echo "✅ FIMO scanning completed for ${sample}"
    else
        echo "  FIMO output exists for ${sample}, skipping."
    fi
    echo "✅ Done ${sample}"
}

export -f process_sample
export genome t motif_db tf_meme 

# Run all BAMs in parallel
parallel -j ${jobs} process_sample ::: bam/*.bam


echo "🎉 ALL DONE"