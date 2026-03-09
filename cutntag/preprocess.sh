#!/bin/bash
#SBATCH -J cutntag
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -p #node#
#SBATCH -D path
#SBATCH -o cutntag_%j.log
#SBATCH -e cutntag_%j.err

set -eo pipefail

dir=#data/dir#
cd $dir
hg38=#genome/dir#
THREADS_PER_SAMPLE=24


mkdir -p trimmed bam bw

echo "Running on $(hostname)"

parallel -j 12 '
sample=$(basename {= s/_R1.fq.gz// =})

echo "=== Processing ${sample} ==="

trim_galore \
-q 25 \
--cores 4 \
--path_to_cutadapt ~/.local/bin/cutadapt \
--paired raw/${sample}_R1.fq.gz raw/${sample}_R2.fq.gz \
-o trimmed/

bowtie2 --very-sensitive --local --no-unal --no-mixed --no-discordant -I 10 -X 700 \
-p '"$THREADS_PER_SAMPLE"' \
-1 trimmed/${sample}_R1_val_1.fq.gz \
-2 trimmed/${sample}_R2_val_2.fq.gz \
-x '"$hg38"' | \
samtools sort \
-@ '"$THREADS_PER_SAMPLE"' \
-o bam/${sample}.bam

samtools index bam/${sample}.bam


bamCoverage \
-b bam/${sample}.bam \
-o bw/${sample}.bw \
--normalizeUsing CPM \
-p '"$THREADS_PER_SAMPLE"'

echo "=== Done ${sample} ==="

' ::: raw/*_R1.fq.gz

echo "ALL DONE!"