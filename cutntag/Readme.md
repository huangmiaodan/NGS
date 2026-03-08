# CUT&Tag Analysis Pipeline

## Preprocessing
A reproducible pipeline for processing CUT&Tag sequencing data, including:
	•	adapter trimming
	•	genome alignment
	•	BAM processing
	•	coverage track generation (bigWig)

The pipeline is optimized for HPC environments using SLURM and GNU Parallel.

⸻

### Overview

This workflow performs the following steps:
	1.	Adapter trimming using Trim Galore
	2.	Alignment to the reference genome using Bowtie2
	3.	Sorting and indexing BAM files using Samtools
	4.	Generating normalized bigWig coverage tracks using deepTools

The pipeline processes multiple samples in parallel to efficiently use HPC resources.

⸻

### Requirements

Software

Install the following tools:
	•	Trim Galore
	•	Bowtie2
	•	Samtools
	•	deepTools
	•	GNU Parallel

Optional (local installation via Conda)
```bash
conda create -n cutntag -c bioconda -c conda-forge python=3.11 \
trim-galore bowtie2 samtools deeptools parallel meme -y
```
Activate the environment:
```bash
conda activate cutntag
```

⸻

### Input Data

The pipeline expects paired-end FASTQ files named as:
```txt
sample_R1.fq.gz
sample_R2.fq.gz
```
All FASTQ files should be placed in:
```txt
raw/
```
Example:
```txt
raw/
├── sample1_R1.fq.gz
├── sample1_R2.fq.gz
├── sample2_R1.fq.gz
├── sample2_R2.fq.gz
```

⸻

### Project Structure

Recommended directory structure:
```txt
project/
│
├── raw/            # raw FASTQ files
├── trimmed/        # trimmed reads
├── bam/            # aligned BAM files
├── bw/             # bigWig coverage tracks
│
├── preprocess.sh   # alignment pipeline
└── README.md
```
Output folders will be automatically created.

⸻

### Running the Pipeline (HPC)

Submit the preprocessing job using SLURM:

sbatch preprocess.sh

The pipeline will:
	•	process 3 samples simultaneously
	•	allocate 12 threads per sample
	•	use 36 total CPUs

⸻

### Pipeline Script

preprocess.sh
```bash
#!/bin/bash
#SBATCH -J cutntag
#SBATCH -N 1
#SBATCH --cpus-per-task=36
#SBATCH -p node
#SBATCH -D path
#SBATCH -o cutntag_%j.log
#SBATCH -e cutntag_%j.err

set -eo pipefail

dir=data/dir
cd $dir

hg38=genome/dir
THREADS_PER_SAMPLE=12

mkdir -p trimmed bam bw

echo "Running on $(hostname)"
echo "Total CPUs: $SLURM_CPUS_PER_TASK"

parallel -j 3 '

sample=$(basename {= s/_R1.fq.gz// =})

echo "=== Processing ${sample} ==="

trim_galore \
-q 25 \
--cores 4 \
--path_to_cutadapt ~/.local/bin/cutadapt \
--paired raw/${sample}_R1.fq.gz raw/${sample}_R2.fq.gz \
-o trimmed/

bowtie2 \
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
```

⸻

### Output Files
```txt
File	Description
trimmed/*.fq.gz	Adapter-trimmed reads
bam/*.bam	Sorted alignment files
bam/*.bam.bai	BAM index
bw/*.bw	CPM-normalized coverage track
```

⸻

### Visualization

The generated bigWig files can be visualized in genome browsers such as:
	•	IGV
	•	UCSC Genome Browser

⸻

### Troubleshooting

Job fails immediately

Check the SLURM log files:
```txt
cutntag_JOBID.log
cutntag_JOBID.err
```
Low alignment rate

Possible causes:
	•	incorrect genome index
	•	adapter contamination
	•	poor read quality

Check trimming reports produced by Trim Galore.

⸻

### Notes
	•	Samples are automatically detected from raw/*_R1.fq.gz.
	•	Ensure sufficient CPU allocation when adjusting parallel -j.

⸻

## Post-analysis