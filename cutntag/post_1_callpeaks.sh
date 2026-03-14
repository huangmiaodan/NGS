#!/usr/bin/env bash
set -euo pipefail

# Prompt for project path if not provided
if [ $# -lt 1 ]; then
    read -p "Enter the project directory path: " dir
else
    dir=$1
fi
cd "${dir}"

# Prompt for genome path if not provided
if [ $# -lt 2 ]; then
    read -p "Enter the genome directory path: " genome
else
    genome=$2
fi


t=6  # threads per job
jobs=4  # number of parallel samples


mkdir -p qc logs marked peaks bw results results/frip

# Function to process a single BAM file
process_sample() {
    local i=$1
    local sample=$(basename "${i%.bam}")

    echo "=== Processing sample: ${sample} ==="
    # QC
    qc_dir="qc/${sample}"
    if [ ! -d "${qc_dir}" ]; then
        echo "  Running Qualimap QC for ${sample}"
        qualimap bamqc \
            -bam "bam/${sample}.bam" \
            -outdir "${qc_dir}" \
            -nt "${t}" \
            --java-mem-size=8G \
            > "logs/${sample}_qualimap.log" 2>&1
    else
        echo "  QC already exists for ${sample}, skipping Qualimap."
    fi

    # Coverage track
    bw_file="bw/${sample}.bw"
    if [ ! -f "${bw_file}" ]; then
        echo "  Generating coverage track for ${sample}"
        bamCoverage \
            -b "bam/${sample}.bam" \
            -o "${bw_file}" \
            --normalizeUsing CPM \
            --binSize 50 \
            -p "${t}" \
            > "logs/${sample}_bw.log" 2>&1
    else
    
    echo "  Coverage track already exists for ${sample}, skipping bamCoverage."
    fi


    mark_file="marked/${sample}_metrics.txt"
    if [ ! -f "${peak_file}" ]; then
        echo "Marking duplicates for ${sample}"
        picard MarkDuplicates \
            I=$i \
            O=/dev/null \
            M=marked/${sample}_metrics.txt \
            CREATE_INDEX=true \
            REMOVE_DUPLICATES=false \
            VALIDATION_STRINGENCY=SILENT
    fi

    # Peak calling: using -f BAMPE is not better for motif calling than -f BAM
    peak_file="peaks/${sample}_peaks.narrowPeak"
    if [ ! -f "${peak_file}" ]; then
        echo "  Calling peaks for ${sample}"
        macs3 callpeak \
            -t $i \
            -f BAM \
            -g hs \
            -n "peaks/${sample}" \
            -q 0.01 \
            -B --SPMR \
            > "logs/${sample}_macs3.log" 2>&1

        bedGraphToBigWig "peaks/${sample}_treat_pileup.bdg" \
            "${genome}/hg38.chrom.sizes" \
            "peaks/${sample}_SPMR_FE.bw"

        sort -k8,8nr peaks/${sample}_peaks.narrowPeak > peaks/${sample}.sorted.narrowPeak
    else
        echo "  Peaks already exist for ${sample}, skipping peak calling."
    fi

    # Convert narrowPeak to BED
    bed_file="peaks/${sample}_peaks.bed"
    if [ ! -f "${bed_file}" ]; then
        cut -f 1-3 peaks/${sample}.sorted.narrowPeak > "${bed_file}"
        echo "✅ Created BED file for ${sample}"
    else
        echo "  BED file already exists for ${sample}, skipping."
    fi


    # Calculating FRiP
    saf_file="peaks/${sample}.saf"
    if [ ! -f "${saf_file}" ]; then 
        awk 'OFS="\t" {print $1"-"$2+1"-"$3, $1, $2+1, $3, "+"}' $peak_file > peaks/${sample}.saf
        # Run featureCounts
        featureCounts -p -a peaks/${sample}.saf -F SAF -o results/frip/${sample}.txt $i
    fi

    # Extract sequences from genome
    fasta_file="peaks/${sample}_peaks.fasta"
    if [ ! -f "${fasta_file}" ]; then
        bedtools getfasta \
            -fi "${genome}/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
            -bed "${bed_file}" \
            -fo "${fasta_file}"
        fasta-unique-names -r "${fasta_file}"
        echo "✅ Extracted FASTA sequences for ${sample}"
    else
        echo "  FASTA file already exists for ${sample}, skipping."
    fi

}

export -f process_sample
export genome t

# Run all BAMs in parallel
parallel -j ${jobs} process_sample ::: bam/*.bam

frip_file="results/combined_index_FRiP.txt"
echo -e "Sample\tTotalReads\tMappedReads\tReadsInPeaks\tFRiP" > $frip_file

for file in results/frip/*.txt.summary; do
    sample=$(basename "$file" .txt.summary)

    assigned=$(awk '/Assigned/ {print $2}' "$file")
    unmapped=$(awk '/Unassigned_Unmapped/ {print $2}' "$file")
    nofeature=$(awk '/Unassigned_NoFeatures/ {print $2}' "$file")

    total=$((assigned + nofeature + unmapped))
    mapped=$((assigned + nofeature))
    frip=$(awk -v a="$assigned" -v m="$mapped" 'BEGIN {printf "%.4f", a/m}')

    echo -e "${sample}\t${total}\t${mapped}\t${assigned}\t${frip}" >> $frip_file
done

echo "🎉 ALL DONE"