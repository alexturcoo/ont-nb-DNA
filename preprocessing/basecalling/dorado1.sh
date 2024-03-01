#!/bin/bash
#SBATCH --mem=140G # MEMORY PER NODE
#SBATCH --gpus-per-node=v100:4 #specify type of GPU and number (only V100 on Beluga so it doesn't matter here)
#SBATCH --cpus-per-task=8
#SBATCH --time=70:00:00 #hrs:mins:secs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexanderturco1@gmail.com #sends me an email once done
#SBATCH --job-name=dorado1
#SBATCH --output=dorado1.o
#SBATCH --error=dorado1.e

#This script is to run oxford nanopores dorado basecaller to generate aligned basecalls from pod5 data
#This step must be done first, prior to requiggling (now known as signal mapping refinement)

module load StdEnv/2023
module load dorado
module load samtools

# Set paths as variables
BASECALL_MODEL="/home/alextu/scratch/data/basecall_models/dna_r9.4.1_e8_hac@v3.3"
REFERENCE_GENOME="/home/alextu/scratch/data/reference_genomes/ref_genome_grch38_ensemble/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"

BASECALLING_POD5S="/home/alextu/scratch/basecalling/hgsvc_pod5s/HG00268/20210903_210825_21-lee-006_PCT0053_2-A1-D1"
BAM_OUTPUT_DIR="/home/alextu/scratch/basecalling/called_aligned_bams_hg38_3/HG00268"

# Loop through the files
files=("$BASECALLING_POD5S"/*.pod5)
for file in "${files[@]}"; do

    filename=$(basename "$file" | sed 's/\.[^.]*$//')  # Extract just the file name without the path
    
    ### BASECALL + ALIGN
    dorado basecaller \
    "$BASECALL_MODEL" \
    "$file" \
    --emit-moves \
    --verbose \
    --reference "$REFERENCE_GENOME" \
    > "$BAM_OUTPUT_DIR/${filename}.bam"

    # Sort basecalled bam file
    samtools sort "$BAM_OUTPUT_DIR/${filename}.bam" > "$BAM_OUTPUT_DIR/${filename}_sorted.bam"

    # Index the sorted BAM file
    samtools index "$BAM_OUTPUT_DIR/${filename}_sorted.bam"

    # Print this to check status
    echo "1 File basecalled and aligned, bam file sorted & Indexed (file: $file)"

done