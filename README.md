# ont-nb-DNA
This repository is home to a variety of scripts utilized for preprocessing and analyzing Oxford Nanopore Sequencing (ONT) Data in the context of Non-B DNA Structures.

## PREPROCESSING
This Image is the overall workflow flow preprocessing... Will be here when I complete.

### Converting Fast5 ONT files to POD5
The FAST5 format is the standard sequencing output for Oxford Nanopore sequencers. Based on the hierarchical data format HDF5 format which enables storage of large and comples data.
Remember, these are binary so you need tools to view the contents of the data (h5py package in python)

For optimal performance, Dorado requires POD5 file input. We first convert directories of fast5s to pod5s prior to basecalling and alignment.\
  * The script can be found at ... ```PUT SCRIPT LOCATION IN REPO HERE```
```bash
module load python/3.10 gcc arrow/11
python -c "import pyarrow"

pod5 convert fast5 \
    /home/alextu/scratch/data/hgsvc_data/HG01457/20211022_211012_21-lee-006_PCT0053_2-A5-D5/fast5_pass/*.fast5 \
    --output /home/alextu/scratch/basecalling/hgsvc_pod5s/HG01457/20211022_211012_21-lee-006_PCT0053_2-A5-D5 \
    --one-to-one /home/alextu/scratch/data/hgsvc_data/HG01457/20211022_211012_21-lee-006_PCT0053_2-A5-D5/fast5_pass # --one-to-one command matches fast5 file name to pod5
```

### Basecalling + Alignment with Dorado
There are a lot of options for basecalling Oxford Nanopore Sequencing Data. Albacore and Guppy are commonly used however Albacore is outdated and Guppy has
been shown to be outperformed by Oxford Nanopore's Newest Open Source Basecaller called Dorado. Dorado can also call modified bases (5mC, 5hmC, etc).

[aws benchmarks, Dorado Vs Guppy](https://aws.amazon.com/blogs/hpc/benchmarking-the-oxford-nanopore-technologies-basecallers-on-aws/#:~:text=Dorado%20delivers%20significantly%20higher%20performance,instance%20type%20tested%2C%20the%20p4d.)

In this code we run Dorado with the model dna_r9.4.1_e8_hac@v3.3.\
r9.4.1 = pore type (flowcell type). ONT's newest flowcell is 10.4.1)\
e8 = Chemistry Type\
hac = balanced model choice (fast, hac, sup are the 3 model choices, hac being a balance between speed and accuracy)\
 * The script can be found at ... ```PUT SCRIPT LOCATION IN REPO HERE```
```bash
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
```
After running, you will have a corresponding bam file for each pod5 file.
You can change this to produce a single bam file for all pod5s directly however I found it
to be faster and easier to keep track of when doing one file at a time. For
Downstream processing you must `merge` the bams, and then `sort` and `index` the merged bam using samtools).

### Preparing non-b windows from the non-b database (USING HG38 REFERENCE GENOME)

1. Download .gff files from [The Non-B DNA Database](https://nonb-abcc.ncifcrf.gov/apps/ftp/browse). HG38 is the most recent annotation they have available
2. `create_nonb_dfs.py` - converts gff files into large df of non-b annotations
3. `fix_windows.py` - fixes windows around non-b motifs in the non-b motif df
4. `fix_windows_opposite.py` - fixes windows around non-b motifs on the reverse strand
5. `create_motif_free_windows.py` - find motif free fregions (B-DNA)
6. `find_nonoverlapping_windows.py` - find windows of B-DNA and Non-B DNA that do not overlap

### Preprocessing plots

## MATCHING ONT READS TO NON-B WINDOWS

### Filter Basecalled + Aligned Reads in range of non-b structures

### Matching ONT Reads to Non-B windows

### Extracting Translocation times of mapped reads

