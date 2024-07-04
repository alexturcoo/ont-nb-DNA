# ont-nb-DNA
This repository is home to a variety of scripts utilized for preprocessing and analyzing Oxford Nanopore Sequencing (ONT) Data in the context of Non-B DNA Structures.

# PREPROCESSING
This Image is the overall preprocessing workflow
![workflow](/imgs/workflow.png)

## Converting Fast5 ONT files to POD5
The FAST5 format is the standard sequencing output for Oxford Nanopore sequencers. Based on the hierarchical data format HDF5 format which enables storage of large and comples data.
Remember, these are binary so you need tools to view the contents of the data (h5py package in python)

For optimal performance, Dorado requires POD5 file input. We first convert directories of fast5s to pod5s prior to basecalling and alignment.
  * The script can be found at ```preprocessing/basecalling/```
```bash
module load python/3.10 gcc arrow/11
python -c "import pyarrow"

pod5 convert fast5 \
    /home/alextu/scratch/data/hgsvc_data/HG01457/20211022_211012_21-lee-006_PCT0053_2-A5-D5/fast5_pass/*.fast5 \
    --output /home/alextu/scratch/basecalling/hgsvc_pod5s/HG01457/20211022_211012_21-lee-006_PCT0053_2-A5-D5 \
    --one-to-one /home/alextu/scratch/data/hgsvc_data/HG01457/20211022_211012_21-lee-006_PCT0053_2-A5-D5/fast5_pass # --one-to-one command matches fast5 file name to pod5
```

## Basecalling + Alignment with Dorado
There are a lot of options for basecalling Oxford Nanopore Sequencing Data. Albacore and Guppy are commonly used however Albacore is outdated and Guppy has
been shown to be outperformed by Oxford Nanopore's Newest Open Source Basecaller called Dorado. Dorado can also call modified bases (5mC, 5hmC, etc).

[aws benchmarks, Dorado Vs Guppy](https://aws.amazon.com/blogs/hpc/benchmarking-the-oxford-nanopore-technologies-basecallers-on-aws/#:~:text=Dorado%20delivers%20significantly%20higher%20performance,instance%20type%20tested%2C%20the%20p4d.)

In this code we run Dorado with the model `dna_r9.4.1_e8_hac@v3.3.`\
`r9.4.1` = pore type (flowcell type). ONT's newest flowcell is 10.4.1)\
`e8` = Chemistry Type\
`hac` = balanced model choice (fast, hac, sup are the 3 model choices, hac being a balance between speed and accuracy)\
 * The script can be found at ```preprocessing/basecalling/```
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

## Preparing non-b windows from the non-b database (USING HG38 REFERENCE GENOME)
Run these scripts in the following order to prepare windows of B-DNA and Non-B DNA. Annotations are based off of non-B database.
These scripts will prepare 100bp windows centered around the non-B motif, and regions in between these non-B motifs are taken as regular B-DNA.
* These scripts can be found at ```preprocessing/nonb_db_preprocessing_hg38/```
1. Download .gff files from [The Non-B DNA Database](https://nonb-abcc.ncifcrf.gov/apps/ftp/browse). HG38 is the most recent annotation they have available
2. `create_nonb_dfs.py` - converts gff files into large df of non-b annotations
3. `fix_windows.py` - fixes windows around non-b motifs in the non-b motif df
4. `fix_windows_opposite.py` - fixes windows around non-b motifs on the reverse strand
5. `create_motif_free_windows.py` - find motif free fregions (B-DNA)
6. `find_nonoverlapping_windows.py` - find windows of B-DNA and Non-B DNA that do not overlap
7. `make_query_bounds.py` - Get ranges of nonB structures (smallest start to largest end) on each chromosome

## Computing Translocation Times
Run these scripts in the following order to extract translocation time metrics from ONT reads that match to the created windows.
We extract a metric called `dwell` which represents the number of sample points assigned to a base. We divide this dwell metric
by the sample rate in order to calculate translocation times (time it took for base to pass through nanopore)
* These scripts can be found at ```preprocessing/translocation_time_extraction/```

1. `filter_query_samtools.py` - Filter Basecalled + Aligned Reads in ranges of non-b structures on each chromosome
2. `create_window_inputs.py` - Match Filtered ONT Reads to Non-B Windows
3. `metric_extraction.py` - Extract translocation time metrics across all reads
4. `extract_tts.py` - Extract specific window translocation times of matched reads

## Preprocessing Plots
Run these scripts to produce plots and tables to analyze the preprocessed data.
* These scripts can be found at ```preprocessing/preprocessing_ploting```

1. `basic_read_plotting.py` - basic signal plotting with Remora (Without signal mapping refinement), simple pod5 raw signal plotting
2. `signal_mapping_refinement.py` - signal plotting with Remora (With Signal Mapping Refinement)
3. `heatmap_differences.py` - display a heatmap showing the differences between annotations/counts (on hg38 vs chm13) of nonB structures on each chromosome 
4. `make_plots_nonb_motiflength_distributions.py` - plot the motif length distribution of non-b DNA from the non-B DNA database
5. `nonb_exploration.py` - Get the # of each type of feature, find longest and shortest non-B motif, plot stacked bar of features per chromosome, individual heatmaps of structures per chrom
6. `compare_tt_quantiles.py` - Compare quantiles of exracted translocation time values for different B vs Non-B DNA structures

Example of script #6 output to show quantiles of translocation time values for specific non-B DNA motifs
![quantile_plot](/imgs/G_Quadruplex_Motif_Control_HG00268_chrY_chr22.png)

## Creating A Basic Feedforward Neural Net To Classify non-B DNA Structures based on Translocation Times
These scripts are used to create the dataand develop a simple Neural Net For Classification of non-B DNA structures
* These scripts can be found at

1. `prepare_datasets_onedirection.py` - script to produce data formatted for model (calculates median translocation time for windows w/ coverage > 5 reads for one type of strand)
2. `model_play.ipynb` - A jupyter notebook containing code to format, scale, build, deploy, and assess a simple feed forward neural network for classification

Highlighting the results of the classification model
![neuralnet_results](/imgsresults_lr_1e-05_epochs_200_batch_32_noweight.png)

