import logging
from pathlib import Path
import pod5 as p5
import pandas as pd
import polars as pl
import plotnine as p9
import remora
from remora import io, refine_signal_map, util
import matplotlib.pyplot as plt
import pysam
import multiprocessing as mp
import sys
import numpy as np
import os

# List of chromosomes to process
chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
               "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
               "chr21", "chr22", "chrX", "chrY"]

def extract_metrics_per_chr(chromosome):

    test_data_root = Path("/home/alextu/scratch/basecalling/hgsvc_pod5s/HG00733/20220223_220216_21-lee-006_PC24B149_3F")
    pod5_dr = p5.DatasetReader(test_data_root)
    bam_fh = pysam.AlignmentFile("/home/alextu/scratch/basecalling/merged_bams/HG00733/20220223_220216_21-lee-006_PC24B149_3F_merged.bam") # Path to your merged BAM file
    chromosome_number = chromosome[3:]

    # Load in kmer Table (from remora repository, not sure if I am using the right one now)
    level_table = "/home/alextu/scratch/compute_windows/inputs/9mer_levels_v1.txt"
    sig_map_refiner = refine_signal_map.SigMapRefiner(
    kmer_model_filename=level_table,
    do_rough_rescale=True,
    scale_iters=0,
    do_fix_guage=True,)

    # Read in Windows and bam info
    windows = pd.read_csv(f"/home/alextu/scratch/compute_windows/matched_windows/HG00733/20220223_220216_21-lee-006_PC24B149_3F/{chromosome}_matched_windows.csv")
    windows_pos = windows[windows['strand'] == "+"] #Get positive windows
    #windows_neg = windows[windows['strand'] == "-"] #Get Negative Windows

    # Find the smallest mapped_start and largest mapped_end
    smallest_start_pos = windows_pos['mapped_start'].min()
    largest_end_pos = windows_pos['mapped_end'].max()

    #smallest_start_neg = windows_neg['mapped_start'].min()
    #largest_end_neg = windows_neg['mapped_end'].max()

    # Define the step size for the loop
    step_size = 10000000  # Adjust the step size as needed

    for start_pos in range(smallest_start_pos, largest_end_pos, step_size):
        end_pos = min(start_pos + step_size - 1, largest_end_pos)  # Adjust end position
        print(start_pos, end_pos)
        # Define the reference region for the chunk
        ref_reg_pos = io.RefRegion(ctg=chromosome_number, strand="+", start=start_pos, end=end_pos)

        # Define the filename for the numpy array
        npy_filename = f"/home/alextu/scratch/compute_windows/extracted_metrics_hg00733_20220223/{ref_reg_pos.ctg}_{ref_reg_pos.start}_{ref_reg_pos.end}_positive.npy"

        # Check if the file already exists
        if not os.path.exists(npy_filename):
            try:
                # Get samples metrics for the chunk
                samples_metrics_pos, all_bam_reads_pos = io.get_ref_reg_samples_metrics(
                ref_reg_pos,
                [(pod5_dr, bam_fh)],
                metric="dwell",
                sig_map_refiner=sig_map_refiner,
                missing_ok=True,  # Allows to continue even if Bam record not found in POD5

            )
            except Exception as e:
                print("An error occurred:", e)
                # If you want to continue execution despite the error, you can add a pass statement
                samples_metrics_pos = False

            if samples_metrics_pos:
                print("metric + completed for positive strand", start_pos, end_pos)

                # Save the numpy array for positive strand with chr, mapped_start, and mapped_end in the filename
                np.save(npy_filename, samples_metrics_pos)

            else:
                print("no reads covering region")
        else:
            print("File already exists:", npy_filename)

i = int(sys.argv[1])
extract_metrics_per_chr(chromosomes[i])