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


def extract_tts_in_window(chrom):
    chromosome_number = chrom[3:]
    # Read in Windows and bam info
    # Define the directory containing the numpy arrays
    windows = pd.read_csv(f"/home/alextu/scratch/compute_windows/matched_windows/HG00733/20220228_220223_21-lee-006_PC24B149_1F/{chrom}_matched_windows.csv")
    windows_pos = windows[windows['strand'] == "-"]
    #windows_pos = windows[(windows['strand'] == "+") & (windows['feature'] != "Control")]
    #windows_neg = windows[windows['strand'] == "-"] #Get Negative Windows


    ### FILTER BASED ON COVERAGE TO ONLY PULL TTS WITH COVERAGE > 5
    # Extract counts of unique ranges
    range_counts = windows_pos.groupby(['mapped_start', 'mapped_end']).size().reset_index(name='count')

    # Filter out combinations with counts less than 5
    filtered_ranges = range_counts[range_counts['count'] >= 5]

    # Create a new DataFrame by filtering out combinations with counts less than 5
    relevant_windows = windows_pos.merge(filtered_ranges, on=['mapped_start', 'mapped_end'], how='inner')


    # Initialize sets to keep track of processed combinations for positive and negative strands
    processed_combinations_pos = set()
    #processed_combinations_neg = set()

    # Initialize an empty DataFrame to store all extracted metrics
    all_extracted_metrics_df_pos = pd.DataFrame()

    directory = "/home/alextu/scratch/compute_windows/extracted_metrics_hg00733_20220228"
    # Iterate over files in the directory
    for filename in os.listdir(directory):
        # Check if the filename starts with '19_'
        if filename.startswith(f"{chromosome_number}_") and filename.endswith("_negative.npy"):
            # Extract the second and third parts of the name
            parts = filename.split('_')
            reg_start = int(parts[1])
            reg_end = int(parts[2])


            # Now you can load the numpy array using the full filename path
            full_path = os.path.join(directory, filename)
            # LOAD IN NUMPY ARRAY of sample metrics
            samples_metrics_pos = np.load(full_path, allow_pickle=True)

            # Filter windows to intersect with the read
            relevant_windows_subset = relevant_windows[(relevant_windows['mapped_start'] >= reg_start) & 
                                                       (relevant_windows['mapped_end'] <= reg_end)]

            # Initialize an empty list to store extracted metrics DataFrames for this numpy array
            extracted_metrics_dfs = []

            # Iterate over the rows of the windows DataFrame for positive strand
            
            for index, row in relevant_windows_subset.iterrows():
                mapped_start_pos = row['mapped_start'] - reg_start
                mapped_end_pos = row['mapped_end'] - reg_start
                # Check if the combination has been processed before for positive strand
                combination_pos = (mapped_start_pos, mapped_end_pos)
                if combination_pos in processed_combinations_pos:
                    continue  # Skip processing if combination has been processed before for positive strand
                else:
                    processed_combinations_pos.add(combination_pos)  # Add combination to processed set for positive strand
                    # Initialize an empty dictionary to hold the extracted metrics for positive strand
                    extracted_metrics_dict_pos = {}
                    # Iterate over the metrics within the specified range for positive strand
                    for key, value in samples_metrics_pos[0].items():
                        if isinstance(value, np.ndarray):
                            # Extract the metrics within the specified range for positive strand
                            metrics_within_range = value[:, mapped_start_pos:mapped_end_pos + 1]
                            # Create a column for each metric within the range for positive strand
                            for i in range(metrics_within_range.shape[1]):
                                column_name = f"translocation_time_{i+1}"  # Create unique column name
                                extracted_metrics_dict_pos[column_name] = metrics_within_range[:, i] / 4000 #sample rate is 4000, divide dwell/sample rate
                    # Create a DataFrame for positive strand from the dictionary
                    extracted_metrics_df_row_pos = pd.DataFrame(extracted_metrics_dict_pos)
                    # Add other necessary columns from the window information for positive strand
                    extracted_metrics_df_row_pos['chr'] = row['chr']
                    extracted_metrics_df_row_pos['strand'] = row['strand']
                    extracted_metrics_df_row_pos['feature'] = row['feature']
                    extracted_metrics_df_row_pos['mapped_start'] = mapped_start_pos + reg_start
                    extracted_metrics_df_row_pos['mapped_end'] = mapped_end_pos + reg_start
                    # Append the extracted metrics DataFrame to the list
                    extracted_metrics_dfs.append(extracted_metrics_df_row_pos)
            # Check if any extracted metrics DataFrames were created before attempting to concatenate
            if extracted_metrics_dfs:
                # Concatenate all extracted metrics DataFrames for this numpy array into a single DataFrame
                extracted_metrics_df_pos_subset = pd.concat(extracted_metrics_dfs, ignore_index=True)
                # Append the subset DataFrame to the accumulated DataFrame
                all_extracted_metrics_df_pos = pd.concat([all_extracted_metrics_df_pos, extracted_metrics_df_pos_subset], ignore_index=True)

           

    # Filter out empty rows for positive strand
    filtered_extracted_metrics_df_pos = all_extracted_metrics_df_pos.dropna()

    # Write the combined DataFrame to the output CSV file
    filtered_extracted_metrics_df_pos.to_csv(f'/home/alextu/scratch/compute_windows/translocations_dfs/HG00733/20220228_220223_21-lee-006_PC24B149_1F/tts_{chrom}_all_nonbstructures_neg.csv', index=False)

i = int(sys.argv[1])
extract_tts_in_window(chromosomes[i])