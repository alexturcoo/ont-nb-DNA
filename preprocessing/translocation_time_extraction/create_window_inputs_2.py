import os
import sys
import pandas as pd
import numpy as np
import multiprocessing as mp

### This script takes all non-overlapping windows of B and non-B DNA
### and matches these windows to ONT reads in the range of the windows
### After running this you should have csv files of reads that map to
### windows of B and non-B dna created from non-b df.

# List of chromosomes to process
chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
               "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
               "chr21", "chr22", "chrX", "chrY"]

def process_chromosome(chrom):
    print("Started")
    flowcell = "20220223_220216_21-lee-006_PC24B149_3F"
    # Read non-B DNA windows
    windows_path = '/home/alextu/scratch/nonb_db_preprocessing/non_overlapping_windows_0_based'
    this_windows_path = os.path.join(windows_path, 'non_overlapping_windows_' + chrom + '.csv')
    windows_per_chr = pd.read_csv(this_windows_path)
    
    windows_per_chr = windows_per_chr.sort_values(by=['feature', 'win_start']).reset_index(drop=True)
    # Read filtered reads for the positive strand
    filtered_reads_path = '/home/alextu/scratch/compute_windows/reads_in_nonb_range/HG00733/20220223_220216_21-lee-006_PC24B149_3F'
    this_filtered_reads_path_pos = os.path.join(filtered_reads_path, f'{flowcell}{chrom}+.bed')
    filtered_reads_per_chr_pos = pd.read_csv(this_filtered_reads_path_pos, names=['chr', 'mapped_start', 'mapped_end', 'read_id', 'score', 'strand'], sep='\t')
    filtered_reads_per_chr_pos = filtered_reads_per_chr_pos[filtered_reads_per_chr_pos['score'] >= 20].reset_index(drop=True) #filter w qscore > 20
    # Read filtered reads for the negative strand
    this_filtered_reads_path_neg = os.path.join(filtered_reads_path, f'{flowcell}{chrom}-.bed')
    filtered_reads_per_chr_neg = pd.read_csv(this_filtered_reads_path_neg, names=['chr', 'mapped_start', 'mapped_end', 'read_id', 'score', 'strand'], sep='\t')
    filtered_reads_per_chr_neg = filtered_reads_per_chr_neg[filtered_reads_per_chr_neg['score'] >= 20].reset_index(drop=True) #filter w qscore > 20
    # List to store reads falling into the windows
    reads_in_windows = []
    #skipped_reads = []
    print("Start Looping")

    switch_direction = {'same': 'opposite', 'opposite': 'same'}
    #other_dir = switch_direction[direction]

    # Process reads from the positive strand
    for _, read in filtered_reads_per_chr_pos.iterrows():
        # Filter windows to intersect with the read
        relevant_windows = windows_per_chr[(windows_per_chr['win_end'] >= read['mapped_start']) & 
                                           (windows_per_chr['win_start'] <= read['mapped_end'])]

        for _, window in relevant_windows.iterrows():
            if read['mapped_end'] >= window['win_end'] and read['mapped_start'] <= window['win_start']:
                start_idx = max(0, window['win_start'] - read['mapped_start'])
                end_idx = min(window['win_end'] - read['mapped_start'] + 1, read['mapped_end'] - read['mapped_start'] + 1)
                # Compute truncated read start and end once
                truncated_read_start = read['mapped_start'] + start_idx
                truncated_read_end = read['mapped_start'] + end_idx - 1

                # Store truncated read information directly without creating intermediate dictionary
                reads_in_windows.append({
                    'chr': read['chr'],
                    'mapped_start': truncated_read_start,
                    'mapped_end': truncated_read_end,
                    'read_id': read['read_id'],
                    'score': read['score'],
                    'strand': "+",
                    'direction': window['direction'],
                    'feature': window['feature'],
                    'length': truncated_read_end - truncated_read_start
                })

    # Process reads from the negative strand
    for _, read in filtered_reads_per_chr_neg.iterrows():
        # Filter windows to intersect with the read
        relevant_windows = windows_per_chr[(windows_per_chr['win_end'] >= read['mapped_start']) & 
                                           (windows_per_chr['win_start'] <= read['mapped_end'])]

        for _, window in relevant_windows.iterrows():
            if read['mapped_end'] >= window['win_end'] and read['mapped_start'] <= window['win_start']:
                start_idx = max(0, window['win_start'] - read['mapped_start'])
                end_idx = min(window['win_end'] - read['mapped_start'] + 1, read['mapped_end'] - read['mapped_start'] + 1)
                # Compute truncated read start and end once
                truncated_read_start = read['mapped_start'] + start_idx
                truncated_read_end = read['mapped_start'] + end_idx - 1

                # Store truncated read information directly without creating intermediate dictionary
                reads_in_windows.append({
                    'chr': read['chr'],
                    'mapped_start': truncated_read_start,
                    'mapped_end': truncated_read_end,
                    'read_id': read['read_id'],
                    'score': read['score'],
                    'strand': "-",
                    'direction': switch_direction[window['direction']],
                    'feature': window['feature'],
                    'length': truncated_read_end - truncated_read_start
                })

    # Convert the list of reads in windows to a DataFrame
    reads_in_windows_df = pd.DataFrame(reads_in_windows)
    # Assuming 'mapped start' is the name of the column you want to sort by
    reads_in_windows_df_sorted = reads_in_windows_df.sort_values(by='mapped_start')
    #skipped_reads_df = pd.DataFrame(skipped_reads)
    # Export the DataFrame to a CSV file
    reads_in_windows_df_sorted.to_csv(f'/home/alextu/scratch/compute_windows/matched_windows/HG00733/{chrom}_matched_windows.csv', index=False)
    #skipped_reads_df.to_csv(f'/home/alextu/scratch/compute_windows/matched_windows_2/{chrom}_skipped_windows.csv', index=False)
    print("processed chrom", chrom)
# Use multiprocessing to process each chromosome in parallel
# Adjust the number of processes based on your system's capabilities
pool = mp.Pool(processes=24)
pool.map(process_chromosome, chromosomes)

