from operator import ne
import sys
import os
import networkx as nx
import pandas as pd
import numpy as np
from subprocess import call
from multiprocessing import Pool

### This script finds windows of b and non-B DNA that don't overlap
### Taken from Gofae-DND github, modified pieces
def get_interval(row):
    if row.motif_len > 100:
        biger_interval = (row.start, row.end)
    else:
        biger_interval = (row.win_start, row.win_end)
    return biger_interval

def find_non_overlapping_windows_pre_chr_no_graph(inp):
    chrom, strand, this_ch_stra, save_path = inp
    this_ch_stra['interval_start'] = this_ch_stra.apply(lambda row: get_interval(row)[0], axis=1)
    this_ch_stra['interval_end'] = this_ch_stra.apply(lambda row: get_interval(row)[1], axis=1)
    this_ch_stra = this_ch_stra.sort_values(by='interval_start').reset_index(drop=True)
    flag = 0
    this_ch_stra_non_overlapping = pd.DataFrame(columns=list(this_ch_stra.columns.values), index=range(len(this_ch_stra)))

    if strand == "+":
        for i in range(len(this_ch_stra) - 1):
            if this_ch_stra.loc[i, 'interval_start'] > flag and \
                    this_ch_stra.loc[i, 'interval_end'] < this_ch_stra.loc[i + 1, 'interval_start']:
                this_ch_stra_non_overlapping.loc[i, :] = this_ch_stra.loc[i, :]
            flag = max(flag, this_ch_stra.loc[i, 'interval_end'])
        this_ch_stra_non_overlapping = this_ch_stra_non_overlapping.dropna(how="all").reset_index(drop=True)
        this_ch_stra_non_overlapping.to_csv(os.path.join(save_path, 'non_overlapping_windows_' + chrom + '+.csv'), index=False)
        print('chromosome', chrom, 'strand', strand, 'done.')

    else:
        for i in range(len(this_ch_stra) - 1):
            if this_ch_stra.loc[i, 'interval_start'] > flag and \
                    this_ch_stra.loc[i, 'interval_end'] < this_ch_stra.loc[i + 1, 'interval_start']:
                this_ch_stra_non_overlapping.loc[i, :] = this_ch_stra.loc[i, :]
            flag = max(flag, this_ch_stra.loc[i, 'interval_end'])
        this_ch_stra_non_overlapping = this_ch_stra_non_overlapping.dropna(how="all").reset_index(drop=True)
        this_ch_stra_non_overlapping.to_csv(os.path.join(save_path, 'non_overlapping_windows_' + chrom + '-.csv'), index=False)
        print('chromosome', chrom, 'strand', strand, 'done.')


#### FINDING NON-OVERLAPPING WINDOWS
# path = '/labs/Aguiar/non_bdna/annotations/all_non_b_windows_0_based.csv'
path = '/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/all_non_b_windows_0_based_both_directions.csv'

save_path_pos = '/home/alextu/scratch/nonb_db_preprocessing/non_overlapping_windows_0_based/pos'
save_path_neg = '/home/alextu/scratch/nonb_db_preprocessing/non_overlapping_windows_0_based/neg'

# Read the windows dataframe
windows = pd.read_csv(path)

chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", 
               "chr6", "chr7", "chr8", "chr9", "chr10", 
               "chr11", "chr12", "chr13", "chr14", "chr15", 
               "chr16", "chr17", "chr18", "chr19", "chr20", 
               "chr21", "chr22", "chrX", "chrY"]

# Initialize inputs list
pos_inputs = []
neg_inputs = []
inputs = []

# Iterate over chromosomes and strands
for chrom in chromosomes:
    for strand in ["-"]:
        control_path = f'/home/alextu/scratch/nonb_db_preprocessing/control_windows/control_conservative_{chrom}{strand}.csv'
        controls = pd.read_csv(control_path)

        #combine non-b windows and B-dna (control) windows
        all_windows = pd.concat([windows, controls]).reset_index(drop=True)

        # Assign direction based on strand and feature
        all_windows.loc[(all_windows['feature'] == 'Control') & (all_windows['strand'] == '-'), 'direction'] = 'opposite'
        all_windows.loc[(all_windows['feature'] == 'Control'), 'motif_seq'] = ''
        all_windows.loc[(all_windows['feature'] == 'Control'), 'motif_id'] = ''

        # Convert all values in 'chr' column to strings
        all_windows['chr'] = all_windows['chr'].astype(str)
        all_windows['strand'] = all_windows['strand'].astype(str)

        this_ch_stra = all_windows[(all_windows['chr'] == chrom) & (all_windows['strand'] == strand)].reset_index(drop=True)
        neg_inputs.append([chrom, strand, this_ch_stra, save_path_neg])

    for strand in ["+"]:
        control_path = f'/home/alextu/scratch/nonb_db_preprocessing/control_windows/control_conservative_{chrom}{strand}.csv'
        controls = pd.read_csv(control_path)

        #combine non-b windows and B-dna (control) windows
        all_windows = pd.concat([windows, controls]).reset_index(drop=True)

        # Assign direction based on strand and feature
        all_windows.loc[(all_windows['feature'] == 'Control') & (all_windows['strand'] == '+'), 'direction'] = 'same'
        all_windows.loc[(all_windows['feature'] == 'Control'), 'motif_seq'] = ''
        all_windows.loc[(all_windows['feature'] == 'Control'), 'motif_id'] = ''

        # Convert all values in 'chr' column to strings
        all_windows['chr'] = all_windows['chr'].astype(str)
        all_windows['strand'] = all_windows['strand'].astype(str)

        this_ch_stra = all_windows[(all_windows['chr'] == chrom) & (all_windows['strand'] == strand)].reset_index(drop=True)
        pos_inputs.append([chrom, strand, this_ch_stra, save_path_pos])

inputs = pos_inputs + neg_inputs

# Print inputs with strand values that are negative
pool = Pool(15)
pool.map(find_non_overlapping_windows_pre_chr_no_graph, inputs)