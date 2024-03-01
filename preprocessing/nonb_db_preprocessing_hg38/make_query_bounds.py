import sys
import os
import networkx as nx
import pandas as pd
import numpy as np
from subprocess import call

### This code is taken from Go-Fae DND repository (just rearranging based on my workflow and outputs)
### This code finds the start of the first widnwo and end of the last window on ALL chromosomes - makes bounds to filter reads in these regions only
### This code uses the csv output from fix_windows_opposite.py

# add_non_b_ann = pd.read_csv('/labs/Aguiar/non_bdna/annotations/all_non_b_windows_0_based_both_directions.csv') # UCHC
add_non_b_ann = pd.read_csv('/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/all_non_b_windows_0_based_both_directions.csv') # beagle
# mapping_file_path = '/labs/Aguiar/non_bdna/annotations/mapping_chr_chromosome.csv' # UCHC
mapping_file_path = '/home/alextu/scratch/nonb_db_preprocessing/chromosome_mappings/mapping_chr_chromosome.csv' # beagle
mapping_file = pd.read_csv(mapping_file_path, index_col=False)
chr_dict = {mapping_file.loc[i, 'chr']: mapping_file.loc[i, 'chromosome'] for i in range(len(mapping_file))}
non_b_g_min = add_non_b_ann.groupby(['chr', 'strand'])['win_start'].min()
non_b_g_max = add_non_b_ann.groupby(['chr', 'strand'])['win_end'].max()
groups = list(non_b_g_min.index)
queries_df = pd.DataFrame(columns=['chr', 'strand', 'start', 'end'], index=range(len(groups)))

i = 0
for gr in groups:
    chrom = gr[0]
    strand = gr[1]
    start = non_b_g_min.loc[gr]
    end = non_b_g_max.loc[gr]
    queries_df.loc[i, :] = chrom, strand, start, end
    i += 1
queries_df['chromosome'] = queries_df['chr'].apply(lambda x: chr_dict[x])
queries_df['start'] = queries_df['start'] + 1
queries_df['end'] = queries_df['end'] + 1
# queries_df.to_csv('/labs/Aguiar/non_bdna/annotations/chr_start_end_queries_1_based.csv', index=False) # UCHC
queries_df.to_csv('/home/alextu/scratch/nonb_db_preprocessing/chromosome_start_end_queries/chr_start_end_queries_1_based.csv', index=False) # beagle