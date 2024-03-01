import sys
import os
import networkx as nx
import pandas as pd
import numpy as np
from subprocess import call

### This code is taken from Go-Fae DND repository (just rearranging based on my workflow and outputs)
### This code will handle fixing windows around the motifs from the non-b dna database (in the opposite direction minus strand)
### This code uses the csv output from fix_windows.py


### This function takes a DNA sequence as input and returns its reverse compliment
def make_reverse_complete(seq):
    complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    com_seq = ''.join([complement[base] if base in complement.keys() else base for base in seq])
    reverse_comp_seq = com_seq[::-1]
    return reverse_comp_seq

def get_opposite_direction(row):
    opposite_strand = {'+': '-', '-': '+'}
    chr, feature, start, end, strand, motif_id, motif_seq, motif_len, win_start, win_end, direction = row
    # row['chr'] = chr
    # row['feature'] = feature
    # row['start'] = start
    # row['end'] = end
    # strand = row['strand']
    row['strand'] = opposite_strand[strand]
    # row['motif_id'] = motif_id
    # motif_seq = row['motif_seq']
    row['motif_seq'] = make_reverse_complete(motif_seq)
    # row['motif_len'] = motif_len
    # row['win_start'] = win_start
    # row['win_end'] = win_end
    row['direction'] = 'opposite'
    return row

"""This code reads the windows and computes the windows on the other strand.
The sequence should be reverse complement"""
# main_folder = '/labs/Aguiar/non_bdna/annotations/' # UCHC
main_folder = '/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/' # beagle
path = os.path.join(main_folder, 'all_non_b_windows_0_based.csv')
save_path = os.path.join(main_folder, 'all_non_b_windows_0_based_both_directions.csv')
non_b_0 = pd.read_csv(path)
convert_dict = {'motif_len': int}
non_b_0 = non_b_0.astype(convert_dict)
non_b_0['direction'] = 'same'
non_b_0_opposite = non_b_0.apply(lambda row: get_opposite_direction(row), axis=1)
non_b_0_opposite.to_csv('all_non_b_windows_0_based_opposite.csv', index=False)
non_b_0_both_dir = non_b_0._append(non_b_0_opposite).reset_index(drop=True)
non_b_0_both_dir = non_b_0_both_dir.sort_values(by=['feature', 'win_start']).reset_index(drop=True)
non_b_0_both_dir.to_csv('all_non_b_windows_0_based_both_directions.csv', index=False)