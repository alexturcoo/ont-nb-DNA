import sys
import os
import networkx as nx
import pandas as pd
import numpy as np
from subprocess import call

### This code is taken from Go-Fae DND repository (just rearranging based on my workflow and outputs)
### This code will handle fixing windows around the motifs from the non-b dna database
### This code uses the csv output from create_nonb_dfs.py

def fix_windows_on_motifs(ms, me):
    """ms: motif start, me: motif end"""
    ml = me - ms + 1
    # motif length is even
    if ml % 2 == 0:
        a = ml/2
        b = 50 - a
        win_s = ms-b
        win_e = me+b
    # ml%2 == 1: # motif length is odd
    else:
        a = np.floor(ml/2)
        b_start = 50 - (a + 1)
        b_end = 50 - a
        win_s = ms - b_start
        win_e = me + b_end
    return [win_s, win_e]

add_non_b_ann = pd.read_csv('/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/all_non_b_motifs_0_based.csv')
add_non_b_ann['win_start'] = add_non_b_ann.apply(lambda row: fix_windows_on_motifs(row.start, row.end)[0], axis=1)
add_non_b_ann['win_end'] = add_non_b_ann.apply(lambda row: fix_windows_on_motifs(row.start, row.end)[1], axis=1)
convert_dict = {'win_start': int, 'win_end': int}
add_non_b_ann = add_non_b_ann.astype(convert_dict)
add_non_b_ann.to_csv('/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/all_non_b_windows_0_based.csv', index=False)