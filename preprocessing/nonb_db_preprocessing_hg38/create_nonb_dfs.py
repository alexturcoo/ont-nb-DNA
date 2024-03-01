import sys
import os
import networkx as nx
import pandas as pd
import numpy as np
from subprocess import call

### This code  is from ONT-GoFAE DND preprocessing.py script
### I downloaded gff files from the non-b DNA database
### This code takes all those files and combines the information from the non-B db into a csv
### After this .csv has been produced, we can move on to fixing the windows around these motifs

"""Note gff file is indexed starting from 1 (the first position on genome is one)
    Input: nonb dna database gff file (1-based)
    Output: Dataframe of necessary information (0-based)"""
main_folder = '/home/alextu/scratch/nonb_db_preprocessing/annotations'
gff_files = [name for name in os.listdir(main_folder) if '.gff' in name and 'chrMT' not in name]
add_non_b_ann = pd.DataFrame(columns=['chr', 'feature', 'start', 'end', 'strand'])
for gf in gff_files:
    gff_df = pd.read_csv(os.path.join(main_folder, gf), names=['chr', 'source', 'feature', 'start', 'end', 'score',
                                                               'strand', 'frame', 'attribute'], skiprows=[0],
                        sep='\t')
    gff_df['motif_id'] = gff_df['attribute'].apply(lambda x: x.split(';')[0].split('=')[1])
    seq_idx = [idx for idx, s in enumerate(gff_df.loc[0, 'attribute'].split(';')) if 'sequence' in s][0]
    gff_df['motif_seq'] = gff_df['attribute'].apply(lambda x: x.split(';')[seq_idx].split('=')[1])
    gff_df['motif_len'] = gff_df.apply(lambda row: row.end - row.start + 1, axis=1)
    gff_df['start'] = gff_df['start'] - 1
    gff_df['end'] = gff_df['end'] - 1
    gff_df = gff_df.drop(columns=['source', 'score', 'frame', 'attribute'])
    gff_df = gff_df.dropna()
    add_non_b_ann = add_non_b_ann._append(gff_df, ignore_index=True).reset_index(drop=True)
add_non_b_ann.to_csv('/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/all_non_b_motifs_0_based.csv', index=False)