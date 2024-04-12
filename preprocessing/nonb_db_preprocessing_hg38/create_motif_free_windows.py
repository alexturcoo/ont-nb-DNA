import sys
import os
import pandas as pd
import numpy as np
from subprocess import call
from multiprocessing import Pool

### This script creates motif free intervals, or windows of B-DNA

def make_control_windows_chr_strand_conservative(inp):
    # conservative
    chr, stra, this_ch_stra, save_path = inp

    motif_intervals_list = list(zip(this_ch_stra.start, this_ch_stra.end))
    windows_intervals_list = list(zip(this_ch_stra.win_start, this_ch_stra.win_end))
    interval_list = motif_intervals_list + windows_intervals_list

    # find unions of regions
    union = []
    for begin, end in sorted(interval_list):
        if union and union[-1][1] >= begin - 1:
            union[-1][1] = max(union[-1][1], end)
        else:
            union.append([begin, end])

    union = sorted(union, key=lambda x: x[1])
    # find between regions
    motif_free_intervals = []
    for i in range(len(union) - 1):
        free = [union[i][1] + 1, union[i + 1][0] - 1]
        if free[1] - free[0] > 140:
            # print(free[1] - free[0])
            motif_free_intervals.append(free)

    # make the windows (bed file)
    control_df = pd.DataFrame(columns=['chr', 'feature', 'start', 'end', 'strand', 'motif_len', 'win_start', 'win_end'],
                              index=range(len(motif_free_intervals)))
    # tt = []
    for i in range(len(motif_free_intervals)):
        start = motif_free_intervals[i][0]
        end = motif_free_intervals[i][1]
        mid = int(np.floor((start + end) / 2))
        win_start = mid - 49
        win_end = mid + 50
        # tt.append((win_start, win_end))
        # print(motif_free_intervals[i],'==>',(win_start, win_end))
        control_df.loc[i, :] = chr, 'Control', win_start, win_end, stra, 100, win_start, win_end

    control_df.to_csv(os.path.join(save_path, 'control_conservative_' + chr + stra + '.csv'))


### CREATING MOTIF FREE INTERVALS
# path = '/labs/Aguiar/non_bdna/annotations/all_non_b_windows_0_based.csv'
path = '/home/alextu/scratch/nonb_db_preprocessing/nonb_windows/all_non_b_windows_0_based_both_directions.csv'
windows = pd.read_csv(path)
save_path = '/home/alextu/scratch/nonb_db_preprocessing/control_windows'

if not os.path.exists(save_path):
    os.mkdir(save_path)
chr_strand_windows_inputs = []
non_b_groups = windows.groupby(['chr', 'strand']).count()
groups = list(non_b_groups.index)

for gr in groups:
    chr, stra = gr
    # print(chr, stra)
    this_ch_stra = windows[(windows['chr'] == chr) & (windows['strand'] == stra)].reset_index(drop=True)
    print(chr, stra, len(this_ch_stra))
    chr_strand_windows_inputs.append([chr, stra, this_ch_stra, save_path])

pool = Pool(10)
# pool.map(make_control_windows_chr_strand, chr_strand_windows_inputs)
pool.map(make_control_windows_chr_strand_conservative, chr_strand_windows_inputs)
control_windows = pd.DataFrame(columns=['chr', 'feature', 'start', 'end', 'strand', 'motif_len', 'win_start',
                                        'win_end'])
control_files = [name for name in os.listdir(save_path) if '.csv' in name]
for cf in control_files:
    cf_df = pd.read_csv(os.path.join(save_path, cf), index_col=0)
    control_windows = control_windows._append(cf_df, ignore_index=True)
control_windows.to_csv(os.path.join(save_path, 'all_controls_conservative_0_based.csv'))


