import sys
import os
import pandas as pd
import numpy as np
from subprocess import call

### This script finds reads that fall within regions (taken from GoFae-DND)
### (regions are start of first non-B window to end of last non-B window on EACH CHROMOSOME)

save_path = '/home/alextu/scratch/compute_windows/reads_in_nonb_range/HG00733/20220223_220216_21-lee-006_PC24B149_3F'
bam_path = '/home/alextu/scratch/basecalling/merged_bams/HG00733/20220223_220216_21-lee-006_PC24B149_3F_merged.bam'
# path = '/mnt/research/aguiarlab/proj/nonBDNA/data/chr_start_end_queries_1_based.csv' # beagle
path = '/home/alextu/scratch/nonb_db_preprocessing/chromosome_start_end_queries/chr_start_end_queries_1_based.csv'
df = pd.read_csv(path, index_col=False)
sh_path = '/home/alextu/scratch/compute_windows/executables/filter_query_samtools.sh'
to_print = '#!/bin/bash\n' \
            '#SBATCH --account=def-sushant\n' \
            '#SBATCH --nodes=1\n' \
            '#SBATCH --cpus-per-task=8\n' \
            '#SBATCH --mem=64G\n' \
            '#SBATCH --time=10:00:00\n' \
            '#SBATCH --mail-type=ALL\n' \
            '#SBATCH --mail-user=alexanderturco1@gmail.com\n' \
            '#SBATCH --job-name=filter_query_samtools\n' \
            '#SBATCH --output=filter_query_samtools.o\n' \
            '#SBATCH --error=filter_query_samtools.e\n' \
            '\n' \
            'echo `hostname`\n' \
            'module load samtools\n' \
            'module load bedtools\n' \
            '\n'
for i in range(len(df)):
    chr = df.loc[i, 'chromosome']
    strand = df.loc[i, 'strand']
    chrno = df.loc[i, 'chr']
    start = df.loc[i, 'start']
    end = df.loc[i, 'end']
    name = chrno + strand
    query = chr + ':' + str(start) + '-' + str(end)
    comm1 = ''
    if strand == '+':
        comm1 = 'samtools view {} {} -F 0xF14 -o {}{}.bam;'.format(bam_path, query, save_path, name)
    elif strand == '-':
        comm1 = 'samtools view {} {} -F 0xF04 -f 0x10 -o {}{}.bam;'.format(bam_path, query, save_path, name)
    comm2 = 'bedtools bamtobed -i {}{}.bam > {}{}.bed;'.format(save_path, name, save_path, name)
    print(comm1)
    print(comm2)
    to_print += comm1 + '\n'
    to_print += comm2 + '\n'
with open(sh_path, 'w') as f:
    f.write(to_print)
# run the command
call(['sbatch', sh_path])