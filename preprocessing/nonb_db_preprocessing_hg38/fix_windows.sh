#!/bin/bash
#SBATCH --account=def-sushant
#SBATCH --nodes=1 #specify number of gpus
#SBATCH --mem=32G # MEMORY PER NODE
#SBATCH --time=01:00:00 #hrs:mins:secs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexanderturco1@gmail.com #sends me an email once done
#SBATCH --job-name=fix_windows
#SBATCH --output=fix_windows.o
#SBATCH --error=fix_windows.e

module load python/3.10
source /home/alextu/scratch/compute_windows/ENV/bin/activate

python /home/alextu/scratch/nonb_db_preprocessing/scripts2_gquadsonly/fix_windows.py
