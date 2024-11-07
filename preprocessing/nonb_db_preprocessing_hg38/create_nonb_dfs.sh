#!/bin/bash
#SBATCH --account=def-sushant
#SBATCH --nodes=1 #specify number of gpus
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G # MEMORY PER NODE
#SBATCH --time=01:00:00 #hrs:mins:secs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexanderturco1@gmail.com #sends me an email once done
#SBATCH --job-name=create_nonb_dfs
#SBATCH --output=create_nonb_dfs.o
#SBATCH --error=create_nonb_dfs.e

module load python/3.10
source /home/alextu/scratch/compute_windows/ENV/bin/activate

python /home/alextu/scratch/chm13_nonb_db_preprocessing/scripts/create_nonb_dfs.py