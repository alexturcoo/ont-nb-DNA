#!/bin/bash
#SBATCH --account=def-sushant
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=700G # MEMORY PER NODE
#SBATCH --time=48:00:00 #hrs:mins:secs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexanderturco1@gmail.com #sends me an email once done
#SBATCH --job-name=extract_tts
#SBATCH --output=extract_tts.o
#SBATCH --error=extract_tts.e

module load StdEnv/2020
module load python/3.10 gcc arrow/11 parasail

source /home/alextu/scratch/compute_windows/ENV/bin/activate

python /home/alextu/scratch/compute_windows/scripts/extract_tts.py