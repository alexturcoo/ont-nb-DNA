#!/bin/bash
#SBATCH --account=def-sushant
#SBATCH --nodes=1 #specify number of gpus
#SBATCH --mem=64G # MEMORY PER NODE
#SBATCH --time=011:00:00 #hrs:mins:secs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexanderturco1@gmail.com #sends me an email once done
#SBATCH --job-name=convert_fast5_pod5
#SBATCH --output=convert_fast5_pod5.o
#SBATCH --error=convert_fast5_pod5.e

module load python/3.10 gcc arrow/11
python -c "import pyarrow"

pod5 convert fast5 \
    /home/alextu/scratch/data/hgsvc_data/HG01457/20211022_211012_21-lee-006_PCT0053_2-A5-D5/fast5_pass/*.fast5 \
    --output /home/alextu/scratch/basecalling/hgsvc_pod5s/HG01457/20211022_211012_21-lee-006_PCT0053_2-A5-D5 \
    --one-to-one /home/alextu/scratch/data/hgsvc_data/HG01457/20211022_211012_21-lee-006_PCT0053_2-A5-D5/fast5_pass
