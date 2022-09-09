#!/bin/bash

#SBATCH --partition=normal 
#SBATCH --job-name=down_raw
#SBATCH --nodes=1 
#SBATCH -o ./Jobs_info/slurm.%j.out
#SBATCH -e ./Jobs_info/slurm.%j.err


module load wget/1.20.3-GCCcore-9.3.0 

# download Pijuan-Sala .fastq files and metadata from ArrayExpress (Downloaded: 10245 files, 927GB)

RAW_DATA_DIR='./Data/Raw_files'

wget -r -l1 -A.gz -nH --cut-dirs 5 -P $RAW_DATA_DIR https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6967/files/ 
wget -P $RAW_DATA_DIR https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6967/E-MTAB-6967.sdrf.txt
 