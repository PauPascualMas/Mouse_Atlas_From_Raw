#!/bin/bash

#SBATCH --partition=normal
#SBATCH --job-name=download
#SBATCH --nodes=1 
#SBATCH -o /homes/users/ppascual/scratch/STAR_Align/Scripts/Output_info/slurm.%j.out
#SBATCH -e /homes/users/ppascual/scratch/STAR_Align/Scripts/Output_info/slurm.%j.err


module load wget/1.20.3-GCCcore-9.3.0 

# download Pijuan-Sala .fastq files from ArrayExpress (Downloaded: 10245 files, 927G)

PIJUAN_DATA_DIR='/path/to/downloaded/files/dir'


wget -r -l1 -A.gz -nH --cut-dirs 5 -P $PIJUAN_DATA_DIR https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6967/files/ 

cd $PIJUAN_DATA_DIR

chmod u+x *.gz
