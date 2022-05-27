#!/bin/bash

#SBATCH --partition=normal 
#SBATCH --job-name=move 
#SBATCH --nodes=1 
#SBATCH -o /homes/users/ppascual/scratch/STAR_Align/Scripts/Output_info/slurm.%j.out
#SBATCH -e /homes/users/ppascual/scratch/STAR_Align/Scripts/Output_info/slurm.%j.err

module load wget/1.20.3-GCCcore-9.3.0 

# download optimized Genome and Annotation https://www.thepoollab.org/resources

DATA_DIR=/homes/users/ppascual/scratch/STAR_Align/Data/


wget -P $DATA_DIR/Genome https://storage.googleapis.com/generecovery/mouse_mm10_optimized_v1.tar.gz

wget -P $DATA_DIR/Annotation https://storage.googleapis.com/generecovery/mouse_mm10_optimized_v1.gtf.gz
