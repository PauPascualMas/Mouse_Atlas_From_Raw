#!/bin/bash

### bash command: sbatch --array=0-11 trim_galore.sh ###

#SBATCH -J star # A single job name for the array
#SBATCH -n 1 # Number of cores 
#SBATCH -N 1 # All cores on one machine
#SBATCH -p haswell # Partition 
#SBATCH --mem 40000 # Memory request (40Gb) 
#SBATCH -o /homes/users/ppascual/scratch/STAR_Align/Scripts/Output_info/slurm.%j.out
#SBATCH -e /homes/users/ppascual/scratch/STAR_Align/Scripts/Output_info/slurm.%j.err


## TRIMMING  ##
#---------  TRIM GALORE----------

interactive
module load Trim_Galore/0.6.7-GCCcore-10.3.0 
module load gzip/1.10-GCCcore-11.2.0 


## Running STARsolo for selected embryonic stage 

E_stage='mixed' 
Sample='Sample_22' 


READS_DIR=/homes/users/ppascual/scratch/STAR_Align/Data/Reads/Pijuan

  
#Default settings: Phred score: 20 (--quality); Discard reads shorter than 20 bp (--length)' 

trim_galore --illumina --output_dir $READS_DIR/${E_stage}/${Sample}/batch_"${SLURM_ARRAY_TASK_ID}" $READS_DIR/${E_stage}/${Sample}/batch_"${SLURM_ARRAY_TASK_ID}"/merged_R1.fastq.gz
