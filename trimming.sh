#!/bin/bash

#SBATCH -J trim # A single job name for the array
#SBATCH -n 1 # Number of cores 
#SBATCH -N 1 # All cores on one machine
#SBATCH -p haswell # Partition 
#SBATCH --mem 20000 # Memory request (40Gb) 
#SBATCH -o /homes/users/ppascual/scratch/STAR_Align/Scripts/Output_info/slurm.%j.out
#SBATCH -e /homes/users/ppascual/scratch/STAR_Align/Scripts/Output_info/slurm.%j.err

### MERGING AND TRIMMING
#Merging all fastqs of each sample (R1, R2 and R3) and creating CB + UMI (R2+R3) read for the later alignment
# Takes cDNA read (R1) and performs adapter trimming and read filtering 

interactive
module load Trim_Galore/0.6.7-GCCcore-10.3.0
module load gzip/1.10-GCCcore-11.2.0 

ALIGN_DATA='/homes/users/ppascual/scratch/STAR_Align/Data'
GENOME_FILE="${ALIGN_DATA}/Genome/genome.fa"
ANNOTATION_FILE="${ALIGN_DATA}/Annotation/mouse_mm10_optimized_v1.gtf.gz"
ALIGN_RESULTS='/homes/users/ppascual/scratch/STAR_Align/Results/Alignment_Pijuan'
READS_DIR="${ALIGN_DATA}/Reads/Pijuan"
UMI_whitelist="${ALIGN_DATA}/BC_whitelist/BC_whitelist_v1/737K-april-2014_rc.txt"

E_stage='E_8.5' 
Sample='Sample_37' 


mkdir $ALIGN_RESULTS/${E_stage}/${Sample}/
cd $READS_DIR/${E_stage}/${Sample}

# merging R1, R2 and R3 files of every batch in sample
cat *_R1*.fastq.gz > merged_R1.fastq.gz
cat *_R2*.fastq.gz > merged_R2.fastq.gz
cat *_R3*.fastq.gz > merged_R3.fastq.gz


### TRIMMING: Trim Galore
#Default settings: Phred score: 20 (--quality); Discard reads shorter than 20 bp (--length)
read_to_trim='merged_R1.fastq.gz'

trim_galore --illumina --output_dir $READS_DIR/${E_stage}/${Sample}/ $READS_DIR/${E_stage}/${Sample}/${read_to_trim}

