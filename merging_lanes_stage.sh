#!/bin/bash


#SBATCH --partition=haswell
#SBATCH --job-name=merge
#SBATCH --nodes=1
#SBATCH --mem 15000

## merging all .fastq.gz file type for every stage

PATH_TO_FILES='PATH/TO/YOUR/READ_FILES'

#E-stage 6.5
cd $PATH_TO_FILES/E_6.5

cat *_R1*.fastq.gz > merged_R1_E_6.5.fastq.gz
cat *_R2*.fastq.gz > merged_R2_E_6.5.fastq.gz
cat *_R3*.fastq.gz > merged_R3_E_6.5.fastq.gz


#E-stage 6.75
cd $PATH_TO_FILES/E_6.75

cat *_R1*.fastq.gz > merged_R1_E_6.75.fastq.gz
cat *_R2*.fastq.gz > merged_R2_E_6.75.fastq.gz
cat *_R3*.fastq.gz > merged_R3_E_6.75.fastq.gz

#E-stage 7.0
cd $PATH_TO_FILES/E_7.0

cat *_R1*.fastq.gz > merged_R1_E_7.0.fastq.gz
cat *_R2*.fastq.gz > merged_R2_E_7.0.fastq.gz
cat *_R3*.fastq.gz > merged_R3_E_7.0.fastq.gz


#E-stage 7.25
cd $PATH_TO_FILES/E_7.25

cat *_R1*.fastq.gz > merged_R1_E_7.25.fastq.gz
cat *_R2*.fastq.gz > merged_R2_E_7.25.fastq.gz
cat *_R3*.fastq.gz > merged_R3_E_7.25.fastq.gz


#E-stage 7.5
cd $PATH_TO_FILES/E_7.5

cat *_R1*.fastq.gz > merged_R1_E_7.5.fastq.gz
cat *_R2*.fastq.gz > merged_R2_E_7.5.fastq.gz
cat *_R3*.fastq.gz > merged_R3_E_7.5.fastq.gz


#E-stage 7.75
cd $PATH_TO_FILES/E_7.75

cat *_R1*.fastq.gz > merged_R1_E_7.75.fastq.gz
cat *_R2*.fastq.gz > merged_R2_E_7.75.fastq.gz
cat *_R3*.fastq.gz > merged_R3_E_7.75.fastq.gz


#E-stage 8.0
cd $PATH_TO_FILES/E_8.0

cat *_R1*.fastq.gz > merged_R1_E_8.0.fastq.gz
cat *_R2*.fastq.gz > merged_R2_E_8.0.fastq.gz
cat *_R3*.fastq.gz > merged_R3_E_8.0.fastq.gz


#E-stage 8.25
cd $PATH_TO_FILES/E_8.25

cat *_R1*.fastq.gz > merged_R1_E_8.25.fastq.gz
cat *_R2*.fastq.gz > merged_R2_E_8.25.fastq.gz
cat *_R3*.fastq.gz > merged_R3_E_8.25.fastq.gz


#E-stage 8.5
cd $PATH_TO_FILES/E_8.5

cat *_R1*.fastq.gz > merged_R1_E_8.5.fastq.gz
cat *_R2*.fastq.gz > merged_R2_E_8.5.fastq.gz
cat *_R3*.fastq.gz > merged_R3_E_8.5.fastq.gz
