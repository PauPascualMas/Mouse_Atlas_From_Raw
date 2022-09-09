#!/bin/bash

#SBATCH -J qc # A job-name
#SBATCH -N 3 # Number of nodes 
#SBATCH -n 16 # number of cores
#SBATCH -p haswell # Partition 
#SBATCH --mem 8000 # Memory request (8Gb) 
#SBATCH -o ./Jobs_info/slurm.%j.out
#SBATCH -e ./Jobs_info/slurm.%j.err

module load MultiQC/1.9-foss-2019b-Python-3.7.4 
module load Trim_Galore/0.6.7-GCCcore-10.3.0

prefix=$1
protocol=$2

SCRIPTS_DIR='/homes/users/ppascual/scratch/Analysis_Mouse_Atlas/Analysis/Pijuan/Alignment'
RAW_FILES_DIR='./Data/Raw_files/'
FASTQC_OUT='./Results/FastQC_Output' 
TRIMMED_FASTQC_OUT='./Results/Trimmed_FastQC_Output'

# MERGE + QC /// TRIMMING + QC 

#1. MERGE + QC (FastQC + MultiQC)
#Pipeline accesses every stage/sample/lane directory and merges its files. 
#Then, FastQC is performed for every merged cDNA file. Finally, MultiQC blends all FastQC reports.

module load FastQC/0.11.9-Java-11 

cd $RAW_FILES_DIR

for stage in */; do
	for sample in "$stage"*; do
		sample_id=$(echo $sample | sed 's:.*/::')
		for lane in "$sample"/*; do
		lane_id=$(echo $lane | sed 's:.*/::')
			for file in "$lane"/*; do
				# merging all R1, R2 and R3 files of each sample and lane
				if [[ ! -e $lane/${prefix}_${sample_id}_${lane_id}_R1.fastq.gz ]]
				then
		    		cat $lane/2*_R1*.fastq.gz > $lane/${prefix}_${sample_id}_${lane_id}_R1.fastq.gz 
		    		cat $lane/2*_R2*.fastq.gz > $lane/${prefix}_${sample_id}_${lane_id}_R2.fastq.gz
		    		cat $lane/2*_R3*.fastq.gz > $lane/${prefix}_${sample_id}_${lane_id}_R3.fastq.gz	
				fi
			done
            #PERFORMING FASTQC FOR EVERY MERGED LANE
            #create output directory
            cd $SCRIPTS_DIR

            if [[ ! -d $FASTQC_OUT ]]; then
                mkdir $FASTQC_OUT
            fi

            #setting cDNA read for QC (depends on 10x protocol)
            if [[ $protocol == v1 ]]; then
                file2qc=${lane}/${prefix}_${sample_id}_${lane_id}_R1.fastq.gz
            fi

            if [[ $protocol == v2 && $protocol == v3 ]]; then
                file2qc=${lane}/${prefix}_${sample_id}_${lane_id}_R2.fastq.gz
            fi

            fastqc -o $FASTQC_OUT $file2qc
            exit
		done
	done
done

#RUNNING MULTIQC ON ALL FASTQC REPORTS
module load MultiQC/1.9-foss-2019b-Python-3.7.4 

mkdir -P ../Results/MultiQC

multiqc Results/FastQC_Output/ -o ../Results/MultiQC/ -n Untrimmed_ 
exit

module unload MultiQC/1.9-foss-2019b-Python-3.7.4 

# 2. TRIMMING + QC (Trim Galore and FastQC + MultiQC)
# cDNA fastq file is trimmed using Trim Galore default parameters 
# Once trimmed, FastQC and MultiQC is performed.

module load Trim_Galore/0.6.7-GCCcore-10.3.0

cd $RAW_FILES_DIR

for stage in */; do
	e_stage=$(echo $stage | sed 's:.*/::')
	for sample in "$stage"*; do
		sample_id=$(echo $sample | sed 's:.*/::')
		for lane in "$sample"/*; do
			lane_id=$(echo $lane | sed 's:.*/::')

			#PERFORMING FASTQC FOR EVERY MERGED LANE
			#create output directory
			cd $SCRIPTS_DIR

			if [[ ! -d $TRIMMED_FASTQC_OUT ]]; then
				mkdir $TRIMMED_FASTQC_OUT
			fi

			#setting cDNA read for Trimming and QC (depends on 10x protocol)
			cd $RAW_FILES_DIR
			if [[ $protocol == v1 ]]; then
				file2trim=${lane}/${prefix}_${sample_id}_${lane_id}_R1.fastq.gz		
			fi

			if [[ $protocol == v2 && $protocol == v3 ]]; then
				file2trim=${lane}/${prefix}_${sample_id}_${lane_id}_R2.fastq.gz
			fi

			# TRIMMING: Trim Galore (default parameters)
			trim_galore --illumina --output_dir $lane $file2trim
			exit
			# Defining cDNA file depending on protocol
			if [[ $protocol == v1 ]]; then
				file2trim_qc=${lane}/${prefix}_${sample_id}_${lane_id}_R1_trimmed.fq.gz
			fi

			if [[ $protocol == v2 && $protocol == v3 ]]; then
				file2trim_qc=${lane}/${prefix}_${sample_id}_${lane_id}_R2_trimmed.fq.gz
			fi

			#RUNNING FASTQC ON TRIMMED FASTQ FILE (CDNA)

			fastqc -o ../../$TRIMMED_FASTQC_OUT $file2trim_qc
			exit
		done
	done
done

module unload Trim_Galore/0.6.7-GCCcore-10.3.0

#RUNNING MULTIQC ON ALL FASTQC REPORTS
module load MultiQC/1.9-foss-2019b-Python-3.7.4 

mkdir -P ../Results/MultiQC

multiqc Results/Trimmed_FastQC_Output/ -o ../Results/MultiQC/ -n Trimmed_ 
exit

module unload MultiQC/1.9-foss-2019b-Python-3.7.4 
module unload FastQC/0.11.9-Java-11 
