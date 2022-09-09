#!/bin/bash

#SBATCH -J star # A job-name
#SBATCH -p haswell # Partition 
#SBATCH -n 8 # Number of cores 
#SBATCH --mem 50000 #(GB per node)
#SBATCH -o ./Jobs_info/slurm.%j.out
#SBATCH -e ./Jobs_info/slurm.%j.err

interactive
module load gzip/1.10-GCCcore-11.2.0 
module load STAR/2.7.10a-GCC-11.2.0

### input paths (to modify by user)

SCRIPTS_DIR='/homes/users/ppascual/scratch/Analysis_Mouse_Atlas/Analysis/Pijuan/Alignment' #relative to your path
RAW_FILES_DIR="./Data/Raw_files" # path to where reads are stored
GENOME_INDEX_DIR="./../References/genome_indices"
GENOME_FILE="'./../References/fasta/genome.fa" # path to genome
ANNOTATION_FILE="./../References/genes/genes.gtf" # path to reference annotation
BC_WHITELIST_DIR="./../References/BC_whitelist" # path to barcode list
ALIGN_RESULTS="./Results/STARsolo_Output" 
ALIGN_QC_DIR="./Results/Alignment_QC"

cd $SCRIPTS_DIR

### check input parameters
if [ $# -ne 2 ]
then
    echo "Please, give:"
    echo "1) prefix name for fastq files (name before _R1.fastq.gz and _R2.fastq.gz/CB_UMI.fastq.gz)"   
    echo "2) 10x chemistry version (v1: R1=cDNA, R2=CBC, R3=UMI, I1=Illumina index /// v2 and v3: R1=CBC+UMI, R2=cDNA, R3=Illumina index)"
    exit
fi

prefix=$1
protocol=$2

# 0. GENOME INDEX GENERATION
# If not previously generated, aligner needs to use the reference genome and annotation file in order to generate an INDEX for the later alignment

if [[ ! -d $GENOME_INDEX_DIR ]]
then
    mkdir -p $GENOME_INDEX_DIR
fi

if [[ -f $ANNOTATION_FILE.gz ]]
then
    gunzip $ANNOTATION_FILE.gz
fi


if [[ ! -f $GENOME_INDEX_DIR/genomeParameters.txt ]] 
then
    STAR \
    --runMode genomeGenerate \
    --runThreadN 40 \
    --genomeDir $GENOME_INDEX_DIR \
    --genomeFastaFiles $GENOME_FILE \
    --sjdbGTFfile $ANNOTATION_FILE \
    --genomeSAindexNbases 12 \
    --genomeSAsparseD 3 
fi 


####### 1. PERFORM ALIGNMENT: STARsolo
cd $RAW_FILES_DIR

# Creating CB+UMI file (R2+R3) in case of using v1 protocol, which has different files for CB and UMI
for stage in */; do
    for sample in "$stage"*; do
        sample_id=$(echo $sample | sed 's:.*/::')
        for lane in "$sample"/*; do
            lane_id=$(echo $lane | sed 's:.*/::')

            ### create cb+umi file in case 10x Chemistry is v1
            if [[ $protocol == v1 && ! -e ${lane}/cb_umi_${prefix}_${sample_id}_${lane_id}.fastq.gz ]]
            then
                cd $lane
                cb=$(ls ${prefix}_${sample_id}_${lane_id}_R2.fastq.gz)
                umi=$(ls ${prefix}_${sample_id}_${lane_id}_R3.fastq.gz)
                paste <(zcat $cb) <(zcat $umi) | awk '{if (NR%4==1) {print $1 " " $2} else if (NR%4==3) {print "+"} else {print $1 $2}}' > cb_umi_${prefix}_${sample_id}_${lane_id}.fastq
                gzip cb_umi_${prefix}_${sample_id}_${lane_id}.fastq
                rm -rf cb_umi_${prefix}_${sample_id}_${lane_id}.fastq
            fi
            cd $SCRIPTS_DIR
            cd $RAW_FILES_DIR
        done
    done
done

if [[ $protocol == v1 ]]
then
    echo 'CB + UMI file created successfully for v1 10x protocol'
fi

# Assigning parameters and running STARsolo and FastQC on output
for stage in *; do
    for sample in $stage/*; do
        sample_id=$(echo $sample | sed 's:.*/::')
        for lane in $sample/*; do
            lane_id=$(echo $lane | sed 's:.*/::')

            ### defining  input file = FW (cDNA) and RV file (cb+umi identifiers), may change between protocols
            ### CB and UMI start and length
            if [[ $protocol == v1 ]] 
            then
                # input read files
                r1=$(ls $lane/${prefix}_${sample_id}_${lane_id}_R1_trimmed.fq.gz)
                cb_umi=$(ls $lane/cb_umi_${prefix}_${sample_id}_${lane_id}.fastq.gz)

                # cb and umi star parameters
                cb_len=14
                umi_start=15
                umi_len=10
            fi

            if [[ $protocol == v2 ]] 
            then
                # input read files
                r1=$(ls $lane/${prefix}_${sample_id}_${lane_id}_R2_trimmed.fq.gz) 
                cb_umi=$(ls $lane/${prefix}_${sample_id}_${lane_id}_R1.fastq.gz)  

                # cb and umi star parameters
                cb_len=16
                umi_start=17
                umi_len=10   
            fi

            if [[ $protocol == v3 ]] 
            then
                # input read files
                r1=$(ls ${lane}/${prefix}_${sample_id}_${lane_id}_R2_trimmed.fq.gz)  
                cb_umi=$(ls ${lane}/${prefix}_${sample_id}_${lane_id}_R1.fastq.gz)

                # cb and umi star parameters
                cb_len=16
                umi_start=17
                umi_len=12  
            fi

            # select which barcode whitelist to use depending on protocol
            if [[ $protocol == v1 ]]
            then
                cb_umi_whitelist=$BC_WHITELIST_DIR/737K-april-2014_rc.txt
            fi

            if [[ $protocol == v2 ]] 
            then
                cb_umi_whitelist=$BC_WHITELIST_DIR/737K-august-2016.txt
            fi

            if [[ $protocol == v3 ]] 
            then
                cb_umi_whitelist=$BC_WHITELIST_DIR/3M-february-2018.txt.gz
            fi

            cd $SCRIPTS_DIR

            # runnning STARsolo aligner, parameters to resemble the most as possible as CellRanger
            STAR \
            --genomeDir $GENOME_INDEX_DIR \
            --genomeLoad LoadAndKeep \
            --runThreadN 80 \
            --soloType CB_UMI_Simple \
            --genomeSAsparseD 3 \
            --readFilesCommand zcat \
            --readFilesIn $RAW_FILES_DIR/${r1} .$RAW_FILES_DIR/${cb_umi} \
            --soloCBwhitelist $cb_umi_whitelist \
            --soloCBlen $cb_len \
            --soloUMIstart $umi_start \
            --soloUMIlen $umi_len \
            --outFileNamePrefix .$ALIGN_RESULTS/${lane}/ \
            --clipAdapterType CellRanger4 \
            --outFilterScoreMin 30 \
            --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloUMIfiltering MultiGeneUMI_CR \
            --soloUMIdedup 1MM_CR \
            --soloStrand Forward \
            --outSAMtype BAM Unsorted \
            --outSAMattributes NH HI AS nM GX GN sM sQ sS \
            --soloFeatures GeneFull Velocyto

            rm -rf $ALIGN_RESULTS/E-MTAB-6967.sdrf.txt

            echo 'Alignment for ' $lane ' has been completed!'

            cd $SCRIPTS_DIR
            cd $RAW_FILES_DIR
        done
    done
done


####### 2. QUALITY CONTROL ON STARSOLO RESULTS: MULTIQC


cd $SCRIPTS_DIR
cd $ALIGN_RESULTS

# # QC on the alignment results (summary.csv file)
for stage in *; do
    for sample in $stage/*; do
        sample_id=$(echo $sample | sed 's:.*/::')
        for lane in $sample/*; do
            lane_id=$(echo $lane | sed 's:.*/::')

            # create results folder
            cd $SCRIPTS_DIR
            mkdir -p ./Results/Alignment_QC

            # copying bam files to folder 
            cp $ALIGN_RESULTS/$lane/Log.final.out $ALIGN_QC_DIR/${sample_id}_${lane_id}_Log.final.out
            cd $ALIGN_RESULTS

        done
    done
done


# MultiQC on all FastQC reports, all samples compared at once
module load MultiQC/1.9-foss-2019b-Python-3.7.4 

# cd to where all FastQC files on the STARsolo outputs are stored
cd $SCRIPTS_DIR

# Runnning MultiQC
multiqc Results/Alignment_QC/ -o ../Results/MultiQC/ -n STARsolo_QC_
exit



########################################################

: ' STARsolo configuration regarding each Chemistry version used

10x v1
Whitelist, 737K-april-2014_rc.txt
CB length, 14
UMI start, 15
UMI length, 10 

10X v2
Whitelist, 737K-august-2016.txt
CB length, 16
UMI start, 17
UMI length, 10

10x v3
Whitelist, 3M-Feb_2018_V3.txt
CB length, 16
UMI start, 17
UMI length, 12 '
