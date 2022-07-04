#!/bin/bash

#SBATCH -J star # A single job name for the array
#SBATCH -n 1 # Number of cores 
#SBATCH -N 1 # All cores on one machine
#SBATCH -p haswell # Partition 
#SBATCH --mem 40000 # Memory request (40Gb) 
#SBATCH -o /homes/users/ppascual/scratch/STAR_Align/Scripts/Output_info/slurm.%j.out
#SBATCH -e /homes/users/ppascual/scratch/STAR_Align/Scripts/Output_info/slurm.%j.err


## ALIGNMENT AND CREATION OF RAW COUNTS MATRIX ##
#---------   STAR 2.7.10a (STARsolo)----------

# STARsolo takes cDNA (R1) read and CB+UMI (R2+R3) read to perform the alignment. 
# Genome (.fa) and annotation file (.gtf) are also provided. 
# Barcode whitelist (.txt) can be provided as guide of which BC are real.

interactive
module load  STAR/2.7.10a-GCC-11.2.0
module load gzip/1.10-GCCcore-11.2.0 

ALIGN_DATA='/homes/users/ppascual/scratch/STAR_Align/Data'
GENOME_FILE="${ALIGN_DATA}/Genome/genome.fa"
ANNOTATION_FILE="${ALIGN_DATA}/Annotation/mouse_mm10_optimized_v1.gtf.gz"
ALIGN_RESULTS='/homes/users/ppascual/scratch/STAR_Align/Results/Alignment_Pijuan'
READS_DIR="${ALIGN_DATA}/Reads/Pijuan"
UMI_whitelist="${ALIGN_DATA}/BC_whitelist/BC_whitelist_v1/737K-april-2014_rc.txt"

## Genome Index generation is preliminary in the alignment by STAR
## In our case, optimized genome and annotation used, as well as index are directly obtained from https://www.thepoollab.org/resources

## If willing to generate your own index, use the following command:
: '
STAR \
--runMode genomeGenerate \
--runThreadN 8 \
--genomeDir /homes/users/ppascual/scratch/STAR_Align/Data/Genome/ \
--genomeFastaFiles $GENOME_FILE \
--sjdbGTFfile $ANNOTATION_FILE \
--genomeSAindexNbases 12 \
--limitGenomeGenerateRAM 10000000000'


### STARsolo

#Select stage and sample to align

E_stage='E_8.5' 
Sample='Sample_37' 

cd $READS_DIR/${E_stage}/${Sample}

#creating fastq file merging cb and umi (in 10x Chemistry v1 cb and umi are in separate files), as required for the alignment

rm -rf CB_UMI_merged.fastq.gz #remove preexisting file, just in case
paste <(zcat merged_R2.fastq.gz) <(zcat merged_R3.fastq.gz) | awk '{if (NR%4==1) {print $1 " " $2} else if (NR%4==3) {print "+"} else {print $1 $2}}' > CB_UMI_merged.fastq
gzip CB_UMI_merged.fastq
chmod u+x CB_UMI_merged.fastq.gz
rm -rf CB_UMI_merged.fastq

### input file = FW / define RV file

fw='merged_R1_trimmed.fq.gz'
cb_umi='CB_UMI_merged.fastq.gz'

cd 

STAR \
--genomeDir $ALIGN_DATA/Genome \
--runThreadN 30 \
--soloType CB_UMI_Simple \
--genomeSAsparseD 3 \
--readFilesCommand zcat \
--readFilesIn $READS_DIR/${E_stage}/${Sample}/merged_R1_trimmed.fq.gz $READS_DIR/${E_stage}/${Sample}/CB_UMI_merged.fastq.gz \
--soloCBwhitelist $UMI_whitelist \
--soloCBlen 14 \
--soloUMIstart 15 \
--soloUMIlen 10 \
--clipAdapterType CellRanger4 \
--outFilterScoreMin 30 \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloUMIfiltering MultiGeneUMI_CR \
--soloUMIdedup 1MM_CR \
--outFileNamePrefix $ALIGN_RESULTS/${E_stage}/${Sample} \
--outSAMtype BAM SortedByCoordinate \
--bamRemoveDuplicatesType UniqueIdenticalNotMulti \
--outSAMattributes NH HI AS nM GX CB UB GN sM sQ \
--soloFeatures GeneFull Velocyto


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

