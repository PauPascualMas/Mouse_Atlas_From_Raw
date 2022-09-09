#!/bin/bash

#SBATCH --partition=normal 
#SBATCH --job-name=download 
#SBATCH --nodes=1 
#SBATCH -o ./Jobs_info/slurm.%j.out
#SBATCH -e ./Jobs_info/slurm.%j.err

module load wget/1.20.3-GCCcore-9.3.0 

if [[ ! -e $BC_WHITELIST_DIR ]]
then
    mkdir -P ./References/BC_whitelist
fi

REF_DIR='./References'
BC_WHITELIST_DIR='./References/BC_whitelist'

# download optimized Genome and Annotation https://www.thepoollab.org/resources
wget -P $REF_DIR https://storage.googleapis.com/generecovery/mouse_mm10_optimized_v1.tar.gz

cd $REF_DIR 
tar xvfz mouse_mm10_optimized_v1.tar.gz
mv home/allan-hermann/Documents/mouse_mm10_optimized_v1/* .
rm -rf home/

# download BC whitelist used (it changes depending on the 10x protocol version used)

wget https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt -O $BC_WHITELIST_DIR/737K-april-2014_rc.txt
wget https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/737K-august-2016.txt -O $BC_WHITELIST_DIR/737K-august-2016.txt
wget https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz -O $BC_WHITELIST_DIR/3M-february-2018.txt.gz
        
