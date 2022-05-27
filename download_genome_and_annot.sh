#!/bin/bash

#SBATCH --partition=normal 
#SBATCH --job-name=dwonload
#SBATCH --nodes=1 


module load wget/1.20.3-GCCcore-9.3.0 

# download optimized Genome and Annotation https://www.thepoollab.org/resources

mkdir /path/to/genome
mkdir /path/to/annotation

GEN_DIR=/path/to/genome
ANNOT_DIR=/path/to/annot

wget -P $GEN_DIR/ https://storage.googleapis.com/generecovery/mouse_mm10_optimized_v1.tar.gz

wget -P $ANNOT_DIR/ https://storage.googleapis.com/generecovery/mouse_mm10_optimized_v1.gtf.gz
