# Mouse_Atlas_From_Raw (v0.1)
This pipeline starts by downloading the single-cell raw reads sequencing data (.fastq) from [Pijuan-Sala et al 2019](https://www.nature.com/articles/s41586-019-0933-9) in order to obtain an assembled count matrix which can be used in downstream analysis.

## 0. Download data
Data can be downloaded by executing `download_pijuan.sh`.

This script will download all raw reads file to the assigned directory together with the sequencing experiment metadata.

Each read is divided into four files as 10x Genomics v.1 Chemistry was used for the sequencing (R1, R2, R3 and I1). For instance:

> 23216_2_TCAGTCAA_S108_L002_R1_001.fastq.gz
> 
> 23216_2_TCAGTCAA_S108_L002_R2_001.fastq.gz
> 
> 23216_2_TCAGTCAA_S108_L002_R3_001.fastq.gz
> 
> 23216_2_TCAGTCAA_S108_L002_I1_001.fastq.gz

**R1** = transcript (98bp); 
**R2** = 14 barcode; 
**R3** = 10bp UMI;
**I1** = sample index; 

## 1. Stage Classification and Merging Lanes
As all files from all the different stages are downloaded together, `stage_clasiifyer.ipynb` moves each file to the correspoding stage folder according to the experiment metadata.

Lanes within each embryonic stage are merged with `merging_lanes_stage.sh`.

## 2. Cell Barcode Filtering (not finished)
Previous to the Trimming process and Alignment, a preliminary exploration and barcode complexity is made with. By this, files with low complexity or not associated to UMIs will be discarded to optimize the downstream pipeline.

Anna's pipeline: `get_cellBC_fastqPijuansala.ipynb`and`concatenatorPijuansala.ipynb`

## Trimming: Trim-Galore (may cause problems with alignment)
Fastq adapter trimming process with [Trim Galore](https://github.com/FelixKrueger/TrimGalore)

## Alignment with STARsolo

Alignment with [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md).

Optimized references: [Enhanced Gene Recovery](https://www.biorxiv.org/content/10.1101/2022.04.26.489449v1.full) and [The Pool Lab](https://www.thepoollab.org/resources)

