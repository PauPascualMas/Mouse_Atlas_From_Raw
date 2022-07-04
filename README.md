# Mouse_Atlas_From_Raw (v0.1)
This pipeline starts by downloading the single-cell raw reads sequencing data (fastq files) from [Pijuan-Sala et al 2019](https://www.nature.com/articles/s41586-019-0933-9) in order to align these reads against a reference genome and annotation to get an assembled count matrix which can be used in downstream analysis.

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

To download the genome and annotation reference used later in the alignment, run `genome_and_annotation_download.sh`
Optimized references: [Enhanced Gene Recovery](https://www.biorxiv.org/content/10.1101/2022.04.26.489449v1.full) and [The Pool Lab](https://www.thepoollab.org/resources)

## 1. Stage Classification and Merging Lanes

As all files from all the different stages are downloaded together, `stage_sample_clasifyer.ipynb` creates folders for the diifferent embryonic stages and samples and moves the read files according to the experiment metadata.

## 2. Cell Barcode Filtering (not finished)

Previous to the Trimming process and Alignment, a preliminary exploration and barcode complexity is made with. By this, files with low complexity or not associated to UMIs will be discarded to optimize the downstream pipeline.

Anna's pipeline: `get_cellBC_fastqPijuansala.ipynb` and `concatenatorPijuansala.ipynb`

## 3. Merging and Trimming: Trim Galore

First, `ḿerging_and_trimming.sh` merges all files of each type (R1, R2, and R3) into a **merged_R-.fastq.gz** file containing all reads for the selected sample (stage and sample have to be input manually as variables: _E_stage_ and _Sample_).

Then, [Trim Galore](https://github.com/FelixKrueger/TrimGalore) (with default settings) performs the trimming process on the cDNA read (R1), whose output is saved as **merged_R1_trimmed.fq.gz** in the selected stage/sample.


## Alignment with STARsolo

Regarding the alignment, `alignment.sh` firstly creates the CB+UMI file needed as the second read input before running [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) (STAR v 2.7.10a). Parameters are set to resemble at its finest CellRanger results, but much more efficiently in terms of time and resources.

The alignment output consists of: 
**matrix.mtx**: sparse count matrix 

**barcodes.tsv**: sorted list of the barcodes

**features.tsv**: sorted list of gene names

Moreover, for Velocyto results matrix.mtx is divided into **spliced, unspliced** and **ambiguous** counts.
