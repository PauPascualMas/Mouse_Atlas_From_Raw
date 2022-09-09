#!/usr/bin/env python
# coding: utf-8

# ### **METADATA CREATION AND OBJECT CONCATENATION**
# #### Script takes the output files from STARsolo and gets them ready for the downstream analysis

import sys, os
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import scipy as sc
from scipy.io import mmread, mminfo, mmwrite 
from scipy import sparse
from scipy.sparse import vstack
from os import listdir
from os.path import isfile, join


#setting working directory
align_results_dir = './Results/STARsolo_Output/'

# #creation of metadata containing cell number, barcode, stage and sample
# stages = next(os.walk(align_results_dir))[1]
# for stage in stages:
#     samples = next(os.walk(align_results_dir + stage + '/'))[1]
#     for sample in samples:
#         lanes = next(os.walk(align_results_dir + stage + '/' + sample + '/'))[1]
#         for lane in lanes:
#             batches = next(os.walk(align_results_dir + stage + '/' + sample + '/' + lane + '/'))[1]
#             batches = [k for k in batches if 'Solo.out' in k]
#             for batch in batches:
#                 features = next(os.walk(align_results_dir + stage + '/' + sample + '/' + lane + '/' + batch))[1]
#                 for feature in features:
#                     if os.path.isfile(align_results_dir + stage + '/' + sample + '/' + lane + '/' + batch + '/' + feature + '/raw/barcodes.tsv'):
#                         bc = pd.read_csv(align_results_dir + stage + '/' + sample + '/' + lane + '/' + batch + '/' + feature + '/raw/barcodes.tsv', header=None)
#                         bc.rename(columns={0: 'barcode'}, inplace=True)
#                         df = bc.drop_duplicates(subset=['barcode'], inplace=False)
#                         df['cell'] = 'cell_' + df.reset_index().index.astype(str)
#                         df['sample'] = sample
#                         df['stage'] = stage
#                         df['stage'] = df['stage'].str.replace('_', '')
#                         df = pd.DataFrame(df)
#                         df.to_csv(align_results_dir + stage + '/' + sample + '/' + lane + '/' + batch + '/' + feature + '/raw/meta.tab', sep="\t", index=False)
#                     else:
#                         print(sample + '_' + lane + ' not found')
#                         break

# # deleting Gene Expression column from features.tsv and metadata concatenation of all samples
# # just on Genefull, as Velocyto has 3 different matrices
# genes = []
# metadata = []
# features = ['GeneFull']
# for stage in stages:
#     samples = next(os.walk(align_results_dir + stage + '/'))[1]
#     for sample in samples:
#         lanes = next(os.walk(align_results_dir + stage + '/' + sample + '/'))[1]
#         for lane in lanes:
#             batches = next(os.walk(align_results_dir + stage + '/' + sample + '/' + lane + '/'))[1]
#             batches = [k for k in batches if 'Solo.out' in k]
#             for batch in batches:
#                 for feature in features:
#                     if os.path.isfile(align_results_dir + stage + '/' + sample + '/' + lane + '/' + batch + '/' + feature + '/raw/barcodes.tsv'):
#                         #Load the data
#                         gene = pd.read_csv(align_results_dir + stage + '/' + sample + '/' + lane + '/' + batch + '/' + feature + '/raw/features.tsv',sep="\t",skiprows=0, header=None)
#                         meta = pd.read_csv(align_results_dir + stage + '/' + sample + '/' + lane + '/' + batch + '/' + feature + '/raw/meta.tab',sep="\t",skiprows=0)
#                         gene = gene.drop(gene.columns[2], axis=1)
#                         metadata.append(meta)
#                     else:
#                         print(sample + '_' + lane + ' file not found')

# df_metadata = pd.concat(metadata, ignore_index=True)

# gene.to_csv(align_results_dir + 'features.tsv', sep='\t', index=False)
# df_metadata.to_csv(align_results_dir + 'metadata_whole.tab', sep='\t', index=False)
# print('Gene Features and Metadata file saved')
                    

# creation of whole matrix object concatenation all samples (high-resource process)
whole_matrix = []
features = ['GeneFull']
stages = next(os.walk(align_results_dir))[1]
for stage in stages:
    samples = next(os.walk(align_results_dir + stage + '/'))[1]
    for sample in samples:
        lanes = next(os.walk(align_results_dir + stage + '/' + sample + '/'))[1]
        for lane in lanes:
            batches = next(os.walk(align_results_dir + stage + '/' + sample + '/' + lane + '/'))[1]
            batches = [k for k in batches if 'Solo.out' in k]
            for batch in batches:
                for feature in features:
                    if os.path.isfile(align_results_dir + stage + '/' + sample + '/' + lane + '/' + batch + '/' + feature + '/raw/barcodes.tsv'):
                        #Load the data
                        mat = mmread(align_results_dir + stage + '/' + sample + '/' + lane + '/' + batch + '/' + feature + '/' + 'raw/matrix.mtx').tocsr().transpose()
                        whole_matrix.append(mat)
                        break
                    else:
                        print(sample + '_' + lane + ' file not found')
                    
whole_matrix = sc.sparse.vstack((whole_matrix))
mmwrite(align_results_dir + 'matrix_whole', whole_matrix)

print('Whole Matrix file saved')




