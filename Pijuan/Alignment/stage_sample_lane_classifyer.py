#!/usr/bin/env python
# coding: utf-8

# ### **FILE CLASSIFICATION INTO STAGE AND SAMPLE**
# #### Creates a directory classification of fastq files into the corresponding **embryonic stage**, **sample** and **sequencing lane**.
# #### Raw fastq files and the sequencing metadata are downloaded from *Pijuan-Sala et al. (2019)*

# In[5]:


import pandas as pd
import os 
import shutil
from os import listdir
from os.path import isfile, join


# #### Listing *fastq files* and loading *metadata*
#path to where files have been downloaded
path_to_files= './Data/Raw_files/'
path_to_metadata='./Data/'

#defining the different stages
stages=['E_6.5','E_6.75','E_7.0','E_7.25','E_7.5','E_7.75','E_8.0','E_8.25','E_8.5', 'mixed']

#list of files in the directory
onlyfiles = [f for f in listdir(path_to_files) if isfile(join(path_to_files, f))]
onlyfiles = [k for k in onlyfiles if 'fastq.gz' in k]

#importing metadata file
met = pd.read_csv(path_to_metadata + 'E-MTAB-6967.sdrf.txt', sep='\t')
met['Source Name']= met.replace(' ', '_', regex=True)
met

#add column containing lane information
file=met['Comment[index1 file]']
lane = [i.split("_")[4] for i in file.values] 
met['Lane']=lane
met

#filename reference from metadata
names = met['Assay Name'].to_list()

# #### Folder creation for **stage / sample / lane** and file classification
for stage in stages:
    stage_dir = os.makedirs(path_to_files + stage + '/', exist_ok=True)
    print('Stage ' + stage + ' folder created')

#classification of fastq.gz files into the different folders
for filename in os.listdir(path_to_files):
    if filename.endswith('.fastq.gz'):
        for filename in onlyfiles:
            for assay in names:
                if filename.startswith(assay):
                    stage = met.loc[met['Assay Name'] == assay, 'Characteristics[developmental stage]'].item()
                    if stage == 'embryonic day 6.5':
                            os.replace(path_to_files + filename, path_to_files + "E_6.5/" + filename)
                    elif stage == 'embryonic day 6.75':
                            os.replace(path_to_files + filename, path_to_files + "E_6.75/" + filename)
                    elif stage == 'embryonic day 7.0':
                            os.replace(path_to_files + filename, path_to_files + "E_7.0/" + filename)
                    elif stage == 'embryonic day 7.25':
                            os.replace(path_to_files + filename, path_to_files + "E_7.25/" + filename)
                    elif stage == 'embryonic day 7.5':
                            os.replace(path_to_files + filename, path_to_files + "E_7.5/" + filename)
                    elif stage == 'embryonic day 7.75':
                            os.replace(path_to_files + filename, path_to_files + "E_7.75/" + filename)
                    elif stage == 'embryonic day 8.0':
                            os.replace(path_to_files + filename, path_to_files + "E_8.0/" + filename)
                    elif stage == 'embryonic day 8.25':
                            os.replace(path_to_files + filename, path_to_files + "E_8.25/" + filename)
                    elif stage == 'embryonic day 8.5':
                            os.replace(path_to_files + filename, path_to_files + "E_8.5/" + filename)
                    elif stage == 'mixed':
                            os.replace(path_to_files + filename, path_to_files + "mixed/" + filename)
                    print('Files succesfully classified into embryonic stages!')   
                    break
else:
    print('Fastq files not found. Maybe classification has already been done.')                                          
    break                              

#moving files according to sample 
for filename in os.listdir(path_to_files + stage):
    if filename.endswith('.fastq.gz'):
        for filename in onlyfiles:
            for stage in stages:
                files_stage=[f for f in listdir(path_to_files + stage) if isfile(join(path_to_files + stage, f))]
                for filename in files_stage:
                    for assay in names:
                        if filename.startswith(assay):
                            sample = met.loc[met['Assay Name'] == assay, 'Source Name'].item()
                            os.makedirs(path_to_files + stage + '/' + sample + '/', exist_ok=True)
                            shutil.move(path_to_files + stage + '/' + filename, path_to_files + stage + '/' + sample + '/' + filename)  
    else:
        print('Fastq files not found. Maybe classification has already been done.')
        break

# moving files according to lane
for filename in os.listdir(path_to_files + 'E_6.5/Sample_1/'): #as a exampple
    if filename.endswith('.fastq.gz'):
        for stage in stages:
            samples = next(os.walk(path_to_files + stage))[1]
            for sample in samples:
                files_sample=[f for f in listdir(path_to_files + stage + '/' + sample) if isfile(join(path_to_files + stage + '/' + sample, f))]
                for filename in files_sample:
                    for assay in names:
                        if filename.startswith(assay):
                            lane = met.loc[met['Assay Name'] == assay, 'Lane'].item()
                            os.makedirs(path_to_files + stage + '/' + sample + '/' + lane + '/', exist_ok=True)
                            shutil.move(path_to_files + stage + '/' + sample + '/' + filename, path_to_files + stage + '/' + sample + '/' + lane + '/' + filename)
    else:
        print('Fastq files not found. Maybe classification has already been done.')
       

