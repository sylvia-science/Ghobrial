#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 15:36:59 2020

@author: sujwary
"""

import os
import pandas as pd
import os.path
from os import path

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = pd.read_excel(filename_metaData)
metaData = metaData[metaData['Run']== 1]

base = '/disk2/Projects/EloRD/Data/Bam/'
os.chdir(base)




for i in range(0,(metaData.shape[0])):
    sample_name = metaData['Sample'].iloc[i]
    print(sample_name)
    
    command = 'freebayes-parallel <(fasta_generate_regions.py /disk2/Projects/EloRD/Data/ReferenceData/refdata-cellranger-GRCh38-1.2.0.fasta.fai  100000) 36 '
    command = command + '-f /disk2/Projects/EloRD/Data/ReferenceData/refdata-cellranger-GRCh38-1.2.0.fasta '
    command = command + sample_name + '.bam >'+ sample_name +'.vcf'

    #runCommand = path.exists(base + sample_name + '.vcf') and not (path.exists(base + sample_name + '.vep.maf'))
    if path.exists(base + sample_name + '.bam') and not (path.exists(base + sample_name + '.vcf')):
        print(command)
        os.system(command)


for i in range(0,(metaData.shape[0])):
    sample_name = metaData['Sample'].iloc[i]
    print(sample_name)
    
    command = 'perl  /disk2/Projects/Code/mskcc-vcf2maf-47c4a18/vcf2maf.pl '
    command = command + '--input-vcf ' + sample_name + '.vcf '
    command = command + '--output-maf ' + sample_name + '.vep.maf '
    command = command + '--ref-fasta /disk2/Projects/EloRD/Data/ReferenceData/refdata-cellranger-GRCh38-1.2.0.fasta '
    command = command + 'perl  /disk2/Projects/Code/mskcc-vcf2maf-47c4a18/vcf2maf.pl '
    command = command + '--filter-vcf 0 --vep-path /disk2/Projects/Code/ensembl-vep/'
    #runCommand = path.exists(base + sample_name + '.vcf') and not (path.exists(base + sample_name + '.vep.maf'))
    if path.exists(base + sample_name + '.vcf'): #and not (path.exists(base + sample_name + '.vep.maf')):
        print(command)
        os.system(command)
        
        
for i in range(0,(metaData.shape[0])):
    sample_name = metaData['Sample'].iloc[i]
    print(sample_name)
    
    command = 'perl  /disk2/Projects/Code/mskcc-vcf2maf-47c4a18/vcf2maf.pl '
    command = command + '--input-vcf ' + sample_name + '.vcf '
    command = command + '--output-maf ' + sample_name + '.vep.maf '
    command = command + '--ref-fasta /disk2/Projects/EloRD/Data/ReferenceData/refdata-cellranger-GRCh38-1.2.0.fasta '
    command = command + 'perl  /disk2/Projects/Code/mskcc-vcf2maf-47c4a18/vcf2maf.pl '
    command = command + '--filter-vcf 0 --vep-path /disk2/Projects/Code/ensembl-vep/'
    #runCommand = path.exists(base + sample_name + '.vcf') and not (path.exists(base + sample_name + '.vep.maf'))
    if path.exists(base + sample_name + '.vcf'): #and not (path.exists(base + sample_name + '.vep.maf')):
        print(command)
        os.system(command)


