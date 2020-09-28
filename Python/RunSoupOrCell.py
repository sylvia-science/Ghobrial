#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 20:32:21 2020

@author: sujwary
"""

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

import gzip
import shutil


filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = pd.read_excel(filename_metaData)
#metaData = metaData[metaData['Run']== 1]
metaData = metaData[metaData['Sample Type']== 'PBMC']


base = '/disk2/Projects/EloRD/Data/Bam/PBMC/'
os.chdir(base)

#sample_name = 'NBM7CD138N'
for i in range(0,(metaData.shape[0])):
    sample_name = metaData['Sample'].iloc[i]
    print(sample_name)
    
    output = base + sample_name + '_demux_data/'


    command = 'singularity exec /disk2/Projects/Code/souporcell.sif souporcell_pipeline.py '
    command = command + '-i ' + sample_name + '.bam '
    command = command + '-b ' + sample_name + '_out_cell_barcodes.csv '
    command = command + '-f ' + base + 'ReferenceData/refdata-cellranger-GRCh38-1.2.0.fasta '
    command = command +  '-t 8 '
    command = command + '-o ' + output + ' '
    command = command + '-k 8 '
    #command = command + '--ignore True'
    
    if not (path.exists(output)):
        print(command)
        os.system(command)


for i in range(0,(metaData.shape[0])):
    sample_name = metaData['Sample'].iloc[i]
    print(sample_name)
    
    output = base + sample_name + '_demux_data/'
    command = 'perl  /disk2/Projects/Code/mskcc-vcf2maf-47c4a18/vcf2maf.pl '
    command = command + '--input-vcf ' + output + 'souporcell_merged_sorted_vcf.vcf '
    command = command + '--output-maf ' +  output + sample_name + '.vep.maf '
    command = command + '--ref-fasta /disk2/Projects/EloRD/Data/ReferenceData/refdata-cellranger-GRCh38-1.2.0.fasta '
    command = command + 'perl  /disk2/Projects/Code/mskcc-vcf2maf-47c4a18/vcf2maf.pl '
    command = command + '--filter-vcf 0 --vep-path /disk2/Projects/Code/ensembl-vep/'
    
    file = base + sample_name + '_demux_data/' + 'souporcell_merged_sorted_vcf.vcf'
    if path.exists(file + '.gz') and not path.exists(file):
        with gzip.open(file + '.gz', 'rb') as f_in:
            with open(file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            
    #runCommand = path.exists(base + sample_name + '.vcf') and not (path.exists(base + sample_name + '.vep.maf'))
    if path.exists(base + sample_name + '_demux_data/' + 'souporcell_merged_sorted_vcf.vcf') and not path.exists(base + sample_name + '_demux_data_harmony/' + sample_name + '.vep.maf'):
        print(command)
        os.system(command)
        


