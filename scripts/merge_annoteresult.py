#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 16:16:43 2021

@author: cavu
"""

import os
import sys
import numpy as np

annorestdir = sys.argv[1] 
#annorestdir = "/home/cavu/DATA/annotation/@1RTP1提取/batch_annotation_of_RTP1/annotation_result3"

annofiles = os.listdir(annorestdir)
annofiles.sort()

anno_amino = [i for i in annofiles if "_aminoseq" in i]
anno_cds = [i for i in annofiles if "_cds" in i]


def anno_merge(files):
    anno_list = []
    for file in files:
        int_file = os.path.join(annorestdir, file)
        with open(int_file) as f:
            linelist = [line for line in f ]
        name = ">" + file + "\n"
        linelist = ''.join(linelist)
        line = linelist.strip().split('\n')
        seq = ""
        for sub_line in line:
            sub_seq = name + ''.join(sub_line)+ "\n"
            seq = seq + sub_seq
        #all = name + seq
        anno_list.append(seq)
    return anno_list
        
        
amin_anno =  anno_merge(anno_amino)  
amin_cds=  anno_merge(anno_cds) 

def exclude_somespec(all_merge,cds):
    species_list = []
    seq_list = []
    intact_infor = []
    intact_cda = []
    for seq,cdsseq in zip(all_merge,cds):
        name = seq.split("txt")[0]
        species = name.split("-")[0]
        seq_contact = seq.strip().strip().split("txt")[1]
        if len(seq_contact) != 0:
            if not(species in species_list):
                species_list.append(species)
                seq_list.append(seq_contact)
                intact_infor.append(seq)
                intact_cda.append(cdsseq)
                '''
            elif species in species_list:
                loc = species_list.index(species)
                if seq_list[loc] != seq_contact:
                    species_list.append(species)
                    seq_list.append(seq_contact)
                    intact_infor.append(seq)
                    intact_cda.append(cdsseq)
                    '''
    return species_list,seq_list,intact_infor,intact_cda
                
                
a,b,amin_anno_exclude,cds_anno_exclude= exclude_somespec(amin_anno,amin_cds)        
        
        

with open(annorestdir + "/mergeamino2.txt","w") as f:
    f.writelines(amin_anno_exclude)

with open(annorestdir + "/mergecds2.txt","w") as f:
    f.writelines(cds_anno_exclude)