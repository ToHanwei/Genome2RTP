#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:37:00 2021

@author: cavu
"""

#!coding:utf-8

import os
import sys
import logging
import multiprocessing as mp
import pandas as pd
hmmdir = sys.argv[1]  #参数存放hmm的路径
genomedir = sys.argv[2]
nhmmerout = sys.argv[3]
EvalueLimit = sys.argv[4]
cpus = sys.argv[5]
annotout = sys.argv[6]
error_out = sys.argv[6]

hmmfiles = os.listdir(hmmdir)
genomefiles = os.listdir(genomedir)  #返回指定的文件夹包含的文件或文件夹的名字的列表
genomefiles.sort()
genomefiles = pd.DataFrame(genomefiles)  #哺乳动物基因组很多重复的，直接做465个，不全部做了

genomefiles["species"], genomefiles["detail"] = genomefiles[0].str.split("-",1).str
species_list = genomefiles["species"].value_counts(sort = False)
species_unique = species_list.index[species_list ==1]
species_unique_pd = genomefiles[genomefiles["species"].isin(species_unique)] #只有一个基因组的物种
species_notunique = species_list.index[species_list !=1]
species_multi_pd= genomefiles[genomefiles["species"].isin(species_notunique)] #有多个基因组的物种
species_multi_pd_GCF = species_multi_pd.groupby('species').apply(lambda t: t[t[0].str.contains("GCF")])#优先选用GCF
species_multi_pd_notGCF= species_multi_pd[-species_multi_pd["species"].isin(species_multi_pd_GCF["species"])]
species_multi_pd_notGCF_one = species_multi_pd_notGCF.groupby('species').head(1)#直接取第一个

sort_genomefiles = list(species_unique_pd[0]) + list(species_multi_pd_GCF[0]) + list(species_multi_pd_notGCF_one[0])


'''
#处理得到想要的genome.fasta
def fasta(n):
    if "func_ORs_dna.fasta" in n:
        return True
genomefiles = list(filter( fasta ,genomefiles))
'''

# number of genomic in inputdir Folder
#nums = len(files)
#inpaths = [os.path.join(indir, fi) for fi in files]  #拼接路径名组件，让路径可达
#npaths = [os.path.join(nhmmerout, fi + ".tblout") for fi in files] #给每个文件命名，使用了列表生成式

def batch_annotate():  
    # batch run nhmmer and annotation program
    
    for hmmfile in hmmfiles:
        pofile = os.path.join(hmmdir, hmmfile)
        nhmmout = os.path.join(nhmmerout, genomefile+'_'+ hmmfile + ".tblout")
        annotationout = os.path.join(annotout, genomefile+'_')
        #prefix = infile.split('.')[0]
        # run nhmmer program
        hmmcom = ("python nhmmer.py "
                  + pofile + " "
                  + genome + " "
                  + nhmmout
                  + " -e " + EvalueLimit
                  + " -c " + cpus
                  )
        #if verbose:
            #hmmcom += " -v"
        os.system(hmmcom)
    
        annotationcom = ("python RTP1annotation.py "
                         + nhmmout + " "
                         + genome + " "
                         + EvalueLimit + " "
                         + annotationout + " "
                         + error_out
                         )
        #if verbose:
        #hmmcom += " -v"
        os.system(annotationcom)

if __name__=='__main__':
    
    num_of_jobs = len(sort_genomefiles)
    for i in range(0, num_of_jobs, int(cpus)):
        batch_jobs = sort_genomefiles[i:i+int(cpus)]
        jobs = []
        for genomefile in batch_jobs:
            genome = os.path.join(genomedir, genomefile)
            p = mp.Process(target=batch_annotate)
            p.start()
            jobs.append(p)
        for job in jobs:
            job.join() 
     

    


