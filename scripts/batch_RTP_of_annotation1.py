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

hmmdir = sys.argv[1]  #参数存放hmm的路径
genomedir = sys.argv[2]
nhmmerout = sys.argv[3]
EvalueLimit = sys.argv[4]
cpus = sys.argv[5]
annotout = sys.argv[6]
error_out = sys.argv[6]

hmmfiles = os.listdir(hmmdir)
genomefiles = os.listdir(genomedir)  #返回指定的文件夹包含的文件或文件夹的名字的列表

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

# batch run nhmmer and FindOR program
for genomefile in genomefiles:
    genome = os.path.join(genomedir, genomefile)
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
    


