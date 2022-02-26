# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 11:13:22 2021

@author: CAVU
"""
import pandas as pd
import numpy as np 
import re

EXON1lengthlimit = 320
EXON2lengthlimit = 610
RTP1seqlen = 750
RTP2seqlen = 650


def proc_nhmmer_out(file, EvalueLimit):#,EXON1lengthlimit,EXON2lengthlimit):
    """
    Function:
        proc_nhmmer_out
        process nhmmer .out file, get sca, fr, to, sign,  & new_fr, new_to
    Parameter:
        file, nhmmer outfile,header like following form
        "target name", "accession", "query name", "accession",
        "hmmfrom", "hmm to", "alifrom", "ali to", "envfrom",
        "env to", "sq len", "strand", "E-value", "score", "bias",
        "description of target"
        EvalueLimit, Sequence similarity threshold
        SeqLengthLimit, Sequence length threshold
    Return:
        return type (hmmout) -> Dictionary
        key   -> gene_name, hit gene name
        value -> a out list
    """

    with open(file) as tbloutf:
        # Filter the header and tail
        linelist = [line for line in tbloutf if line[0:1] != '#']

    hmmout = {}
    for line in linelist:
        temp = line.strip().split()
        # Extract information from NHMMER outfile
        sca = temp[0]
        hmmfr = int(temp[4])
        hmmto = int(temp[5])
        envfr = int(temp[8])
        envto = int(temp[9])
        slen = int(temp[10])
        sign = temp[11]
        evalue = float(temp[12])
        hmmlen = abs(hmmto - hmmfr) + 1
        envlen = abs(envto - envfr) + 1
        gene_name = sca + '_' + str(envfr) + '_' + str(envto) + '_' + sign
        #either exon1 or exon2
        if hmmfr < 20:
            exontype = "exon1"
            #if exon1 length less than exon1lengthlimit, extend it
            extend = EXON1lengthlimit - envlen
            if extend > 0:
                stop_extend = 0
                start_extend = extend 
            else:
                stop_extend = 0
                start_extend = 0
        elif hmmfr > 150:
            exontype = "exon2"
            #if exon1 length less than exon2lengthlimit, extend it
            extend = EXON2lengthlimit - envlen
            if extend > 0:
                stop_extend = extend
                start_extend = 0
            else:
                stop_extend = 0
                start_extend = 0
        # sequence similarity filter
        if evalue > EvalueLimit: continue
        
    
        if sign == '+':
            new_fr = max(envfr - start_extend, 0)
            new_to = min(envto + stop_extend, slen)
        elif sign == '-':
            # The chain of antisense needs to be reversed
            new_fr = max(envto - stop_extend, 0)
            new_to = min(envfr + start_extend, slen)
        else:
            # sign(strand) must be '+' or '-'
            logging.error("NHMMER tool output 'strand' "
                          + "column only '+' or '-'")
            raise StrandError(sign)
        fraglen = abs(new_to - new_fr)+1
        hmmout[gene_name] = [
                sca, hmmfr, hmmto, envfr, envto, new_fr,
                new_to, sign, evalue, slen, hmmlen, envlen, fraglen, exontype
            ]
    return hmmout


#处理配对关系以及序列提取拼接
def pairandmerge(hmmout,genomefile):
    tblout = pd.DataFrame(hmmout)
    tblout = tblout.transpose()
    colnames = [
                "sca", "hmmfr", "hmmto", "envfr", "envto", 
                "new_fr", "new_to", "sign", "evalue", "slen",
                "hmmlen", "envlen", "fraglen", "exontype"
            ]
    tblout.columns = colnames
    #pair relationship
    RTP1 = tblout.groupby('exontype').apply(lambda t: t[t.evalue == t.evalue.min()])
    RTP2 = tblout.groupby('exontype').apply(lambda t: t[t.evalue != t.evalue.min()])
           
    #RTP2可能不会出现 故会重复RTP1
    RTP1_CDS, RTP1_amino , RTP1_error_out= extract_intact(RTP1,RTP1seqlen)
    RTP2_CDS, RTP2_amino , RTP2_error_out= extract_intact(RTP2,RTP2seqlen)
    error_out = "RTP1\n"+RTP1_error_out +"\n"+ "RTP2\n"+RTP2_error_out+"\n"
    
    return RTP1_CDS, RTP1_amino,RTP2_CDS, RTP2_amino, error_out
    
    

def extract_intact(RTP_dict,SeqLengthLimit):
      #RTP_dict.index = RTP_dict.index.droplevel()
    if RTP_dict.shape[0] !=2: 
       error_out = tblout
       all_cds = []
       amino_seq = []
              
  
    elif RTP_dict.shape[0] ==2: 
        error_out = ""
        pair_right(RTP_dict)#打假
        RTP_dict = RTP_dict.T.to_dict('list')
        RTP_seq = extract_cds(RTP_dict,genomefile)
        #exon1_seq extract
        seq = RTP_seq[list(RTP_seq.keys())[0]]
        exon1ends =  find_all('(?=(GCAGG))', seq)
        exon1list = [seq[0:i+5] for i in exon1ends ]
        #exon2_seq extract
        seq = RTP_seq[list(RTP_seq.keys())[1]]
        #exon2start = find_all('(?=(G.AGG))', seq)
        #stops = find_stop_codons(seq, SeqLengthLimit)
        exon2start = find_gcagg(seq)
        exon2list = [ seq[i:] for i in exon2start ]
        #merge
        all_raw_cds= []
        for exon1seq in exon1list:
            for exon2seq in exon2list:
                cds = exon1seq+exon2seq[5:]
                all_raw_cds.append(cds)
    
        all_cdslist = find_cds(all_raw_cds,SeqLengthLimit)
        all_cds, amino_seq = [],[]
        for cds in all_cdslist:
            iatg, istop, cds_seq = cds #location
            cds_seq_retu = cds_seq + "\n"
            all_cds.append(cds_seq_retu)
            amino = dna_translation(cds_seq)
            amino = amino[:-1]
            amino = amino + "\n"
            amino_seq.append(amino)
      

    return all_cds, amino_seq, error_out
    


def find_cds(hmmout_seq, SeqLengthLimit):
    """
    Function:
        find_cds
        try find ATG and STOP codons for each seq,
        put good cds info into fun = {},
    Parameter:
        hmmout, function 'proc_nhmmer_out' output
        hmmout_seq, function 'extract_cds' output
        SeqLengthLimit, artificially set OR's sequence length threshold
    Return:
    return type(funcdict) -> Dictory
    key   -> hit gene name
    value -> assume that OR CDS
    return type(pseudos) -> Dictory
    key   -> hit gene name
    value -> assume that OR pseudogenes
    """
    all_cdslist = []
    for raw_cds in hmmout_seq:
        #raw_cds= all_raw_cds[4]
        len_cds = len(raw_cds)
        starts = find_all('(?=(ATG))', raw_cds)
        #starts = sorted(i for i in starts if i < len_cds - SeqLengthLimit)
        stops = find_stop_codons(raw_cds, SeqLengthLimit)
        cdslist = [(i, j, raw_cds[i:j+3]) for i in starts for j in stops]
        # Determine if all CDS's meet the length limit
        cdslist = cds_length_filter(cdslist, SeqLengthLimit)
        if cdslist:
            # Insert or delete codon (pseudogene)
            cdslist = insert_filter(cdslist)
            if cdslist:
                # Interrupting stop codon (pseudogene)
                cdslist = interrupt_stop_codon(cdslist)
                all_cdslist.extend(cdslist)
                        
    return all_cdslist


def cds_length_filter(cdslist, SeqLengthLimit):
    """
    Function:
        cds_length_filter
        Determine if all CDS's meet the length limit
    Parameter:
        cdslist, cds list from a DNA fragment
    Return:
        return type -> list
        Filtered cdslist
    """
    filterlist = [cds for cds in cdslist if len(cds[2]) >= SeqLengthLimit]
    return filterlist


def insert_filter(cdslist):
    """
    Function:
        insert_filter
        Determine if all CDS's ware inserted or deleted
    Parameter:
        cdslist, cds list from a DNA fragment
    Return:
        return type -> list
        Filtered cdslist
    """
    filterlist = [cds for cds in cdslist if len(cds[2]) % 3 == 0]
    return filterlist

def interrupt_stop_codon(cdslist):
    """
    Function:
        interrupt_stop_codon
        Determine if all CDS ware interrupting stop codon
    Parameter:
        cdslist, cds list from a DNA fragment
    Return:
        return type -> list
        Filtered cdslist
    """
    filterlist = []
    for cds in cdslist:
        codons = re.findall('...', cds[2][:-3])#最终的stopcode不记录在内
        if 'TAG' in codons: continue
        if 'TGA' in codons: continue
        if 'TAA' in codons: continue
        filterlist.append(cds)
    return filterlist


CODON_TABLE = {
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R',
    'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S',
    'AGC': 'S', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'TTA': 'L',
    'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GTT': 'V',
    'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'AAT': 'N', 'AAC': 'N', 'GAT': 'D', 'GAC': 'D',
    'TGT': 'C', 'TGC': 'C', 'CAA': 'Q', 'CAG': 'Q', 'GAA': 'E',
    'GAG': 'E', 'CAT': 'H', 'CAC': 'H', 'AAA': 'K', 'AAG': 'K',
    'TTT': 'F', 'TTC': 'F', 'TAT': 'Y', 'TAC': 'Y', 'ATG': 'M',
    'TGG': 'W', 'TAG': '*', 'TGA': '*', 'TAA': '*'
}



def dna_translation(dna_seq):
    """
    Function:
        dna_translation
        Translate DNA sequence to protein sequence
    Parameter:
        dna_seq, DNA sequence without header
    Return:
        return type -> String
        Translated protein sequence
    """

    protein = ''
    i = 0
    len_dna = len(dna_seq)
    if len_dna % 3 != 0:
        logging.error("The length of CDS({0}) is not divisible by 3."
                      .format(len_dna))
        raise LengthError(len_dna)
    plen = len(dna_seq) / 3
    while i < plen:
        n = dna_seq[3 * i: 3 * (i + 1)]
        if 'N' in n:
            r = '*'
        else:
            r = CODON_TABLE[n]
        i += 1
        protein += r
    return protein




def find_all(substring, string):
    """
    Function:
        find_all
        find all indexes of a substring in a string.
        Overlapping is considered.
    Parameter:
        substring,
        string,
    Return:
        return type -> List
        substring starting position
    """

    indexes = [m.start() for m in re.finditer(substring, string)]
    return indexes


def find_stop_codons(seq, SeqLengthLimit):
    """
    Function:
        find_stop_codons
        Find stop codons in seq, no matter what reading frame is
    Parameter:
        seq, DNA sequence without header
        SeqLengthLimit, artificially set OR's sequence length threshold
    Return:
        return type -> List
        All 'stop codon' of seq(input DNA sequence)
    """

    stop_tag = find_all('(?=(TAG))', seq)
    stop_tga = find_all('(?=(TGA))', seq)
    stop_taa = find_all('(?=(TAA))', seq)
    stops = stop_tag + stop_tga + stop_taa
    # Filter stop long enough
    #stops = sorted(i for i in stops if i > SeqLengthLimit)
    return stops
    
def find_gcagg(seq):
    """
    Function:
        find_stop_codons
        Find stop codons in seq, no matter what reading frame is
    Parameter:
        seq, DNA sequence without header
        SeqLengthLimit, artificially set OR's sequence length threshold
    Return:
        return type -> List
        All 'stop codon' of seq(input DNA sequence)
    """

    stop_g_agg = find_all('(?=(G.AGGTT))', seq)
    stop__cagg = find_all('(?=(.CAGGTT))', seq)
    stop_gcag_ = find_all('(?=(GCAG.TT))', seq)
    GCAGGs = stop_g_agg + stop__cagg + stop_gcag_
    GCAGGs =list(set(GCAGGs))
    # Filter stop long enough
    #GCAGGs = sorted(i for i in GCAGGs if i > SeqLengthLimit)
    return GCAGGs



def pair_right(group):
    if len(list(set(group.sca))) !=1:
        logging.error("sca_error:The loctaion of exon1 and exon2 is not in a scanfold."
                      .format(list(set(group.sca))))
        raise LocationError(list(set(group.sca)))
    if len(list(set(group.sign))) !=1:
        logging.error("sign_error:The loctaion of exon1 and exon2 is in different chromsome."
                      .format(list(set(group.sign))))
        raise LocationError(list(set(group.sign)))
    if group.sign[0] == "+":
        if group.new_fr[group.exontype == "exon2"][0] < group.new_to[group.exontype == "exon1"][0]:
            logging.error("location_error:The loctaion of exon1 and exon2 is wrong."
                      .format(group.newfr[group.exontype == "exon2"]))
            raise LocationError(group.newfr[group.exontype == "exon2"])
    elif group.sign[0] == "-":
        if group.new_fr[group.exontype == "exon2"][0] > group.new_to[group.exontype == "exon1"][0]:
            logging.error("The loctaion of exon1 and exon2 is wrong."
                      .format(group.newfr[group.exontype == "exon2"]))
            raise LocationError(group.newfr[group.exontype == "exon2"])
            


def extract_cds(hmmout, gefile):
    """
    Function:
         extract_cds
         Extract cds from genomic file
    Parameter:
         hmmout, function 'proc_nhmmer_out' output
         gefile, genome sequence file
    Return:
        return type -> Dict
        key   -> hit gene name
        value -> CDS sequence
    """
     
    hmmout_seq, sour_seq = {}, {}
    genomef = open(gefile)
    header, seq_line, flag = '', '', False
    line = genomef.readline()
    if line[0] != '>':
        logging.error("Your genome start with '{0}', FASTA format?"
                      .format(line[0]))
        raise FastaFormatError(line[0])
    while line:
        if line[0] == '>':
            if flag:
                cds  = extract_cds_match(header, hmmout, seq_line)
                hmmout_seq.update(cds)
                #sour_seq.update(sourcds)
                flag = False
            seq_line = ''
            header = line[1:].split()[0]
        else:
            flag = True
            seq_line += line.strip()
        line = genomef.readline()
    # Process last sequence in genome
    cds = extract_cds_match(header, hmmout, seq_line)
    hmmout_seq.update(cds)
    #sour_seq.update(sourcds)

    genomef.close()
    return hmmout_seq#, sour_seq


def extract_cds_match(scaf, hmmout, dna):
    """
    Function:
        extract_cds_match
        Serves ectract_cds function
    Parameter:
        scaf, scaffold name
        hmmout, funtion 'proc_nhmmer_out' output
        dna, scaffold (DNA) squence
    Return:
        return type -> str
        hit, hit gene name from genome
        return type -> str
        cutseq, cut sequence from scaffold
    """

    outdict = {}
    #sourdict = {}
    hmmout = {k: v for k, v in hmmout.items() if v[0] == scaf}
    for hit in hmmout:
        #fr, to = sorted(hmmout[hit][1:3])
        new_fr = hmmout[hit][5]
        new_to = hmmout[hit][6]
        sign = hmmout[hit][7]
        cutseq = dna[new_fr:new_to].upper()
        #sourseq = dna[fr:to].upper()
        seq_replace = ''
        sour_replace = ''
        for s in cutseq:
            # Other non-standard bases are converted to N
            if s in ['A', 'T', 'C', 'G', 'N']:
                seq_replace += s
            else:
                seq_replace += 'N'
        '''for s in sourseq:
            # Other non-standard bases are converted to N
            if s in ['A', 'T', 'C', 'G', 'N']:
                sour_replace += s
            else:
                sour_replace += 'N'
        '''
        cutseq = seq_replace
        #sourseq = sour_replace
        if sign == '-':
            cutseq = reverse_complement(cutseq)
            #sourseq = reverse_complement(sourseq)
        outdict[hit] = cutseq
        #sourdict[hit] = sourseq
    return outdict


def reverse_complement(string):
    """
    Function:
        reverse_complement
        get reverse complement of a DNA string
    Parameter:
        string, a DNA sequence
    Return:
        return type -> String
        Reversed DNA sequence
    """

    str_reverse = string.strip()[::-1]
    try:
        str_comp = ''.join([CompBase[base] for base in str_reverse])
    except KeyError:
        str_replace = ''
        for s in string:
            if s in CompBase.keys():
                str_replace += s
            else:
                str_replace += 'N'
        str_comp = reverse_complement(str_replace)
        logging.error('Illegal letters appear in the genome!')
    return str_comp

# Complementary base
CompBase = {'A': 'T', 'T': 'A', 'C': 'G',
            'G': 'C', 'N': 'N', }


import sys

tblout = sys.argv[1]  
genomefile = sys.argv[2]
EvalueLimit = float(sys.argv[3])
#exon1Limit = int(sys.argv[4])
#exon2Limit = int(sys.argv[5])
annotationout = sys.argv[4]
out_error = sys.argv[5]

hmmout = proc_nhmmer_out(tblout, EvalueLimit)#, exon1Limit, exon2Limit)
rtp1_cds, rtp1_aminoseq,rtp2_cds, rtp2_aminoseq , error_out = pairandmerge(hmmout, genomefile)




'''
file = "/home/cavu/DATA/annotation/RTP/try/tblout/Acinonyx_jubatus-GCA_003709585.1_Aci_jub_2_genomic.fna_RTP1_3.hmm.tblout"
hmmout = proc_nhmmer_out(file, 1e-20)

genomefile = '/home/cavu/DATA/annotation/RTP/try/genome/Acinonyx_jubatus-GCA_003709585.1_Aci_jub_2_genomic.fna'
rtp1_cds, rtp1_aminoseq,rtp2_cds, rtp2_aminoseq, error_out= pairandmerge(hmmout, genomefile)

annotationout = "/home/cavu/DATA/annotation/RTP/try/annotate_res/"
'''
with open(annotationout+"rtp1_cds.txt","w") as f:
    f.writelines(rtp1_cds)

with open(annotationout+"rtp1_aminoseq.txt","w") as f:
    f.writelines(rtp1_aminoseq)
    
with open(annotationout+"rtp2_cds.txt","w") as f:
    f.writelines(rtp2_cds)

with open(annotationout+"rtp2_aminoseq.txt","w") as f:
    f.writelines(rtp2_aminoseq)
    
with open(out_error+"error_out.txt", 'a') as f:
    f.writelines(error_out) 
