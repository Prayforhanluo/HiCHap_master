# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 16:18:01 2018

@author: han-luo
"""

import os, subprocess, multiprocessing, logging, time, bisect, cPickle
from HiCHap.mapping import Getchunks
import numpy as np
import pysam

log = logging.getLogger(__name__)



#==========================About genome enzyme fragment=============================

def LoadFragments(FragPath):
    """
    Load fragment data from fragment file created by HapHiC rebuild module.
    
    return : dic
        
        A dict save the fragment Infomation.
        chromosome for keys.
        enzyme sites position for each chromosome.
        
    """
    Frag_type = np.dtype({'names':['chr','start','end'],
                          'formats':['S4', np.int, np.int]})
    
    tmp_Frag = np.loadtxt(FragPath, dtype = Frag_type)
    chroms = set(tmp_Frag['chr'])
    
    Frags = {}
    for c in chroms:
        if c.isdigit() or c == 'X' or c == 'Y':
            tmp_chro = tmp_Frag[tmp_Frag['chr'] == c]
            Frags[c] = np.array([1] + list(tmp_chro['end']))
                
    return Frags


def FragMid(Frags,read):
    """
    Map pos to which fragment. Return the fragment middle point.
    
    """
    chro = read.reference_name.lstrip('chr')
    pos = read.pos + 1
    Frag = Frags[chro]
    index = bisect.bisect_left(Frag,pos)
    
    return (Frag[index - 1] + Frag[index]) // 2
    



#================================About SNP==========================================

def LoadingSNPs(Snps_file):
    """
    Loading The SNPs Information from SNP file created by HapHiC rebuild module.
    
    """
    Snps = cPickle.load(open(Snps_file,'rb'))
    
    return Snps




def SnpsMatch(read,Snps,Allelic):
    
    """
        Allelic SNP matching for a mate.
        
        return the matching SNP site count.
    """
    chro = read.reference_name.lstrip('chr')
    pos = read.pos + 1
    Len = read.query_length
    Seq = read.seq

    start = bisect.bisect_left(Snps[chro]['pos'],pos)
    end = bisect.bisect_left(Snps[chro]['pos'],pos+Len)
    
    Position = Snps[chro]['pos'][start:end]
    if Allelic == 'Maternal':
        alt = Snps[chro]['m_alt'][start:end]
    else:  
        alt = Snps[chro]['p_alt'][start:end]
    read_alt = [Seq[i - pos] for i in Position]
    
    count = 0
    for i in range(len(Position)):
        if read_alt[i] ==alt[i]:
            count += 1
            
    
    return str(count)


#==============================About Mapping Reads==================================
    

def is_unmapped_read(read):
    """
        Whether Unmapped. Scaffold mapping results will be regard as unmapped.
    """
    if read.is_unmapped:
        return 1
    else:
        chro = read.reference_name.lstrip('chr')
        if chro.isdigit() or chro == 'X' or chro == 'Y':
            return 0
        else:
            return 1

        
def is_unique_read(read,level = 1):
    """
     Whether unique mapping.
    """
    mark = 0
    if not is_unmapped_read(read) and read.has_tag('AS'):    
        if level == 1:
            if read.has_tag('XS'):
                mark = 0
            else:
                mark = 1
        else:
            if read.has_tag('XS'):
                first = read.get_tag('AS')
                second = read.get_tag('XS')
                if first > second:
                    mark = 1
            else:
                mark = 1
    
    return mark

#============================About Non-Allelic Pipeline=============================
def Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags):
    """
        Extract pair Information which has no candidate mate.
    """
    name = '_'.join(mate1.query_name.split('_')[:-1])
    
    F1 = FragMid(Frags,mate1)
    F2 = FragMid(Frags,mate2)
    S1 = 0
    S2 = 0
    
    info = [name]+[mate1.reference_name,mate1.flag,mate1.pos+1,
                   mate1.query_length,mate1.get_tag('AS')]+\
           [F1,S1]+[mate2.reference_name,mate2.flag,mate2.pos+1,
                    mate2.query_length,mate2.get_tag('AS')]+\
           [F2,S2]
        
    info = map(str,info)    
    return info


def Info_Generation_There_mate_Non_Allelic(mate1,mate2,candidate,Frags,mate_mark):
    """
        Extract pair Information which has candidate mate.
    """
    name = '_'.join(mate1.query_name.split('_')[:-1])
    
    F1 = FragMid(Frags,mate1)
    F2 = FragMid(Frags,mate2)
    FC = FragMid(Frags,candidate)
    
    S1 = 0
    S2 = 0
    SC = 0
    
    info = [name]+[mate1.reference_name,mate1.flag,mate1.pos+1,
                   mate1.query_length,mate1.get_tag('AS')]+\
           [F1,S1]+[mate2.reference_name,mate2.flag,mate2.pos+1,
                    mate2.query_length,mate2.get_tag('AS')]+\
           [F2,S2]+[candidate.reference_name,candidate.flag,candidate.pos+1,
                    candidate.query_length,candidate.get_tag('AS')]+\
           [FC,SC,mate_mark]
    
    info = map(str,info)
    return info


def Pair_Integrate_Non_Allelic(tmp,Frags,level):
    """
    """
    if len(tmp) == 2:
        #No cutting sites
        for read in tmp:
            if is_unmapped_read(read):
                return 0              # Unmapped
            elif is_unique_read(read,level):
                continue
            else:
                return 1              #Multi-mapped
        
        mate1 = tmp[0]
        mate2 = tmp[1]
        info = Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags)
        
        return info
        
    elif len(tmp) == 3:
        #Cutting a mate but candidate is not available       
        
        un_count = 0
        for read in tmp:
            if is_unmapped_read(read):
                un_count += 1
        if un_count >= 2:
            return 0                   #Unmapped.
        
        multi_count = 0
        for read in tmp:
            if not is_unique_read(read,level):
                multi_count += 1
        if multi_count >= 2:
            return 1                    #Multi-mapped                       
        
        for read in tmp:
            if is_unmapped_read(read):
                continue
            elif read.query_name[-1] == '1':
                mate1 = read
            elif read.query_name[-1] == '2':
                mate2 = read
        
        info = Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags)
        
        return info    

    elif len(tmp) == 4:
        #one mate cutting model
        name_tag = []
        for read in tmp:
            tag = read.query_name.split('_')[-1]
            name_tag.append(tag)
        
        name_tag.sort()
        
        if name_tag == ['1','11','12','2']:
            #R1 mate cutting
            for read in tmp:
                tag = read.query_name.split('_')[-1]
                if tag == '11':
                    mate11 = read
                elif tag == '12':
                    mate12 = read
                elif tag == '2':
                    mate2 = read
                else:
                    continue
            if is_unmapped_read(mate2):
                return 0                  #Unmapped
            elif is_unmapped_read(mate11) and is_unmapped_read(mate12):
                return 0
            elif not is_unique_read(mate2,level):
                return 1                  #Multi-reads
            elif (not is_unique_read(mate11,level)) and (not is_unique_read(mate12,level)):
                return 1
            else:
                if not is_unique_read(mate11,level):
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        return 0
                    else:
                        info = Info_Generation_Two_mate_Non_Allelic(mate12,mate2,Frags)
                        return info
                elif not is_unique_read(mate12,level):
                    info = Info_Generation_Two_mate_Non_Allelic(mate11,mate2,Frags)
                    return info
                else:
                    F11 = FragMid(Frags,mate11)
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        info = Info_Generation_There_mate_Non_Allelic(mate11,mate2,mate12,Frags,'R2')
                        return info
                    elif F11 == F12:
                        info = Info_Generation_There_mate_Non_Allelic(mate11,mate2,mate12,Frags,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate_Non_Allelic(mate11,mate12,Frags)
                        info2 = Info_Generation_Two_mate_Non_Allelic(mate12,mate2,Frags)
                        
                        return Merge_Candidate_interaction(info1,info2)
                        
            
        elif name_tag == ['1','2','21','22']:
            # R2 cutting model
            for read in tmp:
                tag = read.query_name.split('_')[-1]
                if tag == '1':
                    mate1 = read
                elif tag == '21':
                    mate21 = read
                elif tag == '22':
                    mate22 = read
                else:
                    continue
            if is_unmapped_read(mate1):
                return 0
            elif is_unmapped_read(mate21) and is_unmapped_read(mate22):
                return 0
            elif not is_unique_read(mate1,level):
                return 1
            elif (not is_unique_read(mate21,level)) and (not is_unique_read(mate22,level)):
                return 1
            else:
                if not is_unique_read(mate21,level):
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F22 == F1:
                        return 0
                    else:
                        info = Info_Generation_Two_mate_Non_Allelic(mate1,mate22,Frags)
                        return info
                elif not is_unique_read(mate22,level):
                    info = Info_Generation_Two_mate_Non_Allelic(mate1,mate21,Frags)
                    return info
                else:
                    F21 = FragMid(Frags,mate21)
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F21 == F22:
                        info = Info_Generation_There_mate_Non_Allelic(mate1,mate21,mate22,Frags,'R2')
                        return info
                    elif F22 == F1:
                        info = Info_Generation_There_mate_Non_Allelic(mate1,mate21,mate22,Frags,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate_Non_Allelic(mate1,mate22,Frags)
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_Two_mate_Non_Allelic(mate22,mate21,Frags)
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
        elif name_tag == ['1','1','2','2']:
            #two mate cutting but all candidates are not available
            new_tmp = []
            for read in tmp:
                if read.query_length == 150:
                    continue
                else:
                    new_tmp.append(read)
            
            for read in new_tmp:
                if is_unmapped_read(read):
                    return 0              # Unmapped
                elif is_unique_read(read,level):
                    continue
                else:
                    return 1              #Multi-mapped
        
            mate1 = new_tmp[0]
            mate2 = new_tmp[1]
            info = Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags)
            
            return info
        else:
            info = ''
            return info

    elif len(tmp) == 5:
        #two mate cutting but one candidate is not available
        name_tag = []
        for read in tmp:
            tag = read.query_name.split('_')[-1]
            name_tag.append(tag)
        
        name_tag.sort()
        
        if name_tag == ['1','11','12','2','2']:
            
            for read in tmp:
                tag = read.query_name.split('_')[-1]
                
                if tag == '2' and read.query_length < 150:
                    mate2 = read
                elif tag == '11':
                    mate11 = read
                elif tag == '12':
                    mate12 = read
                else:
                    continue
            
            if is_unmapped_read(mate2):
                return 0                # Unmapped
            elif is_unmapped_read(mate11) and is_unmapped_read(mate12):
                return 0
            elif not is_unique_read(mate2,level):
                return 1                # Multi-mapped   
            elif (not is_unique_read(mate11,level)) and (not is_unique_read(mate12,level)):
                return 1
            else:
                if not is_unique_read(mate11,level):
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        return 0
                    else:
                        info = Info_Generation_Two_mate_Non_Allelic(mate12,mate2,Frags)
                        return info
                elif not is_unique_read(mate12,level):
                    info = Info_Generation_Two_mate_Non_Allelic(mate11,mate2,Frags)
                    return info
                else:
                    F11 = FragMid(Frags,mate11)
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        info = Info_Generation_There_mate_Non_Allelic(mate11,mate2,mate12,Frags,'R2')
                        return info
                    elif F11 == F12:
                        info = Info_Generation_There_mate_Non_Allelic(mate11,mate2,mate12,Frags,'R1')
                        return info
                        
                    else:
                        info1 = Info_Generation_Two_mate_Non_Allelic(mate11,mate12,Frags)
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_Two_mate_Non_Allelic(mate12,mate2,Frags)
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
        
        elif name_tag == ['1','1','2','21','22']:
            
            for read in tmp:
                tag = read.query_name.split('_')[-1]
                
                if tag == '1' and read.query_length < 150:
                    mate1 = read
                elif tag == '21':
                    mate21 = read
                elif tag == '22':
                    mate22 = read
                else:
                    continue
            
            if is_unmapped_read(mate1):
                return 0
            elif is_unmapped_read(mate21) and is_unmapped_read(mate22):
                return 0
            elif not is_unique_read(mate1,level):
                return 1
            elif (not is_unique_read(mate21,level)) and (not is_unique_read(mate22,level)):
                return 1
            else:
                if not is_unique_read(mate21,level):
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F22 == F1:
                        return 0
                    else:
                        info = Info_Generation_Two_mate_Non_Allelic(mate1,mate22,Frags)
                        return info
                elif not is_unique_read(mate22,level):
                    info = Info_Generation_Two_mate_Non_Allelic(mate1,mate21,Frags)
                    return info
                else:
                    F21 = FragMid(Frags,mate21)
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F21 == F22:
                        info = Info_Generation_There_mate_Non_Allelic(mate1,mate21,mate22,Frags,'R2')
                        return info
                    elif F22 == F1:
                        info = Info_Generation_There_mate_Non_Allelic(mate1,mate21,mate22,Frags,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate_Non_Allelic(mate1,mate22,Frags)
                        #info1[0] = info1[0] + '_1'
                        
                        info2 = Info_Generation_Two_mate_Non_Allelic(mate22,mate21,Frags)
                        #info2[0] = info2[0] + '_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
        else:
            info = ''
            return info

    elif len(tmp) == 6:
        #two mate cutting model
        for read in tmp:
            tag = read.query_name.split('_')[-1]
            if tag == '11':
                mate11 = read
            elif tag == '12':
                mate12 = read
            elif tag == '21':
                mate21 = read
            elif tag == '22':
                mate22 = read
            else:
                continue
        
        if is_unmapped_read(mate11) and is_unmapped_read(mate12):
            return 0
        elif is_unmapped_read(mate21) and is_unmapped_read(mate22):
            return 0
        elif (not is_unique_read(mate11,level)) and (not is_unique_read(mate12,level)):
            return 1
        elif (not is_unique_read(mate21,level)) and (not is_unique_read(mate22,level)):
            return 1
        else:
            #pairs can be rescued.
            if not is_unique_read(mate11,level):
                mate1 = mate12
                if not is_unique_read(mate22,level):
                    mate2 = mate21
                    info = Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags)
                    return info
                elif not is_unique_read(mate21,level):
                    mate2 = mate22
                    info = Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags)
                    return info
                else:
                    F21 = FragMid(Frags,mate21)
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F21 == F22:
                        info = Info_Generation_There_mate_Non_Allelic(mate1,mate21,mate22,Frags,'R2')
                        return info
                    elif F22 == F1:
                        info = Info_Generation_There_mate_Non_Allelic(mate1,mate21,mate22,Frags,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate_Non_Allelic(mate1,mate22,Frags)
                        #info1[0] = info1[0] + '_1'
                        
                        info2 = Info_Generation_Two_mate_Non_Allelic(mate22,mate21,Frags)
                        #info2[0] = info2[0] + '_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
            elif not is_unique_read(mate12,level):
                mate1 = mate11
                if not is_unique_read(mate22,level):
                    mate2 = mate21
                    info = Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags)
                    return info
                elif not is_unique_read(mate21,level):
                    mate2 = mate22
                    info = Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags)
                    return info
                else:
                    F21 = FragMid(Frags,mate21)
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F21 == F22:
                        info = Info_Generation_There_mate_Non_Allelic(mate1,mate21,mate22,Frags,'R2')
                        return info
                    elif F22 == F1:
                        info = Info_Generation_There_mate_Non_Allelic(mate1,mate21,mate22,Frags,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate_Non_Allelic(mate1,mate22,Frags)
                        #info1[0] = info1[0] + '_1'
                        
                        info2 = Info_Generation_Two_mate_Non_Allelic(mate22,mate21,Frags)
                        #info2[0] = info2[0] + '_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
            elif not is_unique_read(mate22,level):
                mate2 = mate21
                if not is_unique_read(mate11,level):
                    mate1 = mate12
                    info = Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags)
                    return info
                elif not is_unique_read(mate12,level):
                    mate1 = mate11
                    info = Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags)
                    return info
                else:
                    F11 = FragMid(Frags,mate11)
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        info = Info_Generation_There_mate_Non_Allelic(mate11,mate2,mate12,Frags,'R2')
                        return info
                    elif F11 == F12:
                        info = Info_Generation_There_mate_Non_Allelic(mate11,mate2,mate12,Frags,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate_Non_Allelic(mate11,mate12,Frags)
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_Two_mate_Non_Allelic(mate12,mate2,Frags)
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
            elif not is_unique_read(mate21,level):
                mate2 = mate22
                if not is_unique_read(mate11,level):
                    mate1 = mate12
                    info = Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags)
                    return info
                elif not is_unique_read(mate12,level):
                    mate1 = mate11
                    info = Info_Generation_Two_mate_Non_Allelic(mate1,mate2,Frags)
                    return info
                else:
                    F11 = FragMid(Frags,mate11)
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        info = Info_Generation_There_mate_Non_Allelic(mate11,mate2,mate12,Frags,'R2')
                        return info
                    elif F11 == F12:
                        info = Info_Generation_There_mate_Non_Allelic(mate11,mate2,mate12,Frags,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate_Non_Allelic(mate11,mate12,Frags)
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_Two_mate_Non_Allelic(mate12,mate2,Frags)
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1, info2)
            else:
                F11 = FragMid(Frags,mate11)
                F12 = FragMid(Frags,mate12)
                F22 = FragMid(Frags,mate22)
                F21 = FragMid(Frags,mate21)
                if F11 == F12:
                    if F22 == F21:
                        info1 = Info_Generation_There_mate_Non_Allelic(mate11,mate21,mate22,Frags,'R2')
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_There_mate_Non_Allelic(mate12,mate21,mate22,Frags,'R2')
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1, info2)
                    else:
                        info1 = Info_Generation_There_mate_Non_Allelic(mate11,mate22,mate12,Frags,'R1')
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_There_mate_Non_Allelic(mate12,mate21,mate12,Frags,'R1')
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1, info2)
                else:
                    if F22 == F21:
                        info1 = Info_Generation_There_mate_Non_Allelic(mate11,mate21,mate22,Frags,'R2')
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_There_mate_Non_Allelic(mate12,mate21,mate22,Frags,'R2')
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1, info2)
                    else:
                        if F12 == F22:
                            info1 = Info_Generation_There_mate_Non_Allelic(mate11,mate22,mate12,Frags,'R2')
                            #info1[0] = info1[0]+'_1'
                            
                            info2 = Info_Generation_There_mate_Non_Allelic(mate12,mate21,mate22,Frags,'R1')
                            #info2[0] = info2[0]+'_2'
                            
                            return Merge_Candidate_interaction(info1, info2)
                        else:
                            info1 = Info_Generation_Two_mate_Non_Allelic(mate11,mate12,Frags)
                            #info1[0] = info1[0]+'_1'
                            
                            info2 = Info_Generation_Two_mate_Non_Allelic(mate22,mate21,Frags)
                            #info2[0] = info2[0]+'_2'
                            
                            return Merge_Candidate_interaction(info1, info2)
                            
    else:
        logging.error( "error tag number:%s" %len(tmp))
#        for i in tmp:
#            logging.error('error tag name %s' % i.query_name)
        return ''




def Bam_Intergrate_Non_Allelic(bam1,Rebam1,bam2,Rebam2,OutPath,Frags,level):
    """
        Extract Information from bam file.
        We will generate the merged bam of two mate first, And then we will matching
        the mapping results from the same pairs.This way, non-unique reads will be 
        filtered out.On the other hand, We will search the read seq to count the Allelic
        SNP sites on this reads.Finally We will keep those Information for a line:
        
    column    :    Means
         0    :    Pair Name
         1    :    R1 mate Reference
         2    :    R1 mate Strand
         3    :    R1 mate Position
         4    :    R1 mate Length
         5    :    R1 mate AS score
         6    :    R1 mate Fragment Middle point
         7    :    R1 mate SNP Matching num
         
         8    :    R2 mate Reference
         9    :    R2 mate Strand
        10    :    R2 mate Position
        11    :    R2 mate Length
        12    :    R2 mate AS score
        13    :    R2 mate Fragment Middle point
        14    :    R2 mate SNP Matching num
        
        (candidate mate if it is possible)
        
        15    :    Candidate mate Reference
        16    :    Candidate mate Strand
        17    :    Candidate mate Position
        18    :    Candidate mate Length
        19    :    Candidate mate AS score
        20    :    Candidate mate Fragment Middle point
        21    :    Candidate mate SNP Matching num
        22    :    Candidate Index for which mate.
        
    """
    tmp1 = os.path.split(bam1)[1].split('_')
    tmp2 = os.path.split(bam2)[1].split('_')
    out_bam = '_'.join([tmp1[i] for i in range(len(tmp1)) if tmp1[i] == tmp2[i]])
    out_bam = os.path.join(OutPath,out_bam)
    subprocess.call(' '.join(['samtools','merge','-n','-f',out_bam,bam1,bam2,Rebam1,Rebam2]),
                    shell = True)
    
    out_bed = out_bam.replace('.bam','.bed')
    out = open(out_bed,'w')
    bam_f = pysam.AlignmentFile(out_bam,'rb')
    
    Unmapped_num = 0
    Multi_num = 0
    Total_num = 0
    tag = ''
    tmp = []
    
    for read in bam_f:
        name = read.query_name
        name = '_'.join(name.split('_')[:-1])
        if name != tag and len(tmp) == 0:
            tag = name
            tmp.append(read)
        elif name != tag and len(tmp) != 0:
            Total_num += 1
            tag = name
            info = Pair_Integrate_Non_Allelic(tmp,Frags,level)
            if info == 0:
                Unmapped_num += 1
            elif info == 1:
                Multi_num += 1
            elif info == '':
                Unmapped_num += 1
            else:
                if type(info) == tuple:
                    for i in info:
                        out.writelines('\t'.join(i)+'\n')
                else:
                    out.writelines('\t'.join(info)+'\n')
            tmp = []
            tmp.append(read)
        else:
            tmp.append(read)
    
    info = Pair_Integrate_Non_Allelic(tmp,Frags,level)
    Total_num += 1
    if info == 0:
        Unmapped_num += 1
    elif info == 1:
        Multi_num += 1
    elif info == '':
        Unmapped_num += 1
    else:
        if type(info) == tuple:
            for i in info:
                out.writelines('\t'.join(i)+'\n')
        else:
            out.writelines('\t'.join(info)+'\n')
    out.close()
    
    subprocess.call(['rm', out_bam])
    log.log(21,'%s is completed',out_bed)
    return Total_num, Unmapped_num, Multi_num



def Bam_Extract_Non_Allelic(Bam_Path,Re_Bam_Path,Out_Path,Frag,num,level):
    """
        Merge the two mates results of a pair mapping results Extract the information
        We need and transform To the bed format.
        
        Parameters
        ----------
        Bam_Path : str
            bam Path generated by GlobalMapping module.
        
        Re_Bam_Path : str
            ReMapping Path generated by ReMapping module.
        
        Out_Path : str
            Output Path of bed files.
            
        Frag : str
            Genome Fragment path generate by rebuild genome.
        
        num : int
            threads number.
        
        level : int
            level of filtering Unique pairs.
            
    """
    log.log(21,'Loading Fragments...')
    Frags = LoadFragments(Frag)
    log.log(21,'Loading Done!')
    log.log(21,'Bam information Extraction , Fragments Mapping for reads ...')
    chunks,chunks_num,Cell = Getchunks(Bam_Path)
    Rechunks,Re_chunks_num,Cell = Getchunks(Re_Bam_Path)
    pool = multiprocessing.Pool(num)
    Statics = []
    
    for i in range(chunks_num):
        bam1name = [fil for fil in chunks if 'chunk'+str(i)+'_1' in fil][0]
        bam2name = [fil for fil in chunks if 'chunk'+str(i)+'_2' in fil][0]
        Rebam1name = [fil for fil in Rechunks if 'chunk'+str(i)+'_1' in fil][0]
        Rebam2name = [fil for fil in Rechunks if 'chunk'+str(i)+'_2' in fil][0]
        bam1 = os.path.join(Bam_Path,bam1name)
        bam2 = os.path.join(Bam_Path,bam2name)
        Rebam1 = os.path.join(Re_Bam_Path,Rebam1name)
        Rebam2 = os.path.join(Re_Bam_Path,Rebam2name)
        S = pool.apply_async(Bam_Intergrate_Non_Allelic,args=(bam1,Rebam1,bam2,Rebam2,
                                                              Out_Path,Frags,level))
        Statics.append(S)
        
    time.sleep(2)
    pool.close()
    pool.join()
    
    #Statics
    Total_pairs = 0
    Unmapped_pairs = 0
    Multiple_pairs = 0
    
    for res in Statics:
        tmp = res.get()
        Total_pairs += tmp[0]
        Unmapped_pairs += tmp[1]
        Multiple_pairs += tmp[2]
    
    Statics_txt = ['Non-Allelic Mapping Statics : ',
                   '==============================',
                   'Total    pairs : %d' % Total_pairs,
                   'Unmapped pairs : %d' % Unmapped_pairs,
                   'Multiple pairs : %d' % Multiple_pairs,
                   'Unique   pairs : %d' % (Total_pairs - Unmapped_pairs - Multiple_pairs)]
    logging.log(21,'\n'+'\n'.join(Statics_txt)+'\n')
        


    
#==============================About Allelic Pipeline===============================
def Merge_Candidate_interaction(info1,info2):
    """
        When Candidate Interaction happened.
        Merge the interaction that looks like the same Interaction.
    """
    chr11 = info1[1];chr12 = info1[8]
    chr21 = info2[1];chr22 = info2[8]
    if chr11 == chr21 and chr12 == chr22:
        F11 = int(info1[6]);F12 = int(info1[13])
        F21 = int(info2[6]);F22 = int(info2[13])
        if F11 == F21 and F12 == F22:
            return info1
        else:
            info1[0] = info1[0]+'_1'
            info2[0] = info2[0]+'_2'
            return info1,info2
    else:
        info1[0] = info1[0]+'_1'
        info2[0] = info2[0]+'_2'
        return info1,info2
        
    
def Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic):
    """
        Extract pair Information which has no candidate mate.
    """
    name = '_'.join(mate1.query_name.split('_')[:-1])
    
    F1 = FragMid(Frags,mate1)
    F2 = FragMid(Frags,mate2)
    
    S1 = SnpsMatch(mate1,Snps,Allelic)
    S2 = SnpsMatch(mate2,Snps,Allelic)
    
    info = [name]+[mate1.reference_name,mate1.flag,mate1.pos+1,
                   mate1.query_length,mate1.get_tag('AS')]+\
           [F1,S1]+[mate2.reference_name,mate2.flag,mate2.pos+1,
                    mate2.query_length,mate2.get_tag('AS')]+\
           [F2,S2]
        
    info = map(str,info)    
    return info

def Info_Generation_There_mate(mate1,mate2,candidate,Frags,Snps,Allelic,mate_mark):
    """
        Extract pair Information which has candidate mate.
    """
    name = '_'.join(mate1.query_name.split('_')[:-1])
    
    F1 = FragMid(Frags,mate1)
    F2 = FragMid(Frags,mate2)
    FC = FragMid(Frags,candidate)
    
    S1 = SnpsMatch(mate1,Snps,Allelic)
    S2 = SnpsMatch(mate2,Snps,Allelic)
    SC = SnpsMatch(candidate,Snps,Allelic)
    
    info = [name]+[mate1.reference_name,mate1.flag,mate1.pos+1,
                   mate1.query_length,mate1.get_tag('AS')]+\
           [F1,S1]+[mate2.reference_name,mate2.flag,mate2.pos+1,
                    mate2.query_length,mate2.get_tag('AS')]+\
           [F2,S2]+[candidate.reference_name,candidate.flag,candidate.pos+1,
                    candidate.query_length,candidate.get_tag('AS')]+\
           [FC,SC,mate_mark]
    
    info = map(str,info)
    return info
    


def Pair_Integrate(tmp,Frags,Snps,Allelic,level):
    """
    """
    
    if len(tmp) == 2:
        #No cutting sites
        for read in tmp:
            if is_unmapped_read(read):
                return 0              # Unmapped
            elif is_unique_read(read,level):
                continue
            else:
                return 1              #Multi-mapped
        
        mate1 = tmp[0]
        mate2 = tmp[1]
        info = Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic)
        
        return info

    elif len(tmp) == 3:
        #Cutting a mate but candidate is not available       
        
        un_count = 0
        for read in tmp:
            if is_unmapped_read(read):
                un_count += 1
        if un_count >= 2:
            return 0                   #Unmapped.
        
        multi_count = 0
        for read in tmp:
            if not is_unique_read(read,level):
                multi_count += 1
        if multi_count >= 2:
            return 1                    #Multi-mapped                       
        
        for read in tmp:
            if is_unmapped_read(read):
                continue
            elif read.query_name[-1] == '1':
                mate1 = read
            elif read.query_name[-1] == '2':
                mate2 = read
        
        info = Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic)
        
        return info
        
    elif len(tmp) == 4:
        #one mate cutting model
        name_tag = []
        for read in tmp:
            tag = read.query_name.split('_')[-1]
            name_tag.append(tag)
        
        name_tag.sort()
        
        if name_tag == ['1','11','12','2']:
            #R1 mate cutting
            for read in tmp:
                tag = read.query_name.split('_')[-1]
                if tag == '11':
                    mate11 = read
                elif tag == '12':
                    mate12 = read
                elif tag == '2':
                    mate2 = read
                else:
                    continue
            if is_unmapped_read(mate2):
                return 0                  #Unmapped
            elif is_unmapped_read(mate11) and is_unmapped_read(mate12):
                return 0
            elif not is_unique_read(mate2,level):
                return 1                  #Multi-reads
            elif (not is_unique_read(mate11,level)) and (not is_unique_read(mate12,level)):
                return 1
            else:
                if not is_unique_read(mate11,level):
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        return 0
                    else:
                        info = Info_Generation_Two_mate(mate12,mate2,Frags,Snps,Allelic)
                        return info
                elif not is_unique_read(mate12,level):
                    info = Info_Generation_Two_mate(mate11,mate2,Frags,Snps,Allelic)
                    return info
                else:
                    F11 = FragMid(Frags,mate11)
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        info = Info_Generation_There_mate(mate11,mate2,mate12,Frags,
                                                          Snps,Allelic,'R2')
                        return info
                    elif F11 == F12:
                        info = Info_Generation_There_mate(mate11,mate2,mate12,Frags,
                                                          Snps,Allelic,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate(mate11,mate12,Frags,Snps,Allelic)
                        info2 = Info_Generation_Two_mate(mate12,mate2,Frags,Snps,Allelic)
                        
                        return Merge_Candidate_interaction(info1,info2)
                        
            
        elif name_tag == ['1','2','21','22']:
            # R2 cutting model
            for read in tmp:
                tag = read.query_name.split('_')[-1]
                if tag == '1':
                    mate1 = read
                elif tag == '21':
                    mate21 = read
                elif tag == '22':
                    mate22 = read
                else:
                    continue
            if is_unmapped_read(mate1):
                return 0
            elif is_unmapped_read(mate21) and is_unmapped_read(mate22):
                return 0
            elif not is_unique_read(mate1,level):
                return 1
            elif (not is_unique_read(mate21,level)) and (not is_unique_read(mate22,level)):
                return 1
            else:
                if not is_unique_read(mate21,level):
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F22 == F1:
                        return 0
                    else:
                        info = Info_Generation_Two_mate(mate1,mate22,Frags,Snps,Allelic)
                        return info
                elif not is_unique_read(mate22,level):
                    info = Info_Generation_Two_mate(mate1,mate21,Frags,Snps,Allelic)
                    return info
                else:
                    F21 = FragMid(Frags,mate21)
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F21 == F22:
                        info = Info_Generation_There_mate(mate1,mate21,mate22,Frags,
                                                          Snps,Allelic,'R2')
                        return info
                    elif F22 == F1:
                        info = Info_Generation_There_mate(mate1,mate21,mate22,Frags,
                                                          Snps,Allelic,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate(mate1,mate22,Frags,Snps,Allelic)
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_Two_mate(mate22,mate21,Frags,Snps,Allelic)
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
        elif name_tag == ['1','1','2','2']:
            #two mate cutting but all candidates are not available
            new_tmp = []
            for read in tmp:
                if read.query_length == 150:
                    continue
                else:
                    new_tmp.append(read)
            
            for read in new_tmp:
                if is_unmapped_read(read):
                    return 0              # Unmapped
                elif is_unique_read(read,level):
                    continue
                else:
                    return 1              #Multi-mapped
        
            mate1 = new_tmp[0]
            mate2 = new_tmp[1]
            info = Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic)
            
            return info
        else:
            info = ''
            return info
            
    elif len(tmp) == 5:
        #two mate cutting but one candidate is not available
        name_tag = []
        for read in tmp:
            tag = read.query_name.split('_')[-1]
            name_tag.append(tag)
        
        name_tag.sort()
        
        if name_tag == ['1','11','12','2','2']:
            
            for read in tmp:
                tag = read.query_name.split('_')[-1]
                
                if tag == '2' and read.query_length < 150:
                    mate2 = read
                elif tag == '11':
                    mate11 = read
                elif tag == '12':
                    mate12 = read
                else:
                    continue
            
            if is_unmapped_read(mate2):
                return 0                # Unmapped
            elif is_unmapped_read(mate11) and is_unmapped_read(mate12):
                return 0
            elif not is_unique_read(mate2,level):
                return 1                # Multi-mapped   
            elif (not is_unique_read(mate11,level)) and (not is_unique_read(mate12,level)):
                return 1
            else:
                if not is_unique_read(mate11,level):
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        return 0
                    else:
                        info = Info_Generation_Two_mate(mate12,mate2,Frags,Snps,Allelic)
                        return info
                elif not is_unique_read(mate12,level):
                    info = Info_Generation_Two_mate(mate11,mate2,Frags,Snps,Allelic)
                    return info
                else:
                    F11 = FragMid(Frags,mate11)
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        info = Info_Generation_There_mate(mate11,mate2,mate12,Frags,
                                                          Snps,Allelic,'R2')
                        return info
                    elif F11 == F12:
                        info = Info_Generation_There_mate(mate11,mate2,mate12,Frags,
                                                          Snps,Allelic,'R1')
                        return info
                        
                    else:
                        info1 = Info_Generation_Two_mate(mate11,mate12,Frags,Snps,Allelic)
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_Two_mate(mate12,mate2,Frags,Snps,Allelic)
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
        
        elif name_tag == ['1','1','2','21','22']:
            
            for read in tmp:
                tag = read.query_name.split('_')[-1]
                
                if tag == '1' and read.query_length < 150:
                    mate1 = read
                elif tag == '21':
                    mate21 = read
                elif tag == '22':
                    mate22 = read
                else:
                    continue
            
            if is_unmapped_read(mate1):
                return 0
            elif is_unmapped_read(mate21) and is_unmapped_read(mate22):
                return 0
            elif not is_unique_read(mate1,level):
                return 1
            elif (not is_unique_read(mate21,level)) and (not is_unique_read(mate22,level)):
                return 1
            else:
                if not is_unique_read(mate21,level):
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F22 == F1:
                        return 0
                    else:
                        info = Info_Generation_Two_mate(mate1,mate22,Frags,Snps,Allelic)
                        return info
                elif not is_unique_read(mate22,level):
                    info = Info_Generation_Two_mate(mate1,mate21,Frags,Snps,Allelic)
                    return info
                else:
                    F21 = FragMid(Frags,mate21)
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F21 == F22:
                        info = Info_Generation_There_mate(mate1,mate21,mate22,Frags,
                                                          Snps,Allelic,'R2')
                        return info
                    elif F22 == F1:
                        info = Info_Generation_There_mate(mate1,mate21,mate22,Frags,
                                                          Snps,Allelic,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate(mate1,mate22,Frags,Snps,Allelic)
                        #info1[0] = info1[0] + '_1'
                        
                        info2 = Info_Generation_Two_mate(mate22,mate21,Frags,Snps,Allelic)
                        #info2[0] = info2[0] + '_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
        else:
            info = ''
            return info
        
    elif len(tmp) == 6:
        #two mate cutting model
        for read in tmp:
            tag = read.query_name.split('_')[-1]
            if tag == '11':
                mate11 = read
            elif tag == '12':
                mate12 = read
            elif tag == '21':
                mate21 = read
            elif tag == '22':
                mate22 = read
            else:
                continue
        
        if is_unmapped_read(mate11) and is_unmapped_read(mate12):
            return 0
        elif is_unmapped_read(mate21) and is_unmapped_read(mate22):
            return 0
        elif (not is_unique_read(mate11,level)) and (not is_unique_read(mate12,level)):
            return 1
        elif (not is_unique_read(mate21,level)) and (not is_unique_read(mate22,level)):
            return 1
        else:
            #pairs can be rescued.
            if not is_unique_read(mate11,level):
                mate1 = mate12
                if not is_unique_read(mate22,level):
                    mate2 = mate21
                    info = Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic)
                    return info
                elif not is_unique_read(mate21,level):
                    mate2 = mate22
                    info = Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic)
                    return info
                else:
                    F21 = FragMid(Frags,mate21)
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F21 == F22:
                        info = Info_Generation_There_mate(mate1,mate21,mate22,Frags,
                                                          Snps,Allelic,'R2')
                        return info
                    elif F22 == F1:
                        info = Info_Generation_There_mate(mate1,mate21,mate22,Frags,
                                                          Snps,Allelic,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate(mate1,mate22,Frags,Snps,Allelic)
                        #info1[0] = info1[0] + '_1'
                        
                        info2 = Info_Generation_Two_mate(mate22,mate21,Frags,Snps,Allelic)
                        #info2[0] = info2[0] + '_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
            elif not is_unique_read(mate12,level):
                mate1 = mate11
                if not is_unique_read(mate22,level):
                    mate2 = mate21
                    info = Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic)
                    return info
                elif not is_unique_read(mate21,level):
                    mate2 = mate22
                    info = Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic)
                    return info
                else:
                    F21 = FragMid(Frags,mate21)
                    F22 = FragMid(Frags,mate22)
                    F1 = FragMid(Frags,mate1)
                    if F21 == F22:
                        info = Info_Generation_There_mate(mate1,mate21,mate22,Frags,
                                                          Snps,Allelic,'R2')
                        return info
                    elif F22 == F1:
                        info = Info_Generation_There_mate(mate1,mate21,mate22,Frags,
                                                          Snps,Allelic,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate(mate1,mate22,Frags,Snps,Allelic)
                        #info1[0] = info1[0] + '_1'
                        
                        info2 = Info_Generation_Two_mate(mate22,mate21,Frags,Snps,Allelic)
                        #info2[0] = info2[0] + '_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
            elif not is_unique_read(mate22,level):
                mate2 = mate21
                if not is_unique_read(mate11,level):
                    mate1 = mate12
                    info = Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic)
                    return info
                elif not is_unique_read(mate12,level):
                    mate1 = mate11
                    info = Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic)
                    return info
                else:
                    F11 = FragMid(Frags,mate11)
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        info = Info_Generation_There_mate(mate11,mate2,mate12,Frags,
                                                          Snps,Allelic,'R2')
                        return info
                    elif F11 == F12:
                        info = Info_Generation_There_mate(mate11,mate2,mate12,Frags,
                                                          Snps,Allelic,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate(mate11,mate12,Frags,Snps,Allelic)
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_Two_mate(mate12,mate2,Frags,Snps,Allelic)
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1,info2)
            elif not is_unique_read(mate21,level):
                mate2 = mate22
                if not is_unique_read(mate11,level):
                    mate1 = mate12
                    info = Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic)
                    return info
                elif not is_unique_read(mate12,level):
                    mate1 = mate11
                    info = Info_Generation_Two_mate(mate1,mate2,Frags,Snps,Allelic)
                    return info
                else:
                    F11 = FragMid(Frags,mate11)
                    F12 = FragMid(Frags,mate12)
                    F2 = FragMid(Frags,mate2)
                    if F12 == F2:
                        info = Info_Generation_There_mate(mate11,mate2,mate12,Frags,
                                                          Snps,Allelic,'R2')
                        return info
                    elif F11 == F12:
                        info = Info_Generation_There_mate(mate11,mate2,mate12,Frags,
                                                          Snps,Allelic,'R1')
                        return info
                    else:
                        info1 = Info_Generation_Two_mate(mate11,mate12,Frags,Snps,Allelic)
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_Two_mate(mate12,mate2,Frags,Snps,Allelic)
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1, info2)
            else:
                F11 = FragMid(Frags,mate11)
                F12 = FragMid(Frags,mate12)
                F22 = FragMid(Frags,mate22)
                F21 = FragMid(Frags,mate21)
                if F11 == F12:
                    if F22 == F21:
                        info1 = Info_Generation_There_mate(mate11,mate21,mate22,Frags,
                                                           Snps,Allelic,'R2')
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_There_mate(mate12,mate21,mate22,Frags,
                                                           Snps,Allelic,'R2')
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1, info2)
                    else:
                        info1 = Info_Generation_There_mate(mate11,mate22,mate12,Frags,
                                                           Snps,Allelic,'R1')
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_There_mate(mate12,mate21,mate12,Frags,
                                                           Snps,Allelic,'R1')
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1, info2)
                else:
                    if F22 == F21:
                        info1 = Info_Generation_There_mate(mate11,mate21,mate22,Frags,
                                                           Snps,Allelic,'R2')
                        #info1[0] = info1[0]+'_1'
                        
                        info2 = Info_Generation_There_mate(mate12,mate21,mate22,Frags,
                                                           Snps,Allelic,'R2')
                        #info2[0] = info2[0]+'_2'
                        
                        return Merge_Candidate_interaction(info1, info2)
                    else:
                        if F12 == F22:
                            info1 = Info_Generation_There_mate(mate11,mate22,mate12,Frags,
                                                               Snps,Allelic,'R2')
                            #info1[0] = info1[0]+'_1'
                            
                            info2 = Info_Generation_There_mate(mate12,mate21,mate22,Frags,
                                                               Snps,Allelic,'R1')
                            #info2[0] = info2[0]+'_2'
                            
                            return Merge_Candidate_interaction(info1, info2)
                        else:
                            info1 = Info_Generation_Two_mate(mate11,mate12,Frags,Snps,Allelic)
                            #info1[0] = info1[0]+'_1'
                            
                            info2 = Info_Generation_Two_mate(mate22,mate21,Frags,Snps,Allelic)
                            #info2[0] = info2[0]+'_2'
                            
                            return Merge_Candidate_interaction(info1, info2)
                            
    else:
        logging.error( "error tag number:%s" %len(tmp))
        for i in tmp:
            logging.error('error tag name %s' % i.query_name)
        return ''


def Bam_Intergrate(bam1,Rebam1,bam2,Rebam2,OutPath,Frags,Snps,Allelic,level):
    """
        Extract Information from bam file.
        We will generate the merged bam of two mate first, And then we will matching
        the mapping results from the same pairs.This way, non-unique reads will be 
        filtered out.On the other hand, We will search the read seq to count the Allelic
        SNP sites on this reads.Finally We will keep those Information for a line:
        
    column    :    Means
         0    :    Pair Name
         1    :    R1 mate Reference
         2    :    R1 mate Strand
         3    :    R1 mate Position
         4    :    R1 mate Length
         5    :    R1 mate AS score
         6    :    R1 mate Fragment Middle point
         7    :    R1 mate SNP Matching num
         
         8    :    R2 mate Reference
         9    :    R2 mate Strand
        10    :    R2 mate Position
        11    :    R2 mate Length
        12    :    R2 mate AS score
        13    :    R2 mate Fragment Middle point
        14    :    R2 mate SNP Matching num
        
        (candidate mate if it is possible)
        
        15    :    Candidate mate Reference
        16    :    Candidate mate Strand
        17    :    Candidate mate Position
        18    :    Candidate mate Length
        19    :    Candidate mate AS score
        20    :    Candidate mate Fragment Middle point
        21    :    Candidate mate SNP Matching num
        22    :    Candidate Index for which mate.
        
    """
    tmp1 = os.path.split(bam1)[1].split('_')
    tmp2 = os.path.split(bam2)[1].split('_')
    out_bam = '_'.join([tmp1[i] for i in range(len(tmp1)) if tmp1[i] == tmp2[i]])
    out_bam = os.path.join(OutPath,out_bam)
    subprocess.call(' '.join(['samtools','merge','-n','-f',out_bam,bam1,bam2,Rebam1,Rebam2]),
                    shell = True)
    
    out_bed = out_bam.replace('.bam','.bed')
    
    out = open(out_bed,'w')
    bam_f = pysam.AlignmentFile(out_bam,'rb')
    
    Unmapped_num = 0
    Multi_num = 0
    Total_num = 0
    tag = ''
    tmp = []
    for read in bam_f:
        name = read.query_name
        name = '_'.join(name.split('_')[:-1])
        if name != tag and len(tmp) == 0:
            tag = name
            tmp.append(read)
        elif name != tag and len(tmp) != 0:
            Total_num += 1
            tag = name
            info = Pair_Integrate(tmp,Frags,Snps,Allelic,level)
            if info == 0:
                Unmapped_num += 1
            elif info == 1:
                Multi_num += 1
            elif info == '':
                Unmapped_num += 1
            else:
                if type(info) == tuple:
                    for i in info:
                        out.writelines('\t'.join(i)+'\n')
                else:
                    out.writelines('\t'.join(info)+'\n')
            tmp = []
            tmp.append(read)
        else:
            tmp.append(read)
    info = Pair_Integrate(tmp,Frags,Snps,Allelic,level)
    Total_num += 1
    if info == 0:
        Unmapped_num += 1
    elif info == 1:
        Multi_num += 1
    elif info == '':
        Unmapped_num += 1
    else:
        if type(info) == tuple:
            for i in info:
                out.writelines('\t'.join(i)+'\n')
        else:
            out.writelines('\t'.join(info)+'\n')
    out.close()
    
    subprocess.call(['rm',out_bam])
    log.log(21,'%s is completed',out_bed)
    return Total_num, Unmapped_num, Multi_num


def Bam_Extract(Bam_Path,Re_Bam_Path,Out_Path,Maternal_Frag,Paternal_Frag,Snp_Path,num,level):
    """
        Merge the two mates results of a pair mapping results Extract the information
        We need and transform To the bed format.
        
        Parameters
        ----------
        Bam_Path : str
            bam Path generated by GlobalMapping module.
        
        Re_Bam_Path : str
            ReMapping Path generated by ReMapping module.
        
        Out_Path : str
            Output Path of bed files.
            
        Maternal_Frag : str
            Maternal Genome Fragment path generate by rebuild genome.
        
        Paternal_Frag : str
            Paternal Genome Fragment path generate by rebuild genome.
        
        Snp_path : str
            SNPs Path generate by rebuild genome.
        
        num : int
            threads number.
        
        level : int
            level of filtering Unique pairs.

    """
    log.log(21,'Loading Snps ...')
    #print 'Loading Snps ...'
    Snps = LoadingSNPs(Snp_Path)
    log.log(21,'Loading Fragments...')
    #print 'Loading Fragments ...'
    M_Frags = LoadFragments(Maternal_Frag)
    P_Frags = LoadFragments(Paternal_Frag)
    log.log(21,'Loading Done!')
    #print 'Loading Done!'
    log.log(21,'Bam information Extraction , Fragments Mapping and SNP matcing for reads ...')
    #print 'Bam information Extraction , Fragments Mapping and SNP matcing for reads ...'
    chunks,chunks_num,Cell = Getchunks(Bam_Path)
    Rechunks,Re_chunks_num,Cell = Getchunks(Re_Bam_Path)
    pool = multiprocessing.Pool(num)
    Statics_M = []
    Statics_P = []
    
    for i in range(chunks_num):
        bam1name = [fil for fil in chunks if 'chunk'+str(i)+'_1_Maternal' in fil][0]
        bam2name = [fil for fil in chunks if 'chunk'+str(i)+'_2_Maternal' in fil][0]
        Rebam1name = [fil for fil in Rechunks if 'chunk'+str(i)+'_1_Maternal' in fil][0]
        Rebam2name = [fil for fil in Rechunks if 'chunk'+str(i)+'_2_Maternal' in fil][0]
        bam1 = os.path.join(Bam_Path,bam1name)
        bam2 = os.path.join(Bam_Path,bam2name)
        Rebam1 = os.path.join(Re_Bam_Path,Rebam1name)
        Rebam2 = os.path.join(Re_Bam_Path,Rebam2name)
        M = pool.apply_async(Bam_Intergrate,args=(bam1,Rebam1,bam2,Rebam2,
                                                  Out_Path,M_Frags,Snps,'Maternal',level))
        
        time.sleep(1)
        bam1name = [fil for fil in chunks if 'chunk'+str(i)+'_1_Paternal' in fil][0]
        bam2name = [fil for fil in chunks if 'chunk'+str(i)+'_2_Paternal' in fil][0]
        Rebam1name = [fil for fil in Rechunks if 'chunk'+str(i)+'_1_Paternal' in fil][0]
        Rebam2name = [fil for fil in Rechunks if 'chunk'+str(i)+'_2_Paternal' in fil][0]
        bam1 = os.path.join(Bam_Path,bam1name)
        bam2 = os.path.join(Bam_Path,bam2name)
        Rebam1 = os.path.join(Re_Bam_Path,Rebam1name)
        Rebam2 = os.path.join(Re_Bam_Path,Rebam2name)
        P = pool.apply_async(Bam_Intergrate,args=(bam1,Rebam1,bam2,Rebam2,
                                                  Out_Path,P_Frags,Snps,'Paternal',level))
        time.sleep(1)
        
        Statics_M.append(M)
        Statics_P.append(P)
        
    time.sleep(2)
    pool.close()
    pool.join()
    
    #Statics
    Total_pair_M = 0;Total_pair_P = 0
    Unmapped_pair_M = 0;Unmapped_pair_P = 0
    Multiple_pair_M = 0;Multiple_pair_P = 0
    
    #Maternal Statics
    for res in Statics_M:
        tmp = res.get()
        Total_pair_M += tmp[0]
        Unmapped_pair_M += tmp[1]
        Multiple_pair_M += tmp[2]
    
    #Paternal Statics
    for res in Statics_P:
        tmp = res.get()
        Total_pair_P += tmp[0]
        Unmapped_pair_P += tmp[1]
        Multiple_pair_P += tmp[2]
    
    Statics_txt = ['Maternal Mapping Statics : ',
                   '===========================',
                   'Total    pair : %d' % Total_pair_M,
                   'Unmapped pair : %d' % Unmapped_pair_M,
                   'Multiple pair : %d' % Multiple_pair_M,
                   'Unique   pair : %d' % (Total_pair_M - Unmapped_pair_M - Multiple_pair_M),
                   '\n',
                   'Paternal Mapping Statics : ',
                   '===========================',
                   'Total    pair : %d' % Total_pair_P,
                   'Unmapped pair : %d' % Unmapped_pair_P,
                   'Multiple pair : %d' % Multiple_pair_P,
                   'Unique   pair : %d' % (Total_pair_P - Unmapped_pair_P - Multiple_pair_P)]
    log.log(21,'\n'+'\n'.join(Statics_txt)+'\n')
    #print '\n'+'\n'.join(Statics_txt)+'\n'






 