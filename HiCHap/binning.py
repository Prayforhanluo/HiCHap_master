# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 19:31:20 2018

@author: han-luo
"""

from __future__ import division
import numpy as np
import os, logging, subprocess

log = logging.getLogger(__name__)


def properU(res):
    """
        Express a genomic position in a proper unit (KB, MB, or both) 

    """
    
    i_part = int(res) // 1000000 # Integer Part
    d_part = (int(res) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])



def Merge_beds(bed_lst):
    """
        Merging the  beds To Generate the Matraix.
    """
    Merge_cmd = ['cat'] + list(bed_lst)
    
    return Merge_cmd


    
def Load_Genome(genomeSize):
    """
        Loading the Genome Size
        Return the genomeSize dict.
    """
    
    genome = {}
    with open(genomeSize,'r') as f:
        for line in f:
            line = line.strip().split()
            c = line[0].lstrip('chr')
            if (not c.isdigit() and c != 'X'):
                continue
            else:
                genome[c] = int(line[1])
    
    return genome



def Sort_Chromosomes(chro_lst):
    """
        Sort the chromosome.
    """
    chro_lst = [i.lstrip('chr') for i in chro_lst]
    
    Num_chro = []
    Str_chro = []
    for i in chro_lst:
        try:
            i = int(i)
            Num_chro.append(i)
        except:
            Str_chro.append(i)
    Num_chro.sort()
    Num_chro = [str(j) for j in Num_chro]
    Str_chro.sort()
        
    return Num_chro + Str_chro



def Get_Chro_Bins(genomeSize,Resolution):
    """
    """
    ChromBins = Load_Genome(genomeSize)
    Bins = {}
    Chro = Sort_Chromosomes(ChromBins)
    for c,l in ChromBins.items():
        ChromBins[c] = l // Resolution
    
    for index,chro in enumerate(Chro):
        if index == 0:
            Bins[chro] = (0,ChromBins[chro])
        else:
            Bins[chro] = (Bins[Chro[index - 1]][1] +1, ChromBins[chro] + Bins[Chro[index - 1]][1]+1)
    
    Bins_num = Bins[Chro[-1]][1]
    
    return Bins, Bins_num + 1


def NormalMartixBuilding_Non_Allelic(bed_IO, genomeSize, Whole_Res, Local_Res):
    """
        Building the contact Matrix for whole Genome with Non-Allelic Way.
        
    """
    
    Bins_Pos, Bins_Sum = Get_Chro_Bins(genomeSize,Whole_Res)
    genome = Load_Genome(genomeSize)
    Bychro_Matrix = {}
    
    for c,l in genome.items():
        bin_size = (l // Local_Res) + 1
        Bychro_Matrix[c] = np.zeros((bin_size,bin_size),dtype = int)
        
    Whole_Matrix = {}
    Whole_Matrix['Bins'] = Bins_Pos
    Matrix = np.zeros((Bins_Sum,Bins_Sum),dtype = np.int)
    
    count = 0
    for line in bed_IO:
        count += 1
        if count % 20000000 == 0:
            #print "%d pairs completed" % count
            log.log(21,"%d pairs completed" % count)
        line = line.strip().split()
        
        #WholeGenome Matrix Building
        c1 = line[1].lstrip('chr')
        if (not c1.isdigit() and c1 != 'X'):
            continue
        s1 = int(line[6]) // Whole_Res
        bin1 = s1 + Bins_Pos[c1][0]
        
        c2 = line[8].lstrip('chr')
        if (not c2.isdigit() and c2 != 'X'):
            continue
        s2 = int(line[13]) // Whole_Res
        bin2 = s2 + Bins_Pos[c2][0]
        
        if bin1 != bin2:
            Matrix[bin1][bin2] += 1
            Matrix[bin2][bin1] += 1
        else:
            Matrix[bin1][bin2] += 1
        
        #ByChromosome Matrix Building
        if c1 == c2:
            pos1 = int(line[6]) // Local_Res
            pos2 = int(line[13]) // Local_Res
            if pos1 != pos2:
                Bychro_Matrix[c1][pos1][pos2] += 1
                Bychro_Matrix[c1][pos2][pos1] += 1
            else:
                Bychro_Matrix[c1][pos1][pos2] += 1
    
    Whole_Matrix['Matrix'] = Matrix
    
    return Whole_Matrix, Bychro_Matrix


def Bed_To_Normal_Matrix_Non_Allelic(OutPath, genomeSize, Whole_Res, Local_Res, Rep_Path):
    """
        Merge the Replicates and Created the Interaction Matrix.
        
        Rep_Path : Replicates bed Path list.        
    """
    log.log(21,"Building Replicates Matrix respectively")
    
    Merge_All_bed = []
    for rep_p in Rep_Path:
        files = [i for i in os.listdir(rep_p) if '_Valid.bed' in i]
        prefix = files[0].split('Valid')[0]
        files = [os.path.join(rep_p,fil) for fil in files]
        
        Merge_All_bed += files
        
        Merge_cmd = Merge_beds(files)
        
        pread = subprocess.Popen(Merge_cmd, stdout=subprocess.PIPE, bufsize= -1)
        
        instream = pread.stdout
        
        whole_M, local_M = NormalMartixBuilding_Non_Allelic(instream,genomeSize,Whole_Res,Local_Res)
        
        outfil_whole = os.path.join(OutPath,prefix+'Whole_Genome_Matrix.npz')
        outfil_local = os.path.join(OutPath,prefix+'Local_Chromosome_Matrix.npz')
        
        np.savez_compressed(outfil_whole,**whole_M)
        np.savez_compressed(outfil_local,**local_M)
        
        instream.close()
        
        log.log(21,"%s completed" % outfil_whole)
        log.log(21,"%s completed" % outfil_local)
        
        del whole_M, local_M
        
    log.log(21,'Merge Replicates and Building Matraix')
    
    prefix = 'Merged_Reps_'
    #whole_genome
    whole_M = {}
    Reps_fil = [fils for fils in os.listdir(OutPath) if 'Whole_Genome_Matrix' in fils]
    for fil in Reps_fil:
        fil = os.path.join(OutPath,fil)
        Lib = np.load(fil)
        if 'Bins' not in whole_M.keys():
            whole_M['Bins'] = Lib['Bins'][()]
        if 'Matrix' not in whole_M.keys():
            whole_M['Matrix'] = Lib['Matrix']
        else:
            whole_M['Matrix'] += Lib['Matrix']    
    outfil_whole = os.path.join(OutPath,prefix+'Whole_Genome_Matrix.npz')
    np.savez_compressed(outfil_whole,**whole_M)
    log.log(21,"%s completed" % outfil_whole)
    
    #local_chromosome
    local_M = {}
    Reps_fil = [fils for fils in os.listdir(OutPath) if 'Local_Chromosome_Matrix' in fils]
    for fil in Reps_fil:
        fil = os.path.join(OutPath,fil)
        Lib = np.load(fil)
        for c in Lib.keys():
            if c not in local_M.keys():
                local_M[c] = Lib[c]
            else:
                local_M[c] += Lib[c]
    
    outfil_local = os.path.join(OutPath,prefix+'Local_Chromosome_Matrix.npz')
    np.savez_compressed(outfil_local,**local_M)
    
    log.log(21,"%s completed" % outfil_local)
    
    del whole_M,local_M
    
    log.log(21,"Matrix Building Done!")



def NormalMartixBuilding(bed_IO,genomeSize,Whole_Res,Local_Res):
    """
        Bulding the Contact Matrix for Whole Genome.
    """
    
    
    Bins_Pos, Bins_Sum = Get_Chro_Bins(genomeSize,Whole_Res)
    genome = Load_Genome(genomeSize)
    Bychro_Matrix = {}
    
    for c,l in genome.items():
        bin_size = (l // Local_Res) + 1
        Bychro_Matrix[c] = np.zeros((bin_size,bin_size),dtype = int)
        
    Whole_Matrix = {}
    Whole_Matrix['Bins'] = Bins_Pos
    Matrix = np.zeros((Bins_Sum,Bins_Sum),dtype = np.int)
    
    count = 0
    for line in bed_IO:
        count += 1
        if count % 20000000 == 0:
            #print "%d pairs completed" % count
            log.log(21,"%d pairs completed" % count)
        line = line.strip().split()
        
        #WholeGenome Matrix Building
        c1 = line[0].lstrip('chr')
        if (not c1.isdigit() and c1 != 'X'):
            continue
        s1 = int(line[1]) // Whole_Res
        bin1 = s1 + Bins_Pos[c1][0]
        
        c2 = line[2].lstrip('chr')
        if (not c2.isdigit() and c2 != 'X'):
            continue
        s2 = int(line[3]) // Whole_Res
        bin2 = s2 + Bins_Pos[c2][0]
        
        if bin1 != bin2:
            Matrix[bin1][bin2] += 1
            Matrix[bin2][bin1] += 1
        else:
            Matrix[bin1][bin2] += 1
        
        #ByChromosome Matrix Building
        if c1 == c2:
            pos1 = int(line[1]) // Local_Res
            pos2 = int(line[3]) // Local_Res
            if pos1 != pos2:
                Bychro_Matrix[c1][pos1][pos2] += 1
                Bychro_Matrix[c1][pos2][pos1] += 1
            else:
                Bychro_Matrix[c1][pos1][pos2] += 1
    
    Whole_Matrix['Matrix'] = Matrix
    
    return Whole_Matrix, Bychro_Matrix



def AllelicMatrixBuilding(bed_IO,genomeSize,Resolution):
    """
        Building the Contact Matrix for Genome by Chromosome.
        

    """
    
    genome = Load_Genome(genomeSize)
    One_Mate_IF = {}
    Both_Mate_IF = {}
    for c, l in genome.items():
        bin_size = (l // Resolution) + 1
        One_Mate_IF[c] = np.zeros((bin_size,bin_size),dtype = int)
        Both_Mate_IF[c] = np.zeros((bin_size,bin_size),dtype = int)
    
    count = 0
    for line in bed_IO:
        count += 1
        if count % 20000000 == 0:
            #print "%d pairs completed" % count
            log.log(21,"%d pairs completed" % count)
        line = line.split()
        if line[0] != line[2]:
            continue
        else:
            c = line[0].lstrip('chr')
            if (not c.isdigit() and c != 'X'):
                continue
            bin1 = int(line[1]) // Resolution
            bin2 = int(line[3]) // Resolution
            mark = line[-1]
            if bin1 != bin2:
                if mark == 'Both':
                    One_Mate_IF[c][bin1][bin2] += 1
                    One_Mate_IF[c][bin2][bin1] += 1
                    Both_Mate_IF[c][bin1][bin2] += 1
                    Both_Mate_IF[c][bin2][bin1] += 1
                elif mark == 'R1':
                    One_Mate_IF[c][bin1][bin2] += 1
                elif mark == 'R2':
                    One_Mate_IF[c][bin2][bin1] += 1
                else:
                    raise Exception ("Unkown Allelic Mark %s" % mark)
            else:
                if mark == 'Both':
                    One_Mate_IF[c][bin1][bin2] += 1
                    Both_Mate_IF[c][bin1][bin2] += 1
                elif mark == 'R1':
                    One_Mate_IF[c][bin1][bin2] += 1
                elif mark == 'R2':
                    One_Mate_IF[c][bin2][bin1] += 1
                else:
                    raise Exception ("Unkown Allelic Mark %s" % mark)
            
    
    return One_Mate_IF, Both_Mate_IF


def Check_Bed(bed_lst):
    """
        Check whether the necessary bed files are available.
    """
    mark = True
    if len([i for i in bed_lst if 'Bi_Allelic.bed' in i]) > 0:
        pass
    else:
        return False, 'Bi_Allelic'
    
    if len([i for i in bed_lst if 'M_M.bed' in i]) > 0:
        pass
    else:
        return False, 'M_M'
    
    if len([i for i in bed_lst if 'P_P.bed' in i]) > 0:
        pass
    else:
        return False, 'P_P'
    
    if len([i for i in bed_lst if 'M_P.bed' in i]) > 0:
        pass
    else:
        return False, 'M_P'
    
    if len([i for i in bed_lst if 'P_M.bed' in i]) > 0:
        pass
    else:
        return False, 'P_M'
    
    return mark, ''



def Bed_To_Allelic_Matrix(OutPath, genomeSize, Whole_Res, Local_Res, Rep_Path):
    """
        Merge the Replicates and Created the Interaction Matrix.
        
        Rep_Path : Replicates bed Path list.
    """
    
    #-----------------------Replicate Matrix Building-------------------------------   
    #print "Building Replicates Matrix respectively."
    log.log(21,"Building Replicates Matrix respectively")
    
    Merge_All_bed = []
    Merge_M_bed = []
    Merge_P_bed = []
    for rep_p in Rep_Path:
        files = [i for i in os.listdir(rep_p) if 'Bi_Allelic.bed' in i or \
                'M_M.bed' in i or 'M_P.bed' in i or 'P_P.bed' in i or 'P_M.bed' in i]
        
        files.sort()
        prefix = files[0].split('Valid')[0]
        
        if len(files) != 5:
            mark,tmp = Check_Bed(files)
            if mark == True:
                pass
            else:
                raise Exception ("Missing file %s.bed in %s" % (tmp,rep_p))
        
        files = [os.path.join(rep_p,fil) for fil in files]
        
        Merge_All_bed += files
        Merge_cmd = Merge_beds(files)
        #print '\t'.join(Merge_cmd)
        pread = subprocess.Popen(Merge_cmd, stdout=subprocess.PIPE, bufsize = -1)
        
        instream = pread.stdout
        whole_M, local_M = NormalMartixBuilding(instream,genomeSize,Whole_Res,Local_Res)
        
        outfil_whole = os.path.join(OutPath,prefix+'Whole_Genome_Matrix.npz')
        outfil_local = os.path.join(OutPath,prefix+'Local_Chromosome_Matrix.npz')
        
        np.savez_compressed(outfil_whole,**whole_M)
        np.savez_compressed(outfil_local,**local_M)
        instream.close()
        #print "%s completed" % outfil_whole
        #print "%s completed" % outfil_local
        log.log(21,"%s completed" % outfil_whole)
        log.log(21,"%s completed" % outfil_local)
        
        del whole_M, local_M
        
        M_fil = [i for i in files if 'M_M.bed' in i]
        
        Merge_M_bed += M_fil
        Merge_cmd = Merge_beds(M_fil)
        #print '\t'.join(Merge_cmd)
        pread = subprocess.Popen(Merge_cmd, stdout=subprocess.PIPE, bufsize = -1)
        
        instream = pread.stdout
        One_M, Both_M = AllelicMatrixBuilding(instream,genomeSize,Local_Res)
        
        outfil_One = os.path.join(OutPath,prefix+'One_Mate_Maternal.npz')
        outfil_Both = os.path.join(OutPath,prefix+'Both_Mate_Maternal.npz')
        
        np.savez_compressed(outfil_One,**One_M)
        np.savez_compressed(outfil_Both, **Both_M)
        instream.close()
        #print "%s completed" % outfil_One
        #print "%s completed" % outfil_Both
        log.log(21,"%s completed" % outfil_One)
        log.log(21,"%s completed" % outfil_Both)
        del One_M, Both_M
                          
        
        P_fil = [i for i in files if 'P_P.bed' in i]
        
        Merge_P_bed += P_fil
        Merge_cmd = Merge_beds(P_fil)
        #print '\t'.join(Merge_cmd)
        pread = subprocess.Popen(Merge_cmd, stdout=subprocess.PIPE, bufsize = -1)
        
        instream = pread.stdout
        One_M, Both_M = AllelicMatrixBuilding(instream,genomeSize,Local_Res)
        
        outfil_One = os.path.join(OutPath,prefix+'One_Mate_Paternal.npz')
        outfil_Both = os.path.join(OutPath,prefix+'Both_Mate_Paternal.npz')
        
        np.savez_compressed(outfil_One,**One_M)
        np.savez_compressed(outfil_Both, **Both_M)
        instream.close()
        #print "%s completed" % outfil_One
        #print "%s completed" % outfil_Both
        log.log(21,"%s completed" % outfil_One)
        log.log(21,"%s completed" % outfil_Both)
        del One_M, Both_M
    
    #print "Merge Replicates and Building Matrix."
    log.log(21,"Merge Replicates ana Building Matraix")

    
    #--------------------------Tranditinal Matrix-----------------------------------
    if len(Rep_Path) == 1:
        log.log(21,'Only one Replicate was found...')
        log.log(21,'Skip and Down.')
        return
    prefix = 'Merged_Reps_'   
    
    #whole_genome
    whole_M = {}
    Reps_fil = [fils for fils in os.listdir(OutPath) if 'Whole_Genome_Matrix' in fils]
    for fil in Reps_fil:
        fil = os.path.join(OutPath,fil)
        Lib = np.load(fil)
        if 'Bins' not in whole_M.keys():
            whole_M['Bins'] = Lib['Bins'][()]
        if 'Matrix' not in whole_M.keys():
            whole_M['Matrix'] = Lib['Matrix']
        else:
            whole_M['Matrix'] += Lib['Matrix']    
    outfil_whole = os.path.join(OutPath,prefix+'Whole_Genome_Matrix.npz')
    np.savez_compressed(outfil_whole,**whole_M)
    log.log(21,"%s completed" % outfil_whole)
    
    #local_chromosome
    local_M = {}
    Reps_fil = [fils for fils in os.listdir(OutPath) if 'Local_Chromosome_Matrix' in fils]
    for fil in Reps_fil:
        fil = os.path.join(OutPath,fil)
        Lib = np.load(fil)
        for c in Lib.keys():
            if c not in local_M.keys():
                local_M[c] = Lib[c]
            else:
                local_M[c] += Lib[c]
    
    outfil_local = os.path.join(OutPath,prefix+'Local_Chromosome_Matrix.npz')
    np.savez_compressed(outfil_local,**local_M)
    
    log.log(21,"%s completed" % outfil_local)
    
    del whole_M,local_M

    #--------------------------Allelic Maternal Matrix------------------------------
    
    One_M = {}
    Reps_fil = [fils for fils in os.listdir(OutPath) if 'One_Mate_Maternal.npz' in fils]
    for fil in Reps_fil:
        fil = os.path.join(OutPath,fil)
        Lib = np.load(fil)
        for c in Lib.keys():
            if c not in One_M.keys():
                One_M[c] = Lib[c]
            else:
                One_M[c] += Lib[c]      
    outfil_One = os.path.join(OutPath,prefix+'One_Mate_Maternal.npz')
    np.savez_compressed(outfil_One,**One_M)
    log.log(21,"%s completed" % outfil_One)
    
    Both_M = {}
    Reps_fil = [fils for fils in os.listdir(OutPath) if 'Both_Mate_Maternal' in fils]
    for fil in Reps_fil:
        fil = os.path.join(OutPath,fil)
        Lib = np.load(fil)
        for c in Lib.keys():
            if c not in Both_M.keys():
                Both_M[c] = Lib[c]
            else:
                Both_M[c] += Lib[c]    
    outfil_Both = os.path.join(OutPath,prefix+'Both_Mate_Maternal.npz')
    np.savez_compressed(outfil_Both, **Both_M)
    log.log(21,"%s completed" % outfil_Both)
    
    del One_M, Both_M
    
    
    #--------------------------Allelic Paternal Matrix------------------------------
    One_M = {}
    Reps_fil = [fils for fils in os.listdir(OutPath) if 'One_Mate_Paternal.npz' in fils]
    for fil in Reps_fil:
        fil = os.path.join(OutPath,fil)
        Lib = np.load(fil)
        for c in Lib.keys():
            if c not in One_M.keys():
                One_M[c] = Lib[c]
            else:
                One_M[c] += Lib[c]      
    outfil_One = os.path.join(OutPath,prefix+'One_Mate_Paternal.npz')
    np.savez_compressed(outfil_One,**One_M)
    log.log(21,"%s completed" % outfil_One)
    
    Both_M = {}
    Reps_fil = [fils for fils in os.listdir(OutPath) if 'Both_Mate_Paternal' in fils]
    for fil in Reps_fil:
        fil = os.path.join(OutPath,fil)
        Lib = np.load(fil)
        for c in Lib.keys():
            if c not in Both_M.keys():
                Both_M[c] = Lib[c]
            else:
                Both_M[c] += Lib[c]    
    outfil_Both = os.path.join(OutPath,prefix+'Both_Mate_Paternal.npz')
    np.savez_compressed(outfil_Both, **Both_M)
    log.log(21,"%s completed" % outfil_Both)
    
    del One_M, Both_M
    
        
    
    
    
    