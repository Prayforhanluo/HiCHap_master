# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 19:31:12 2018

@author: han-luo
"""

from __future__ import  division
#from sinkhorn_knopp import sinkhorn_knopp as skp
import numpy as np
from mirnylib.numutils import completeIC
import os, time, logging

log = logging.getLogger(__name__)


def Loading_Matrix(NPZ_Path):
    """
    """
    
    return np.load(NPZ_Path)
    

    
def ICE_Normalize(M):
    """
     ICE Normalize from mirnylib
    """
    cM, bias = completeIC(M, returnBias = True)
    
    return cM, bias


def bias_handle(bias):
    """
        
    """
    bias = bias.reshape(bias.shape[0],)
    
    bias[np.isnan(bias)] = 1.0
    
    return bias


def Coverage_M(Matrix):
    """
    """
    cover = []
    for i in Matrix:
        tmp = 1 - ((i == 0).sum() / float(len(i)))
        cover.append(tmp)
    
    return np.array(cover)


def Gap_defined(Matrix):
    """
    """
    Cover = Coverage_M(Matrix)
    
    threshold = np.percentile(Cover[np.nonzero(Cover)],25)
    if threshold > 0.3:
        threshold = 0.3
    gap = []
    for index,i in enumerate(Matrix):
        tmp = 1 - (i == 0).sum() / float(len(i))
        if tmp < threshold:
            gap.append(index)
                        
    return np.array(gap)


def Non_Gap_Defined(N,gap):
    """
        Defined the chromosome State.
        chromosome -> Gap region
                      Non Gap region                
    """ 
    Non_Gap = []
    for i in range(N):
        if i not in gap:
            Non_Gap.append(i)
    
    return np.array(Non_Gap)

def Trans2symmetry(Matrix,gap):
    """
    """
    if gap.size == 0:
        upper_M = np.triu(Matrix) + np.tril(Matrix,-1).T
        Sym_M = np.triu(upper_M,1).T + upper_M
        
        return Sym_M
        
    N = Matrix.shape[0]
    New_Matrix = np.zeros(Matrix.shape)
    Non_Gap = Non_Gap_Defined(N,gap)
      
    # Gap <-> NonGap  &   Gap <-> Gap
    for i in gap:
        for j in range(N):
            if i == j:
                New_Matrix[i][j] = Matrix[i][j]
            else:
                value = max(Matrix[i][j], Matrix[j][i])
                New_Matrix[i][j] = value
                New_Matrix[j][i] = value

    # NonGap <-> NonGap    
    for i in Non_Gap:
        for j in range(N):
            if i == j:
                New_Matrix[i][j] = Matrix[i][j]
            else:
                value = (Matrix[i][j] + Matrix[j][i]) / 2.0
                New_Matrix[i][j] = value
                New_Matrix[j][i] = value
    
    return New_Matrix
        

def Correct_VC(X, alpha):
    x = np.array(X,float)
    s1 = np.sum(x, axis = 1)
    s1 = s1 ** alpha
#    s /= np.mean(s[s!=0])
    s1[s1 == 0] = 1
    s2 = np.sum(x, axis = 0)
    s2 = s2 ** alpha
#    s2 /= np.mean(s2[s2 != 0])
    s2[s2 == 0] = 1
    return x / (s2[None, :] * s1[:, None])


def Normal_VC_Correct(NPZ):
    """
    """
    Raw_Lib = np.load(NPZ)
    Nor_Lib = {}
    for c in Raw_Lib.keys():
        Raw_M = Raw_Lib[c]
        Nor_M = Correct_VC(Raw_M, 2/3)
        
        #Recall
        Re_factor = Raw_M.mean() / Nor_M.mean()
        Nor_M = Nor_M * Re_factor
        Nor_Lib[c] = Nor_M
    
    return Nor_Lib




def Alpha_Allelic_Correct(Raw_NPZ,M_NPZ,P_NPZ,chroms = ['#','X']):
    """
    """
    Raw_Lib = np.load(Raw_NPZ)
    M_Lib = np.load(M_NPZ)
    P_Lib = np.load(P_NPZ)
    Nor_M = {}
    Nor_M_Matrix = {}
    Nor_M_Gap = {}
    Nor_P = {}
    Nor_P_Matrix = {}
    Nor_P_Gap = {}
    for chro in M_Lib.keys():
        if (not chroms) or (chro.isdigit() and '#' in chroms) or (chro in chroms):
            Raw_M = Raw_Lib[chro]
            M_M = M_Lib[chro]
            P_M = P_Lib[chro]
            alpha = []
            N = Raw_M.shape[0]
            Gap_M = Gap_defined(M_M)
            Gap_P = Gap_defined(P_M)
            Non_Gap_M = Non_Gap_Defined(N,Gap_M)
            Non_Gap_P = Non_Gap_Defined(N,Gap_P)
        
            for i in range(N):
                alpha.append((M_M[i].sum() + P_M[i].sum()) / (Raw_M[i].sum() + 1))
           
            alpha = np.array(alpha)
            # remove Gap effect
            Non_Gap_index = list(set(Non_Gap_M) | set(Non_Gap_P))
        
            #SNP Correcting
            alpha /= np.max(alpha[Non_Gap_index])
            alpha[alpha == 0] = 1
            threshold = np.percentile(alpha[Non_Gap_index],30)
            alpha[alpha < threshold] = threshold
        
            S_M_M = M_M / alpha[:, None]
            S_P_M = P_M / alpha[:, None]
        
            #Non-Symmetry Maxtrinx To Symmetry Matrix
            Sym_M_M = Trans2symmetry(S_M_M,Gap_M)
            Sym_M_P = Trans2symmetry(S_P_M,Gap_P)
        
            #HiC bias Correcting
            Cor_Sym_M = Correct_VC(Sym_M_M,2/3)
            Cor_Sym_P = Correct_VC(Sym_M_P,2/3)

            #Recall Factor
            M_R_F = M_M.mean() / Cor_Sym_M.mean()
            P_R_F = P_M.mean() / Cor_Sym_P.mean()        
        
            Nor_M_Matrix[chro] = Cor_Sym_M * M_R_F
            Nor_M_Gap[chro] = Gap_M
            Nor_P_Matrix[chro] = Cor_Sym_P * P_R_F
            Nor_P_Gap[chro] = Gap_P        
            print chro
        
    Nor_M['Matrix'] = Nor_M_Matrix
    Nor_M['Gap'] = Nor_M_Gap
    Nor_P['Matrix'] = Nor_P_Matrix
    Nor_P['Gap'] = Nor_P_Gap
    
    return Nor_M, Nor_P 
    

def Whole_Normalize(OutPath, NPZ_Path,chroms = ['#','X']):
    """
    """
    whole_M_file = [i for i in os.listdir(NPZ_Path) if 'Whole_Genome' in i]
    prefix_lst = [i.split('Whole')[0] for i in whole_M_file]
    
    for fil in whole_M_file:
        outfil = os.path.join(OutPath,'Correct_'+fil)
        Correct_Matrix = {}        
        lib = Loading_Matrix(os.path.join(NPZ_Path,fil))
        Correct_Matrix['Bins'] = lib['Bins'][()]
        M = lib['Matrix'][()]
        
        Correct_M, bias = ICE_Normalize(M)
        Correct_Matrix['Matrix'] = Correct_M
        Correct_Matrix['bias'] = bias
        
        np.savez_compressed(outfil, **Correct_Matrix)
        #print "%s completed" % outfil
        log.log(21,'%s completed' % outfil)
        
        del Correct_Matrix
        time.sleep(3)
    
    local_M_file = [i for i in os.listdir(NPZ_Path) if 'Local_Chromosome' in i]
    
    for fil in local_M_file:
        outfil = os.path.join(OutPath,'Correct_'+fil)
        Correct_Matrix = {}
        Correct_Matrix['Matrix'] = {}
        Correct_Matrix['bias'] = {}
        
        lib = Loading_Matrix(os.path.join(NPZ_Path,fil))
        
        for c in lib.keys():
            if (not chroms) or (c.isdigit() and '#' in chroms) or (c in chroms):                
                M = lib[c]
                Correct_M, bias = ICE_Normalize(M)
                Correct_Matrix['Matrix'][c] = Correct_M           
                Correct_Matrix['bias'][c] = bias
        np.savez_compressed(outfil, **Correct_Matrix)
        #print "%s completed" % outfil
        log.log(21,'%s completed' % outfil)
        
        del Correct_Matrix
        time.sleep(3)
        
    Allelic_M_file = [i for i in os.listdir(NPZ_Path) if '_Mate_' in i]
    
    for prefix in prefix_lst:
        # Both Mate Matrix Correct
        for fil in Allelic_M_file:
            if 'Both' in fil and prefix in fil:
                outfil = os.path.join(OutPath,'Correct_' + fil)
                raw_npz = os.path.join(NPZ_Path,fil)
                Nor_Lib = Normal_VC_Correct(raw_npz)
                np.savez_compressed(outfil,**Nor_Lib)
                log.log(21,'%s completed', outfil)
        
        # One Mate Matrix Correct
        Maternal_npz = ''
        Paternal_npz = ''
        Raw_npz = ''
        
        for fil in Allelic_M_file:
            if 'One' in fil and prefix in fil and 'Maternal' in fil:
                Maternal_npz = os.path.join(NPZ_Path,fil)
                Correct_M_fil = 'Correct_' + fil
            if 'One' in fil and prefix in fil and 'Paternal' in fil:
                Paternal_npz = os.path.join(NPZ_Path,fil)
                Correct_P_fil = 'Correct_' + fil
        for fil in local_M_file:
            if prefix in fil:
                Raw_npz = os.path.join(NPZ_Path,fil)
        
        if Maternal_npz == '' or Paternal_npz == '' or Raw_npz == '':
            raise Exception ("some files lost, pls do not change its name")
        
        
        Correct_M_Matrix, Correct_P_Matrix = Alpha_Allelic_Correct(Raw_npz,Maternal_npz,Paternal_npz)
        np.savez_compressed(os.path.join(OutPath,Correct_M_fil),**Correct_M_Matrix)
        np.savez_compressed(os.path.join(OutPath,Correct_P_fil),**Correct_P_Matrix)

        log.log(21,"%s completed", os.path.join(OutPath,Correct_M_fil))
        log.log(21,"%s completed" % os.path.join(OutPath,Correct_P_fil))

        time.sleep(3)

def Allelic_Single_VC_Correcting(npz,out):
    """
    """
    Correct_fil = os.path.join(out,'Correct_'+os.path.split(npz)[-1])
    Normal_lib = Normal_VC_Correct(npz)
    np.savez_compressed(Correct_fil,**Normal_lib)
    log.log(21,'%s completed' % Correct_fil)


def Allelic_Single_Alpha_Correcting(raw,Maternal,Paternal,out):
    """
    """
    Correct_M_fil = os.path.join(out,'Correct_'+os.path.split(Maternal)[-1])
    Correct_P_fil = os.path.join(out,'Correct_'+os.path.split(Paternal)[-1])
    
    Correct_M_Matrix, Correct_P_Matrix = Alpha_Allelic_Correct(raw,Maternal,Paternal)
    np.savez_compressed(Correct_M_fil,**Correct_M_Matrix)
    np.savez_compressed(Correct_P_fil,**Correct_P_Matrix)
    log.log(21,'%s completed' % Correct_M_fil)
    log.log(21,'%s completed' % Correct_P_fil)
    
          
def Non_Allelic_Single_Correcting(npz,out,Mode,chroms = ['#','X'], save_sparse = True):
    """
    """
    lib = np.load(npz)
    outfil = os.path.join(out,'Correct_'+os.path.split(npz)[1])
    Correct_Matrix = {}
    
    dtype = np.dtype({'names':['bin1','bin2','IF'],
                      'formats':[np.int, np.int, np.float]})
    
    if Mode == 'whole':
        Correct_Matrix['Bins'] = lib['Bins'][()]
        M = lib['Matrix'][()]
        
        Correct_M, bias = ICE_Normalize(M)
        
#        Correct_Matrix['Matrix'] = Correct_M
        Correct_Matrix['bias'] = bias
        if save_sparse:
            Triu = np.triu(Correct_M)
            x, y = np.nonzero(Triu)
            values = Triu[x, y]
            tmp = np.zeros(values.size, dtype = dtype)
            tmp['bin1'] = x
            tmp['bin2'] = y
            tmp['IF'] = values
            Correct_Matrix['Matrix'] = tmp
        else:
            Correct_Matrix['Matrix'] = Correct_M
        
        
    else:
        Correct_Matrix['Matrix'] = {}
        Correct_Matrix['bias'] = {}
        
        for c in lib.keys():
            if (not chroms) or (c.isdigit() and '#' in chroms) or (c in chroms):
                M = lib[c]
                Correct_M, bias = ICE_Normalize(M)
                Correct_Matrix['bias'][c] = bias
                if save_sparse:
                    Triu = np.triu(Correct_M)
                    x, y = np.nonzero(Triu)
                    values = Triu[x, y]
                    tmp = np.zeros(values.size,dtype = dtype)
                    tmp['bin1'] = x
                    tmp['bin2'] = y
                    tmp['IF'] = values
                    Correct_Matrix['Matrix'][c] = tmp
                else:
                    Correct_Matrix['Matrix'][c] = Correct_M
    
    np.savez_compressed(outfil, **Correct_Matrix)
    log.log(21,'%s completed' % outfil)
