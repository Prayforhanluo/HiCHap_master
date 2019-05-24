# -*- coding: utf-8 -*-
"""
Created on Mon May 20 16:21:15 2019

@author: han-luo
"""

from __future__ import division
from scipy import stats
from scipy.stats import ttest_rel
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np
import os, math, bisect, cooler


class LoopAllelicSpecificity(object):
    """
        Calculate the reliability of Allelic Speicficity for input Loops.
    
        The Input Loops can be self-defined.
        Further More:
            U can calculate every single interaction's Allelic Specificity based on
            this class. The level of reliability is based on the all interaction sites.
            So, Change the Input, The reliability of Allelic Specificity for 
            each interaction site may be different.
    
        :: code-block:: python
            

            
    """
    
    def __init__(self, cooler_uri, Loop_file, Res):
        """
            __init__ method
            
        Initializes the Dataset includes M/P haplotype Interaction Matrix and
        Candidate Loops for Allelic Specificity calculation
        
        Input Data
        ----------
        
        
        Loop_file : Loops position file. format as :
                n * rows
                        every rows represents for a candidate loop.
                3 * cols
                        cols 0 : chromosome ID
                        cols 1 : Loop loci 1 for Maternal contact
                        cols 2 : Loop loci 2 for Maternal contact
                        cols 3 : Loop loci 1 for Paternal contact
                        cols 4 : Loop loci 2 for Paternal contact
        
        Res : resolution(int)
            

        """
        self.cooler_fil = '{}::{}'.format(cooler_uri, Res)
        self.cooler = cooler.Cooler(self.cooler_fil)
        self.Loop_file = Loop_file
        self.Res = Res
    
    def DataSetBuilding(self):
        """
             DataSet preparing for model.
        """
        
        
        loop_type = np.dtype({'names':['chr', 'start1', 'end1', 'start2', 'end2'],
                              'formats':['S8', np.int, np.int, np.int, np.int]})
        
        loops = np.loadtxt(self.Loop_file, dtype = loop_type, usecols = [0,1,2,3,4])
        chroms = list(set(loops['chr']))
        Path = os.path.split(self.Loop_file)[0]
        prefix = os.path.split(self.Loop_file)[-1]
        M_Matrix = {}
        P_Matrix = {}
        for chro in chroms:
            M_M = self.cooler.matrix(balance = False).fetch('M'+chro)
            P_M = self.cooler.matrix(balance = False).fetch('P'+chro)
            M_Matrix[chro] = M_M
            P_Matrix[chro] = P_M
            
        
        
        Outfile = os.path.join(Path, 'Allelic_Specificity_'+prefix)
        self.outfile = Outfile        
        
        self.Data = []
        for lp in loops:
            chro = lp['chr']
            start1 = lp['start1'] // self.Res
            end1 = lp['end1'] // self.Res
            start2 = lp['start2'] // self.Res
            end2 = lp['end2'] // self.Res
            M_IF = M_Matrix[chro][start1][end1]
            P_IF = P_Matrix[chro][start2][end2]
            self.Data.append((chro,lp['start1'],lp['end1'],lp['start2'],lp['end2'],M_IF,P_IF))
        
        DataSet_type = np.dtype({'names':['chr','start1','end1','start2','end2','M_IF','P_IF'],
                                 'formats':['S8', np.int, np.int, np.int, np.int, np.float, np.float]})
        
        self.Data = np.array(self.Data, dtype = DataSet_type)
        
        
    def calculate_two_group_stat(self, count, nobs):
        """
            Return stats of two property group test
        """
        p1 = count[0] / nobs[0]
        p2 = count[1] / nobs[1]
        
        p_ = (nobs[0] * p1 + nobs[1] * p2) / (nobs[0] + nobs[1])
        z = (p1 -p2) / math.sqrt((p_ * (1 - p_))*((1 / nobs[0])+(1 / nobs[1])))
        
        return z
    
    def calculate_single_stat(self, p, count, nobs):
        """
            Return stats of Single group property test
        """
        if count == 0 or (nobs - count) == 0:
            return 'NA'
        
        p_ = count / nobs
        
        if p * nobs < 5 or (1 - p) * nobs < 5:
            return 'NA'
            
        elif p * nobs >= 30 and (1 - p) * nobs >= 30:
            t = (nobs * p_ - nobs * p) / math.sqrt(nobs * p * (1 - p))
            
        else:
            t = (abs(nobs * p_ - nobs * p) - 0.5) / math.sqrt(nobs * p * (1 - p))
        
        return t
    
    def calculate_Pvalue(self, stat):
        """
            Return the P-value of z-test of two side
            
        """
        if stat == 'NA':
            return 'NA'
        
        return stats.norm.sf(abs(stat)) * 2
    
    def BH_adjust_pvalue(self, P_value, alpha = 0.05):
        """
            Benjamini-Hochbery procedure for P values.
        """
        Qvalue = multipletests(P_value, alpha = alpha, method = 'fdr_bh')
        
        return Qvalue[1]
    
    def Distribution(self):
        """
            Distribution of Mean.
        """
        Mean = (self.Data['M_IF'] + self.Data['P_IF']) // 2
        Mean = Mean[Mean.nonzero()]
        Mean.sort()
        
        self.Mean = Mean
    
    def BackGround(self):
        """
            BackGround of model
        """
        
        print "Filtering the candidate loops ..."
        self.Distribution()
        nonzero = np.nonzero(self.Mean)
        vmax = np.percentile(nonzero,95)
        mask = ((self.Data['M_IF']+self.Data['P_IF']) / 2 <= vmax) & \
                (self.Data['M_IF'] != 0) & (self.Data['P_IF'] != 0)
        
        self.Data = self.Data[mask]
        
        self.sum_M = self.Data['M_IF'].sum()
        self.sum_T = self.Data['M_IF'].sum() + self.Data['P_IF'].sum()
    
    
    def Running(self):
        """
            Calculate the Allelic Specificity.
        """
        print "=================Allelic Specificity=============="
        print "DataSet preparing ..."
        self.DataSetBuilding()
        print "Background preparing ..."
        self.BackGround()
        
        self.p = self.sum_M / self.sum_T
        print "Maternal total ratio %s." % self.p
        
        print "Modeling ..."
        out = open(self.outfile,'w')
        
        head = ['chr','startM','endM','startP','endP', 'M_IF', 'P_IF', 'QR', 'Log2(FC)', 'stat', 'P_value']
        out.writelines('\t'.join(head)+'\n')
        
        Ratio_P_lst = []
        FC_lst = []
        zstat_lst = []
        pvalue_lst = []
        
        for loop in self.Data:
            loop_M = loop['M_IF']
            loop_T = loop['P_IF'] + loop['M_IF']
            
            stat = self.calculate_single_stat(self.p, loop_M, loop_T)
            pvalue = self.calculate_Pvalue(stat)
            
            if stat == 'NA':
                Ratio_P = 'NA'
                FC = 'NA'
            else:
                loop_mean = loop_T // 2
                Ratio_P = bisect.bisect_left(self.Mean, loop_mean) / len(self.Mean)
                FC = np.log2((loop_M) / (loop_T-loop_M))
            
            Ratio_P_lst.append(Ratio_P)
            FC_lst.append(FC)
            zstat_lst.append(stat)
            pvalue_lst.append(pvalue)
        
        for i, loop in enumerate(self.Data):
            Info = [loop['chr'],loop['start1'],loop['end1'],loop['start2'],loop['end2'],
                    loop['M_IF'],loop['P_IF'],Ratio_P_lst[i],FC_lst[i],zstat_lst[i],pvalue_lst[i]]
                    
            info = map(str,Info)
            
            out.writelines('\t'.join(info)+'\n')
        
        out.close()
        
        print "Done!"



class BoundaryAllelicSpecificity(object):
    """
    """
    
    def __init__(self, cooler_fil, Boundary_fil, Res, offset = 10):
        """
        __init__ method
            
        Initializes the Dataset includes M/P haplotype Interaction Matrix and
        Candidate Loops for Allelic Specificity calculation
        
        Input Data
        ----------
        
        
        Boundary_fil : boundary position file. format as :
                n * rows
                        every rows represents for a candidate loop.
                3 * cols
                        cols 0 : chromosome ID
                        cols 1 : Maternal Boundary
                        cols 2 : Paternal Boundary
        
        Res : resolution(int)
        """
        coolerfil = '{}::{}'.format(cooler_fil, Res)
        self.cooler = cooler.Cooler(coolerfil)
        self.Res = Res
        self.offset = offset
        self.BoundaryFile = Boundary_fil
    
    def LoadingBoundary(self):
        """
        """
        dtype = np.dtype({'names':['chr','pos1','pos2'],
                          'formats':['S4', np.int, np.int]})
        
        self.boundarys = np.loadtxt(self.BoundaryFile, dtype = dtype)
    
    def DataSetBuilding(self):
        """
        """
        print "Loading Boundaries"
        self.LoadingBoundary()
        print "Loading Matrix"
        chroms = list(set(self.boundarys['chr']))
        self.M_Matrix = {}
        self.P_Matrix = {}
        for chro in chroms:
            self.M_Matrix[chro] = self.cooler.matrix(balance=False).fetch('M'+chro)
            self.P_Matrix[chro] = self.cooler.matrix(balance=False).fetch('P'+chro)
        
    def SampleGenerate(self, M, b):
        """
        """
        up = b - self.offset
        down = b + self.offset
        
        upstream = M[up:b, up:b]
        downstream = M[b:down, b:down]
        middlestream = np.tril(M[up:b, b:down])
        
        up_non = upstream[np.nonzero(upstream)]
        down_non = downstream[np.nonzero(downstream)]
        middle_non = middlestream[np.nonzero(middlestream)]
        
        BG = (up_non.sum()+down_non.sum()+middle_non.sum())/(len(up_non)+len(down_non)+len(middle_non))
        
        middlestream = middlestream / BG
        N = middlestream.shape[0] * middlestream.shape[1]
        
        sample = middlestream.reshape((N,))
        
        return sample
    
    
    def removeGap(self, M_S, P_S):
        """
        """
        bool_mask = np.zeros((M_S.shape[0],), dtype = bool)
        for i in range(M_S.shape[0]):
            if M_S[i] != 0 and P_S[i] != 0:
                bool_mask[i] = True
        
        return M_S[bool_mask], P_S[bool_mask]
    
    
    
    def Calculating(self):
        """
        """
        self.DataSetBuilding()
        print "Calculating ..."
        Info = []
        p_value = []
        for bound in self.boundarys:
            chro = bound['chr']
            M_M = self.M_Matrix[chro]
            P_M = self.P_Matrix[chro]
            M_M = M_M - np.diag(M_M.diagonal())
            P_M = P_M - np.diag(P_M.diagonal())
            
            M_pos = bound['pos1'] // self.Res
            P_pos = bound['pos2'] // self.Res
            
            if M_pos == P_pos:
            
                M_S = self.SampleGenerate(M_M, M_pos)
                P_S = self.SampleGenerate(P_M, P_pos)
                
                M_mean = M_S.mean()
                P_mean = P_S.mean()
                if (M_S==0).sum() / len(M_S) >= 0.85:
                    print "Chromosome %s Boundary M:%s -- P:%s surrounded by too much zero bins.. skip ..."  % (chro,bound['pos1'],bound['pos2'])
                    continue
                if (P_S==0).sum() / len(P_S) >= 0.85:
                    print "Chromosome %s Boundary M:%s -- P:%s surrounded by too much zero bins.. skip ..."  % (chro,bound['pos1'],bound['pos2'])
                    continue
                
                M_S, P_S = self.removeGap(M_S, P_S)
                stats, p = ttest_rel(M_S, P_S)
                Info.append((chro, bound['pos1'], bound['pos2'],M_mean, P_mean, stats, p))
                p_value.append(p)
            else:
                M_S1 = self.SampleGenerate(M_M, M_pos)
                P_S1 = self.SampleGenerate(P_M, M_pos)
                M_S2 = self.SampleGenerate(M_M, P_pos)
                P_S2 = self.SampleGenerate(P_M, P_pos)
                if (M_S1==0).sum() / len(M_S1) >= 0.85 or (P_S1==0).sum() / len(P_S1) >= 0.85:
                    if (M_S2==0).sum() / len(M_S2) >= 0.85 or (P_S2==0).sum() / len(P_S2) >= 0.85:
                        print "Chromosome %s Boundary M:%s -- P:%s surrounded by too much zero bins.. skip ..."  % (chro,bound['pos1'],bound['pos2'])
                    else:
                        M_S2, P_S2 = self.removeGap(M_S2, P_S2)
                        stat, p = ttest_rel(M_S2, P_S2)
                        Info.append((chro, bound['pos1'], bound['pos2'],M_mean, P_mean, stats, p))
                        p_value.append(p)
                else:
                    if (M_S2==0).sum() / len(M_S2) >= 0.85 or (P_S2==0).sum() / len(P_S2) >= 0.85:
                        
                        M_S1, P_S1 = self.removeGap(M_S1, P_S1)
                        stat, p = ttest_rel(M_S1, P_S1)
                        Info.append((chro, bound['pos1'], bound['pos2'],M_mean, P_mean, stats, p))
                        p_value.append(p)
                    else:
                        M_S1, P_S1 = self.removeGap(M_S1, P_S1)
                        stat1,p1 = ttest_rel(M_S1, P_S1)
                        M_S2, P_S2 = self.removeGap(M_S2, P_S2)
                        stat2,p2 = ttest_rel(M_S2, P_S2)
                        if p1 < p2:
                            Info.append((chro, bound['pos1'], bound['pos2'],M_S1.mean(), P_S1.mean(), stat1, p1))
                            p_value.append(p1)
                        else:
                            Info.append((chro, bound['pos1'], bound['pos2'],M_S2.mean(), P_S2.mean(), stat2, p2))
                            p_value.append(p2)
                    
        print "BH fdr processing ..."
        Q_value = multipletests(p_value, alpha = 0.05, method = 'fdr_bh')[1]
        AllelicSpecificity = []
        for i in range(len(Info)):
            AllelicSpecificity.append(tuple(list(Info[i]) + [Q_value[i]]))
        
        dtype = np.dtype({'names':['chr','boundary1','boundary2','M_mean','P_mean','stat','p_value','q_value'],
                         'formats':['S4',np.int, np.int, np.float, np.float,np.float,np.float,np.float]})
        
        self.results = np.array(AllelicSpecificity, dtype = dtype)
    
    def OutputTXT(self, outfil):
        """
        """
        print "Output ..."
        head = ['chr','boundaryM','boundaryP','M_mean','P_mean','stat','p_value','q_value']
        
        with open(outfil, 'w') as o:
            o.writelines('\t'.join(head)+'\n')
            
            for line in self.results:
                line = map(str, line)
                o.writelines('\t'.join(line)+'\n')
    
    def Running(self, Outfil):
        """
        """
        
        self.Calculating()
        self.OutputTXT(Outfil)
        
        print 'Done !!!'
        


class CompartmentAllelicSpecificity(object):
    """
    """
    def __init__(self, Maternal_PC, Paternal_PC, Res):
        """
        """
        self.M_PC_file = Maternal_PC
        self.P_PC_file = Paternal_PC
        self.Res = Res
    
    def Loading_PC(self, fil):
        """
        """
        PC = {}
        with open(fil,'r') as f:
            for line in f:
                line = line.strip().split()
                try:
                    PC[line[0]].append(float(line[1]))
                except KeyError:
                    PC[line[0]] = []
                    PC[line[0]].append(float(line[1]))
        
        for k, v in PC.items():
            PC[k] = np.array(v)
        
        return PC
    
    def BackGround(self):
        """
        """
        BG = []
        M_candidate = []
        P_candidate = []
        chroms = self.M_PC.keys()
        for chro in chroms:
            M_PC = self.M_PC[chro]
            P_PC = self.P_PC[chro]
            if np.corrcoef(M_PC, P_PC)[0][1] < 0:
                M_PC = - M_PC
            for i in range(len(M_PC)):
                if M_PC[i] * P_PC[i] < 0:
                    M_candidate.append(M_PC[i])
                    P_candidate.append(P_PC[i])
        
        for i in M_candidate:
            for j in P_candidate:
                BG.append(i - j)
            
        BG.sort()
        
        self.BG = BG
        
        return BG
    
    def Calculating(self):
        """
        """
        print "Loading PC1 values ..."
        self.M_PC = self.Loading_PC(self.M_PC_file)
        self.P_PC = self.Loading_PC(self.P_PC_file)
        
        chroms = self.M_PC.keys()
        Info = []
        P_Value = []
        BG = self.BackGround()
        for chro in chroms:
            print "chroms %s start" % chro
            M_PC = self.M_PC[chro]
            P_PC = self.P_PC[chro]
            if np.corrcoef(M_PC, P_PC)[0][1] < 0:
                M_PC = - M_PC
            for i in range(len(M_PC)):
                if M_PC[i] * P_PC[i] >= 0:
                    continue
                diff = M_PC[i] - P_PC[i]
                
                index_forward = bisect.bisect_left(BG, diff)
                index_reverse = len(BG) - index_forward
                
                index = min(index_forward, index_reverse)
                
                pvalue = index / len(BG)
                    
                Info.append((chro, i *self.Res,M_PC[i], P_PC[i], diff, pvalue))
                P_Value.append(pvalue)
        
        print "BH fdr processing ..."
        Q_value = multipletests(P_Value, method = 'fdr_bh')[1]
        
        AllelicSpecificity = []
        for i in range(len(Info)):
            AllelicSpecificity.append(tuple(list(Info[i]) + [Q_value[i]]))
        
        dtype = np.dtype({'names':['chr', 'pos', 'pc-m','pc-p','diff','p_value', 'Q_value'],
                          'formats':['S4', np.int, np.float, np.float, np.float, np.float, np.float]})
        
        AllelicSpecificity = np.array(AllelicSpecificity, dtype = dtype)
        
        self.results = AllelicSpecificity
    
    def OutputTXT(self, outfil):
        """
        """
        print "Output ..."
        
        head = ['chr','position','PC-M','PC-P','diff','P_Value','Q_Value']
        with open(outfil, 'w') as out:
            out.writelines('\t'.join(head)+'\n')
            
            for line in self.results:
                line = map(str, line)
                out.writelines('\t'.join(line)+'\n')
    
    def Running(self, Outfil):
        """
        """
        self.Calculating()
        self.OutputTXT(Outfil)
        
        
                
                
            
            
   
        