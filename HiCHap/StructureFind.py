# -*- coding: utf-8 -*-
"""
Created on Sat Feb 23 16:08:13 2019

@author: han-luo
"""

from __future__ import division
import matplotlib
matplotlib.use('Agg')
from scipy import sparse
from sklearn import  isotonic
from sklearn.decomposition import PCA
from collections import OrderedDict
from scipy.stats import poisson
from statsmodels.sandbox.stats.multicomp import multipletests
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from itertools import  islice
import os, math, random, ghmm, bisect
import cooler, copy
import numpy as np


# A Customized  Chromotin Structure Finding Class
class StructureFind(object):
    """
    This is  a module class for Chromosome 3D Structure Analysis based on Hi-C Data.
    The basis class can fix the haplotype and non-haplotype structures.
    We can use this class to Get :
    
                                  Compartments
                                  TADs.
                                  Chromatin Loops.
        
    This class includes many self-defined methods. But if U do not want fussy details,
    Only some key Functions API are needed.
    
    Input Data
    ----------
    
    This module follows the Matrix building by HiCHap sub-command (matrix).
    The results of Interaction Matrix will be saved as cooler format.
        
    
    
    ---------------------------------------------------------------------------------------------
    
    More API information can be found in methods.
    
    code::block 
    
    	>>> from HiCHap.StructureFind import StructureFind
	
	>>> #============= Compartment==============
	>>> ## For traditional Hi-C
	>>> GM_T_PC = StructureFind(cooler_fil = 'Merged_Traditional_Multi.cool', Res = 500000, Allelic = False)
	>>> GM_T_PC.run_Compartment(OutPath = 'Traditonal_PC', plot = True, MS = 'IF', SA = False)
	 
	>>> ## For haplotype-resolved Hi-C 
	>>> GM_M_PC = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 500000, Allelic = 'Maternal')
	>>> GM_M_PC.run_Compartment(OutPath = 'Maternal_PC', plot = True, MS = 'IF', SA = False)
	
	>>> GM_P_PC = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 500000, Allelic = 'Paternal')
	>>> GM_P_PC.run_Compartment(OutPath = 'Paternal_PC', plot = True, MS = 'IF', SA = False)

	
	>>> #============= TADs calling=============
	>>> ## For traditional Hi-C 
	>>> GM_tads_T = StructureFind(cooler_fil = 'Merged_Traditional_Multi.cool', Res = 40000, Allelic = False)
	>>> GM_tads_T.run_TADs(OutPath = 'Traditional_TADs', plot = True)
	>>>
	>>> ## For haplotype-resolved Hi-C
	>>> GM_tads_M = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic = 'Maternal')
	>>> GM_tads_M.run_TADs(OutPath = 'Maternal_TADs', plot = True)

	>>> GM_tads_P = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic = 'Paternal')
	>>> GM_tads_P.run_TADs(OutPath = 'Paternal_TADs', plot = True)
	

	>>> #============= Loops calling=============
	>>> ## For traditonal Hi-C
	>>> GM_Loop_T = StructureFind(cooler_fil = 'Merged_Traditional_Multi.cool', Res = 40000, Allelic = False)
	>>> GM_Loop_T.run_Loops(OutPath = 'Traditional_Loops', plot = True)
	
	>>> ## For haplotype-resolved Hi-C
	>>> GM_Loop_M = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic = 'Maternal')
	>>> GM_Loop_M.run_Loops(OutPath = 'Maternal_Loops', plot = True)
	
	>>> GM_Loop_P = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic = 'Paternal')
	>>> GM_Loop_P.run_Loops(OutPath = 'Paternal_Loops', plot = True)
	
 
    
    """
    def __init__(self, cooler_fil, Res, Allelic, GapFile = None, 
                 Loop_ratio = 0.6, Loop_strength = 16):
        
        #------------Initialization the HiC Matrix Dictionary----------------
        self.cooler_fil = '{}::{}'.format(cooler_fil, Res)
        self.Res = Res
        self.Allelic = Allelic
        self.Gap_file = GapFile
        self.ratio = Loop_ratio
        self.LoopStrength = Loop_strength

#===============================Public Functions====================================


    def Sig_update(self, sigs):
        """
            Up date the Sigs to Plot.
            
        Parameter
        ---------
        sigs : array
            sigs array.
        """
        New_Sig = []
        New_Index = []
        
        for index in range(len(sigs) - 1):
            if sigs[index] * sigs[index + 1] < 0:
                New_Sig.append(sigs[index])
                New_Sig.append(0)
                New_Index.append(index)
                New_Index.append(index + 0.5)
            else:
                New_Sig.append(sigs[index])
                New_Index.append(index)
        
        return np.array(New_Index), np.array(New_Sig)
    
    
    def Cmap_Setting(self, types = 2, **kwargs):
        """
            Color bar setting. hexadecimal notation(#16) must be input.
        
        Parameters
        ----------
        start_Color : str
            start_Color (default : #FFFFFF)        
        
        end_Color : str
            end_Color (default : #CD0000)
        """
        
        start_Color = kwargs.get('start_Color','#FFFFFF')
        middle_Color = kwargs.get('middle_Color','#FFFFFF')
        end_Color = kwargs.get('end_Color','#CD0000')
        if types == 2:
            return LinearSegmentedColormap.from_list('interactions',[start_Color,end_Color])
        elif types == 3:
            return LinearSegmentedColormap.from_list('interactions',[start_Color,middle_Color,end_Color])
            


    def properU(self, pos):
        """
          Express a genomic position in a proper unit (KB, MB, or both).
          
        """
        i_part = int(pos) // 1000000    # Integer Part
        d_part = (int(pos) % 1000000) // 1000   # Decimal Part
        
        if (i_part > 0) and (d_part > 0):
            return ''.join([str(i_part), 'M', str(d_part), 'K'])
        elif (i_part == 0):
            return ''.join([str(d_part), 'K'])
        else:
            return ''.join([str(i_part), 'M'])
            
            
    def caxi_H(self,ax):
        """
            Axis Control for HeatMaps
            
        """
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params(axis = 'both', labelsize = 10, length = 5, pad = 7)
    
    def fig_settings(self):
        """
            About figure.
        """
        self.figure_Dict = {}
        self.figure_Dict['size'] = (10,9)
        self.figure_Dict['width'] = 0.618
        self.figure_Dict['Left'] = (1 - self.figure_Dict['width']) / 2
        self.figure_Dict['HB'] = 0.1
        self.figure_Dict['HH'] = self.figure_Dict['width'] * self.figure_Dict['size'][0] / self.figure_Dict['size'][1]
        self.figure_Dict['SB'] = self.figure_Dict['HB'] + self.figure_Dict['HH']
        
    
#============================Compartment Analysis modules===========================


    
    def Distance_Decay(self, M, G_array):
        """
            Get the distance decay line.
            
        Parameter
        ---------
        M : 2D HiC Matrix of Intra-Chromosome
        
        G_array : Gap if not available, Get it.
        
        """
        
        size = M.shape[0]
        bin_arange = np.arange(size)        
        
        if G_array == None:
            Gap_ratio = 0.05
            nonzero_mask = (M != 0).sum(axis = 0) / float(size)
            gap_mask = np.where(nonzero_mask <= Gap_ratio,True,False)
            G_array = bin_arange[gap_mask]
            NG_array = bin_arange[~gap_mask]
        else:
            NG_array = np.array([item for item in bin_arange if item not in G_array])
            
        Matrix_mask = []
        dtype = [('bin1',np.int),('bin2',np.int),('IF',np.float)]
        bin1, bin2 = np.nonzero(M)
        IF = M[bin1,bin2]
        
        Sparse_M = np.zeros((len(IF),),dtype = dtype)
        Sparse_M['bin1'] = bin1
        Sparse_M['bin2'] = bin2
        Sparse_M['IF'] = IF
        
        for item in G_array:
            Matrix_mask.extend(np.where(Sparse_M['bin2'] == item)[0])
        
        Matrix_mask_set = np.array(list(set(Matrix_mask)))
        Mask_True = np.ones((len(Sparse_M),),dtype = bool)
        
        if Matrix_mask_set.shape[0] >= 1:
            Mask_True[Matrix_mask_set] = False
        
        Sparse_M_Less = Sparse_M[Mask_True]
        
        weight = np.array(Sparse_M_Less['IF'])
        weight = np.hstack((weight,np.array([0])))

        distance = np.array(Sparse_M_Less['bin2'] - Sparse_M_Less['bin1'], dtype = np.int)
        distance_abs = np.array([abs(x) for x in distance], dtype = np.int)
        distance_size_abs = np.hstack((distance_abs,np.array([size],dtype = np.int)))
        distance_bin = np.bincount(distance_size_abs,weight)
        
        #-------------------------Remove the Gap effect-----------------------------
        for i in range(size):
            if i == 0:
                gap_num_start = 0
                gap_num_end = sum((0 <= G_array) & (G_array <= size - 1))
                gap_num = gap_num_start + gap_num_end
                bin_num = float(size - i) - gap_num
            else:
                gap_num_start = sum((0 <= G_array) & (G_array <= size - 1 - i))
                gap_num_end = sum((i <= G_array) & (G_array <= size - 1))
                gap_num = gap_num_start + gap_num_end
                bin_num = float(size - i) * 2 - gap_num
            
            if bin_num > 0:
                distance_bin[i] = float(distance_bin[i] / bin_num)
        distance_bin = distance_bin[: size]
    
        return distance_bin, G_array, NG_array
        
    
    def Sliding_Approach(self, M, decline, window):
        """
            Sliding Approach by Ren Bin.  may be a little difference

        """
        step = window // self.Res // 2
        
        N = M.shape[0]
        
        OE_matrix_bigger = np.zeros(M.shape)
        
        for i in range(N):
            for j in range(N):
                if i < step or j < step:
                    OE_matrix_bigger[i][j] = M[i][j] / decline[abs(i-j)]
                elif i > N - step - 1 or j > N - step - 1:
                    OE_matrix_bigger[i][j] = M[i][j] / decline[abs(i-j)]
                else:
                    O_Sum = M[i-step:i+step+1, j-step:j+step+1].sum()
                    E_Sum = 3 * decline[abs(i-j)] + 2 * decline[abs(i-j-1)] + \
                            2 *decline[abs(i-j+1)] + decline[abs(i-j-2)] + \
                            decline[abs(i-j+2)]
                    
                    OE_matrix_bigger[i][j] = O_Sum / E_Sum
        
        return OE_matrix_bigger
                    
        
    def Get_PCA(self, distance_bin, M, NG_array, SA = False):
        """
            Based on Distance_Decay 
            Get the:
                Correlation Matrix
                Obersevered / Expected Matrix (OE Matrix)
                PCA component
        
        Parameters
        ----------
        distance_bin : distance_bin created by Distance_Decay
        
            M        : 2D HiC Matrix of Intra-Chromosome
            
        NG_array     : NG_array created by Distance_Decay
        
        """
        
        decline = distance_bin
        decline[decline == 0] = decline[np.nonzero(decline)].min()
        
        OE_matrix_bigger = np.zeros(M.shape)
        
        if SA == False:
            for i in range(M.shape[0]):
                for j in range(M.shape[1]):
                    if M[i][j] != 0:
                        OE_matrix_bigger[i][j] = M[i][j] / decline[abs(i-j)]
        else:
            OE_matrix_bigger = self.Sliding_Approach(M, decline, window=600000)
        
        OE_matrix = OE_matrix_bigger[:,NG_array]
        
        Cor_matrix = np.corrcoef(OE_matrix,rowvar = False)
        Cor_matrix[np.isnan(Cor_matrix)] = 0
        Cor_matrix[np.isinf(Cor_matrix)] = 1
        pca = PCA(n_components = 3)
        pca.fit(Cor_matrix)
        pca_components=np.vstack(pca.components_) # pca_components=np.vstack((pca.components_[0],pca.components_[1],pca.components_[2]))
        
        return pca_components, Cor_matrix, OE_matrix
        
    
    def Select_PC(self, Cor_matrix, pca_components):
        """
            Based On Get_PCA
            Select the PC 
        
        Parameters
        ----------
        Cor_matrix : Cor_matrix created by Get_PCA
        
        OE_matrix  : OE_matrix created by Get_PCA
        
        """
        select_k = 0
        corr_coef = 0
        for i in range(pca_components.shape[0]):
            tmp_coef = np.array([np.corrcoef(pca_components[i], item)[0,1] for item in Cor_matrix])
            tmp_coef[np.isnan(tmp_coef)] = 0
            tmp_coef[np.isinf(tmp_coef)] = 1
            if abs(tmp_coef).sum() > corr_coef:
                corr_coef = abs(tmp_coef).sum()
                select_k = i
                direction = 1
                if tmp_coef.sum() < 0:
                    direction = -1
        
        pc_selected = pca_components[select_k] * direction
        
        return pc_selected
    
    def Select_PC_new(self, Cor_M, OE_M, pca):
        def means_minus(matrix, pc, eps=1e-5):
            lenth = len(pc)
            locis = np.arange(lenth)
            mask_a = pc > 0 
            mask_b = pc < 0
            locis_a = locis[mask_a]
            locis_b = locis[mask_b]
            if locis_a.shape[0] == 0 or locis_b.shape[0] == 0:
                return 0
            size_a = locis_a.max() - locis_a.min()
            size_b = locis_b.max() - locis_b.min()
            lens = max(locis_a.max(),locis_b.max()) - min(locis_a.min(),locis_b.min())
            matrix_a = matrix[mask_a][:,mask_a]
            matrix_b = matrix[mask_b][:,mask_b]
            matrix_ab = matrix[mask_a][:,mask_b]
            mask_a1 = matrix_a > -1
            mask_a2 = matrix_a < 1 - eps
            mask_b1 = matrix_b > -1
            mask_b2 = matrix_b < 1 - eps
            mask_ab1 = matrix_ab > -1
            mask_ab2 = matrix_ab < 1
            value_a = matrix_a[mask_a1 * mask_a2]
            value_b = matrix_b[mask_b1 * mask_b2]
            value_ab = matrix_ab[mask_ab1 * mask_ab2]
            value_same = np.hstack((value_a,value_b))
            if value_ab.shape[0] == 0 or value_ab.mean() == 0 or value_ab.mean() == -1 or size_a <= lens/2 or size_b <= lens/2:
               return 0
            return value_same.mean() - value_ab.mean()
        def select_ab(oe, compartment):
                matrix = oe
                compt = compartment
                mask_a = compt > 0
                mask_b = compt < 0
                matrix_a = matrix[mask_a][:,mask_a]
                matrix_b = matrix[mask_b][:,mask_b]
                values_a = sparse.coo_matrix(matrix_a).data.mean()
                values_b = sparse.coo_matrix(matrix_b).data.mean()
                if values_b > values_a :
                    compartment = compartment.copy() * -1
                return compartment
        nums = 0
        values = 0
        for i in range(len(pca)):
            minus = means_minus(Cor_M, pca[i])
            if minus > values:
                values = minus
                nums = i    
        pc_selected = select_ab(OE_M, pca[nums])
        return pc_selected 
    
    
    def Loading_Tranditional_PC(self, fil):
        """
        """
        PC_Dic = {}
        with open(fil,'r') as f:
            for line in f:
                line = line.strip().split()
                if line[0] not in PC_Dic.keys():
                    PC_Dic[line[0]] = []
                    PC_Dic[line[0]].append(line[-1])
                else:
                    PC_Dic[line[0]].append(line[-1])
        
        for chro, v in PC_Dic.items():
            v = np.array(v, dtype = float)
            PC_Dic[chro] = v
        
        return PC_Dic
    
    
    def Select_Allelic_PC(self, pca_components, Tranditional_PC, eps = 0.7):
        """
            Select the principal component by  supervised manner in Allelic matrix
            
        """
        PCC = []
        for pc in pca_components:
            pcc = abs(np.corrcoef(pc, Tranditional_PC)[0][1])
            PCC.append(pcc)
        if np.max(PCC) < eps:
            print "    PCC too low for this chromosome, check it if possible!"

        index = np.argmax(PCC)
    
        return pca_components[index]
     
     
    def Refill_Gap(self, M1, M2, NonGap, dtype):
        """
            Refill zero value into Gap in OE Matrix and Correlation Matrix.
            
        """
        Refilled_M = np.zeros(M1.shape, dtype = np.float)
        if dtype == 'Cor':
            tmp = np.zeros((M1.shape[0], M2.shape[0]), dtype = np.float)
            for i in range(len(NonGap)):
                insert_pos = NonGap[i]
                tmp[insert_pos] = M2[i]

            tmp = tmp.T
        
            for i in range(len(NonGap)):
                insert_pos = NonGap[i]
                Refilled_M[insert_pos] = tmp[i]

        elif dtype == 'OE':   
            M2 = M2.T
            for i in range(len(NonGap)):
                insert_pos = NonGap[i]
                Refilled_M[insert_pos] = M2[i]
                
                Refilled_M = Refilled_M.T
                
        return Refilled_M
        
    def Compartment(self, SA = False, Tranditional_PC_file = None):
        """
            Compute the Compartment Structure for each Chromosome
        """
        
        print "Compartment Analysis start ..."
        print "Resolution is %s ..." % self.properU(self.Res)
        
        self.cooler = cooler.Cooler(self.cooler_fil)
        
        if self.Allelic == False:
            chroms = self.cooler.chromnames
        elif self.Allelic == 'Maternal':
            chroms = [i for i in self.cooler.chromnames if i.startswith('M')]
        elif self.Allelic == 'Paternal':
            chroms = [i for i in self.cooler.chromnames if i.startswith('P')]
        else:
            raise Exception ('Unkonwn key word %s, Only Maternal, Paternal, False allowed' % self.Allelic)
        
        self.chroms = chroms
        Matrix_Lib = {}
        for chro in chroms:
            Matrix_Lib[chro] = self.cooler.matrix(balance = False).fetch(chro)
        
        self.Matrix_Dict = Matrix_Lib
        self.Cor_Martrix_Dict = {}
        self.OE_Matrix_Dict = {}
        self.Compartment_Dict = {}
        
        if self.Allelic != False:
            Tranditional_PC = self.Loading_Tranditional_PC(Tranditional_PC_file)
        
        self.RawPCA = {}
        for chro in chroms:
            print "Chromosome %s start" % chro
            M = Matrix_Lib[chro]
            self.RawPCA[chro] = []
            distance_bin, Gap, NonGap = self.Distance_Decay(M=M, G_array=None)
            pca, Cor_M, OE_M = self.Get_PCA(distance_bin=distance_bin, M=M, NG_array=NonGap, SA = SA)
            
            if self.Allelic == False:
                pc_select = self.Select_PC_new(Cor_M, OE_M[NonGap], pca) # pc_select = self.Select_PC(Cor_matrix=Cor_M, pca_components=pca)
                Compartment_zeros = np.zeros((M.shape[0],), dtype = np.float)
                Compartment_zeros[NonGap] = pc_select
                self.Compartment_Dict[chro] = Compartment_zeros
                
            else:
                for i in range(len(pca)):
                    tmp = np.zeros((M.shape[0],), dtype = np.float)
                    tmp[NonGap] = pca[i]
                    self.RawPCA[chro].append(tmp.tolist())

                self.RawPCA[chro] = np.array(self.RawPCA[chro])
                Compartment_zeros = np.zeros((M.shape[0],), dtype = np.float)
                tran_chro = chro[1:]
                pc_select = self.Select_Allelic_PC(self.RawPCA[chro], Tranditional_PC[tran_chro])
                Compartment_zeros[NonGap] = pc_select[NonGap]
                self.Compartment_Dict[chro] = Compartment_zeros
                
            OE_M = self.Refill_Gap(M, OE_M, NonGap, dtype = 'OE')
            Cor_M = self.Refill_Gap(M, Cor_M, NonGap, dtype = 'Cor')

            self.Cor_Martrix_Dict[chro] = Cor_M
            self.OE_Matrix_Dict[chro] = OE_M
   

    def OutPut_PC_To_txt(self, out):
        """
            Output the PC component to a text file
        
        Parameter
        ---------
        out : str
            out file
        """
        f = open(out, 'w')
        if self.Allelic == False:
            for chro, pc in self.Compartment_Dict.items():
                for value in pc:
                    f.writelines(chro+'\t'+str(value)+'\n')
        else:
            for chro, pc in self.Compartment_Dict.items():
                for value in pc:
                    f.writelines(chro[1:]+'\t'+str(value)+'\n')
        
        f.close()
    
    
    def Plot_Compartment(self, out, MS = 'IF'):
        """
            Output the Compartment Structure to a  PDF
        
        Parameter
        --------
        out : str
            out PDF file.
        
        Matrix :  Matrix Str type
        
                'IF' : Interaction Matrix.
                'OE' : Observed / Expected Matrix.
                'Cor' : Correlation Matrix.
                
        """
        pp = PdfPages(out)
        self.fig_settings()
        if MS == 'IF':
            Matrix_Dict = self.Matrix_Dict
            cmap = self.Cmap_Setting()
        elif MS == 'OE':
            Matrix_Dict = self.OE_Matrix_Dict
            cmap = self.Cmap_Setting(types=3,start_Color='#0000FF')
        elif MS == 'Cor':
            Matrix_Dict = self.Cor_Martrix_Dict
            cmap = self.Cmap_Setting(types=3,start_Color='#0000FF')
        else:
            raise Exception("Unknown Matrix String %s, Only IF, OE, Cor strings allowed." % MS)
        
        for chro in self.chroms:
            Matrix = Matrix_Dict[chro]
            Sigs = self.Compartment_Dict[chro]
            N = Matrix.shape[0]
            
            if MS == 'IF':
                nonzero = Matrix[np.nonzero(Matrix)]
                vmax = np.percentile(nonzero,95)
                vmin = 0
            elif MS == 'OE':
                nonzero = Matrix[np.nonzero(Matrix)]
                vmax = np.percentile(nonzero, 90)
                vmin = 2 - vmax
            elif MS == 'Cor':
                nonzero = Matrix[np.nonzero(Matrix)]
                vmax = np.percentile(nonzero, 90)
                vmin = -vmax
            
            fig = plt.figure(figsize = self.figure_Dict['size'])
            ax = fig.add_axes([self.figure_Dict['Left'], self.figure_Dict['HB'],
                               self.figure_Dict['width'], self.figure_Dict['HH']])
                               
            sc = ax.imshow(Matrix, cmap = cmap, aspect = 'auto', interpolation = 'none',
                           extent = (0, N, 0, N), vmin = vmin, vmax = vmax, origin = 'lower')
            
            ticks = list(np.linspace(0, N, 5).astype(int))
            pos = [t * self.Res for t in ticks]
            labels = [self.properU(p) for p in pos]
            
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)
            if self.Allelic == False:
                ax.set_xlabel('Chr'+chro, size = 14)
            else:
                ax.set_xlabel('Chr'+chro[1:], size = 14)
            
            
            
            ax = fig.add_axes([self.figure_Dict['Left']+self.figure_Dict['width']+0.02,
                               self.figure_Dict['HB'],0.01,self.figure_Dict['HH']])
            
            fig.colorbar(sc, cax = ax)
            
            index, sigs = self.Sig_update(Sigs)
            
            ax2 = fig.add_axes([self.figure_Dict['Left'],self.figure_Dict['SB'],
                                self.figure_Dict['width'],self.figure_Dict['HB']])

            for spine in ['right', 'top', 'left']:
                ax2.spines[spine].set_visible(False)
                
            ax2.fill_between(index, sigs, where = sigs <= 0, color = '#7093DB')
            ax2.fill_between(index, sigs, where = sigs >= 0, color = '#E47833')
            ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
            
            ax2.set_xlim(0,len(Sigs))
            ax2.set_ylabel('PC', size = 12)
            
            pp.savefig(fig)
            plt.close(fig)
        
        pp.close()                    
        
        
    def run_Compartment(self, OutPath, plot = True, MS = 'IF', 
                        SA = False, Tranditional_PC_file = None):
        """
            Main function to Get Compartment
        
        Parameters
        ----------
        OutPath : str
            Out Put Path
        """
        OutPath.rstrip('/')
        if os.path.exists(OutPath):
            pass
        else:
            os.mkdir(OutPath)
        
        prefix = os.path.split(OutPath)[-1]
        res = self.properU(self.Res)
        
        fil = os.path.join(OutPath, prefix+'_Compartment_'+res+'.txt')
        pdf = os.path.join(OutPath, prefix+'_Compartment_'+MS+'_'+res+'.pdf')
        
        self.Compartment(SA = SA, Tranditional_PC_file=Tranditional_PC_file)
        self.OutPut_PC_To_txt(fil)
        if plot:
            self.Plot_Compartment(pdf, MS)
        

#==============================TAD Analysis module==================================



    def TAD_parameter_init(self, minTAD, maxTAD, state_num, window, test_type):
        """
            Initial the prior parameters
        """
        
        self.minTAD = minTAD
        self.maxTAD = maxTAD
        self.state_num = state_num
        self.window = window
        self.test_type = test_type


    def Get_Gap(self, M):
        """
            Get the Gap for TAD calling.
        
        Parameters
        ----------
        M : array
            Matrix Data
        
        minTAD : int
            the min TAD size
            
        """
        ratio = 0
        N = M.shape[0]
        col_distribution = (np.array(M) != 0).sum(axis = 0) / float(N)
        gap_flag = np.where(col_distribution < ratio, True, False)
        gap_arange = np.arange(N)
        gap_array_s = np.zeros((N,),dtype = np.bool)
        gap_array_e = np.zeros((N,),dtype = np.bool)
        
        local_bin = int(self.minTAD / self.Res)
        
        t = 2 * local_bin * 0.8
        for i in xrange(N):
            col_num_e = sum(M[i - local_bin:i+local_bin, i] != 0) < t if ((local_bin<=i)&(i<=N-1-local_bin)) else True
            gap_array_e[i] = col_num_e
        
        gap_array = gap_arange[gap_flag | (gap_array_s | gap_array_e)]
        
        return gap_array
    
    def Gap_Filter(self, Gap, M):
        """
            Gap Filtering.
        
        Prameters
        ---------
        Gap : creadted by TAD_Gap
        
        M : Matrix Data
        
        """
        if Gap.shape[0] <= 1:
            Gap_filtered = []
            return Gap_filtered
        
        Gap_Dict = OrderedDict()
        count_start = Gap[0]
        count_end = Gap[0]
        Gap_Len = Gap.shape[0]
        for i in xrange(1, Gap_Len):
            if (Gap[i] - Gap[i - 1] == 1)&(Gap_Len - 1 == i):
                count_end = Gap[i] + 1
                Gap_Dict[(count_start,count_end)] = count_end - count_start
            elif Gap[i] - Gap[i - 1] == 1:
                count_end = Gap[i] + 1
            else:
                Gap_Dict[(count_start,count_end)] = count_end - count_start
                count_start = Gap[i]
                count_end = Gap[i] + 1
        Key_list = Gap_Dict.keys()
        Key_list.sort()
        Gap_Len = [Gap_Dict[item] for item in Key_list]
        
        Gap_mean = np.mean(Gap_Len)
        
        Gap_Dict_filtered = {item:Gap_Dict[item] for item in Key_list if Gap_Dict[item] >= min([10,Gap_mean])}
        
        Key_list = Gap_Dict_filtered.keys()
        Key_list.sort()
        Gap_filtered = []

        for key in Key_list:
            Gap_filtered.extend(range(key[0],key[1]))
        
        if 0 not in Gap_filtered:
            Gap_filtered.insert(0,0)
        if M.shape[0] - 1 not in  Gap_filtered:
            Gap_filtered.append(M.shape[0] - 1)
        
        return Gap_filtered

    def Get_DI(self, M, Gap, window_bin):
        """
            calculate the DI value.
        """
        DI = []
        N = M.shape[0]
        for j_ind in xrange(N):
            window_binj=window_bin[j_ind]
            if j_ind in Gap:
                DI.append(0)
            elif(j_ind<window_binj)|(j_ind>N-window_binj-1):
                DI.append(0)
            elif(window_binj<=j_ind)&(j_ind<=N-window_binj-1):
                up=M[j_ind-window_binj:j_ind,j_ind]
                up=up[::-1]
                down=M[j_ind+1:j_ind+window_binj+1,j_ind]
                bias=0
                if self.test_type=='ttest':
                    up_mean=up.mean();down_mean=down.mean()
                    up_denominator=np.sum((up-up_mean)**2/(up.size*(up.size-1)))
                    down_denominator=np.sum((down-down_mean)**2/
                                            (down.size*(down.size-1)))
                    denominator_sum=np.sqrt(up_denominator+down_denominator)
                    if(denominator_sum!=0):
                        bias=(down_mean-up_mean)/denominator_sum
                elif self.test_type=='chitest':
                    up_sum=up.sum();down_sum=down.sum()
                    expect_con=float(up_sum+down_sum)/2.0
                    if((up_sum!=down_sum)&(expect_con!=0)):
                        bias=(float(down_sum-up_sum)/(abs(down_sum-up_sum))*
                              ((up_sum-expect_con)**2/expect_con+
                              (down_sum-expect_con)**2/expect_con))
                DI.append(bias)
        DI = np.array(DI)
        
        return DI
    
    
    def Data_preprocess(self):
        """
            Data preprocess for ghmm model.
        
        """
        Matrix_Dict = {}
        self.cooler = cooler.Cooler(self.cooler_fil)

        if self.Allelic == False:
            chroms = self.cooler.chromnames
            for chro in chroms:
                M = self.cooler.matrix(balance = True).fetch(chro)
                M = np.nan_to_num(M)
                Matrix_Dict[chro] = M
        elif self.Allelic == 'Maternal':
            chroms = [i for i in self.cooler.chromnames if i.startswith('M')]
            for chro in chroms:
                M = self.cooler.matrix(balance = False).fetch(chro)
                Matrix_Dict[chro] = M
        elif self.Allelic == 'Paternal':
            chroms = [i for i in self.cooler.chromnames if i.startswith('P')]
            for chro in chroms:
                M = self.cooler.matrix(balance = False).fetch(chro)
                Matrix_Dict[chro] = M
        else:
            raise Exception ('Unkonwn key word %s, Only Maternal, Paternal, False allowed' % self.Allelic)
            
        window_bin = int(self.window / self.Res)
        width = 7
        Gap_all = {}
        DI_dict = {}
        DI_all_train = {}
        
        for chro in Matrix_Dict.keys():
            matrix = Matrix_Dict[chro]
            Gap = self.Get_Gap(matrix)
            
            # check the first and last bin
            tmp = list(Gap)
            if 0 not in tmp:
                tmp.insert(0,0)
            if matrix.shape[0] - 1 not in tmp:
                tmp.append(matrix.shape[0] - 1)
            
            Gap = np.array(tmp)
            
            Gap_desity_t = float(Gap.size) / matrix.shape[0] / 2.0
            window_bins = np.ones(matrix.shape[0], dtype = np.int) * window_bin
            
            DI_sub = self.Get_DI(matrix, Gap, window_bins)
            Gap_all[chro] = Gap
            Gap_fitered = self.Gap_Filter(Gap, matrix)
            DI_dict[chro] = DI_sub
            DI_train_dict = {}
            
            for i in xrange(1, len(Gap_fitered)):
                #filter the short DI seq
                if Gap_fitered[i] - Gap_fitered[i - 1] <= width:
                    continue
                
                #filter the DI seq with the Gap density > Gap_desity_t
                if((sum((Gap_fitered[i-1] < Gap) & (Gap < Gap_fitered[i])) / 
                    float(Gap_fitered[i] - Gap_fitered[i-1]-1)) > Gap_desity_t):
                        continue
                
                DI_train_dict[(Gap_fitered[i - 1] + 1, Gap_fitered[i])] = DI_sub[Gap_fitered[i - 1] + 1:Gap_fitered[i]]
            
            DI_all_train[chro] = DI_train_dict
        
        self.DI_all_train = DI_all_train
        self.DI_dict = DI_dict
        self.Gap_all = Gap_all
        self.Matrix_Dict = Matrix_Dict
        self.chroms = chroms
            
        
    def init_parameter_state3(self):
        """
        initial the parameter of ghmm
        3 states
        A:initial states transition probability distribution
        0 -- downstream bias
        1 -- no bias
        2 -- upstream bias
        don't select: '10' -- no bias --> downstream bias: start
        don't select: '21' -- upstream bias --> no bias: end
        '20' -- upstream bias --> downstream bias: start & end
        A means:
            1 downstream more likely change to downstream,and less likely change to upstream
            2 no bias more likely change to no bias,and less likely change to upstream
            3 upstream bias more likely change to upstream and less likely change to no bias            
        """
        
        A = [[0.85, 0.15, 0.00],
             [0.05, 0.80, 0.15],
             [0.19, 0.01, 0.80]]
        #pi: initial states probability distribution
        pi = [0.40, 0.30, 0.30]
        #numdists:Three-distribution Gaussian Mixtures
        numdists = 3
        #W: weight of every gaussian distribution
        W = 1.0 / numdists
        #Var: variance
        num = 6.0
        var = num / (numdists - 1)
        #means : average
        means = [[],[],[]]
        for i in range(numdists):
            means[0].append((i+1)*var)
            means[1].append((i-1)*var)
            means[2].append((i-2)*var)
        B=[[means[i],[var for j in range(numdists)],[W for k in 
              range(numdists)]] for i in range(len(pi))]

        return A, B, pi

    def init_parameter_state5(self):
        """
        initial the parameter of ghmm
        5 states
        A:initial states transition probability distribution
        0 -- start
        1 -- downstream bias
        2 -- no bias
        3 -- upstream bias
        4 -- end
        '20' -- no bias --> start
        '40' -- end --> start
        '42' -- end --> no bias
        start can only change to downstream bias
        downstream bias can only change to downstream bias and no bias
        no bias can only change to start/no bias/upstream bias
        upstream bias can only change to upstream bias end
        end can only change to start/no bias；end can only be changed from upstream
        """
        
        A=[[0.00,1.00,0.00,0.00,0.00],
           [0.00,0.50,0.50,0.00,0.00],
           [0.33,0.00,0.34,0.33,0.00],
           [0.00,0.00,0.00,0.50,0.50],
           [0.50,0.00,0.50,0.00,0.00]]        
        #pi: initial states probability distribution 
        pi=[0.05,0.3,0.3,0.3,0.05]
        #numdists:Three-distribution Gaussian Mixtures
        numdists=3
        #W:weight of every gaussian distribution
        W=1./numdists
        #var:variance
        num=6.0
        var=num/(numdists-1)
        #means: average
        means=[[],[],[],[],[]]                 
        for i in range(numdists):
            means[0].append((i+1) * var)
            means[1].append((i) * var)
            means[2].append((i-1) * var)
            means[3].append((i-2) * var)
            means[4].append((i-3) * var)
        B=[[means[i],[var for j in range(numdists)],[W for k in 
              range(numdists)]] for i in range(len(pi))]
        
        return A, B, pi
    
    
    def init_parameter_state6(self):
        """
        initial the parameter of ghmm
        6 states
        A:initial states transition probability distribution
        0 -- start
        1 -- downstream bias
        2 -- no bias
        3 -- upstream bias
        4 -- end
        5 -- gap
        """
        
        A = [[0.00,1.00,0.00,0.00,0.00,0.00],
             [0.00,0.75,0.20,0.00,0.00,0.05],
             [0.00,0.00,0.60,0.35,0.00,0.05],
             [0.00,0.00,0.00,0.93,0.02,0.05],
             [0.20,0.60,0.20,0.00,0.00,0.00],
             [0.00,0.22,0.06,0.22,0.00,0.50]]
        
        #pi : initial states probability distribution
        pi = [0.01,0.29,0.20,0.10,0.05,0.35]
        #numdists: Three-distribution Gaussian Mixtures
        numdists = 3
        #W:weight of every gaussian distribution
        W = 1.0 / numdists
        #var:variance
        num = 4.2
        var = num / (numdists - 1)
        #mean: average
        means = [[],[],[],[],[],[]]
        for i in range(numdists):
            means[5].append(0)
            means[4].append((i+1)*var)
            means[3].append((i)*var)
            means[2].append((i-1)*var)
            means[1].append((i-2)*var)
            means[0].append((i-3)*var)        
        B=[[means[i],[var for j in range(numdists)],[W for k in 
              range(numdists)]] for i in range(len(pi))]
        
        B[5][1]=[0.0001,0.0001,0.0001]     
        
        return A, B, pi
        
    
    def updateParameter(self, DI_train, A, B, pi):
        """
            Update the parameter
        """
        F = ghmm.Float()
        chromosome = DI_train.keys()
        random.shuffle(chromosome)
        seqs = []
        for chrom in chromosome:
            domain_list = DI_train[chrom].keys()
            domain_list.sort()
            random.shuffle(domain_list)
            temp = []
            for item in domain_list:
                temp.append(DI_train[chrom][item])
            seqs.extend(temp)
        trainset = ghmm.SequenceSet(F, seqs)
        
        model = ghmm.HMMFromMatrices(F, ghmm.GaussianMixtureDistribution(F), A, B, pi)
        model.baumWelch(trainset,nrSteps=10**8)
        
        state_num = len(B)
        gauss_num = len(B[0][0])
        AA = np.zeros((state_num, state_num))
        BB = np.zeros((state_num, 3, gauss_num))
        pipi = np.zeros(state_num)
        
        for i in xrange(state_num):
            for j in xrange(state_num):
                AA[i, j] = model.getTransition(i, j)
            for j in xrange(gauss_num):
                temp = model.getEmission(i, j)
                for k in range(3):
                    BB[i, k, j] = temp[k]
            pipi[i] = model.getInitial(i)
        
        return AA, BB, pi, model
    
    
    def modelTrain(self):
        """
            train the ghmm model.
        """
        self.Data_preprocess()
        
        if self.state_num == 3:
            A, B, pi = self.init_parameter_state3()
        elif self.state_num == 5:
            A, B, pi = self.init_parameter_state5()
        elif self.state_num == 6:
            A, B, pi = self.init_parameter_state6()
        else:
            raise Exception("Error state number, Only 3, 5, 6 will be accepeted")
        
        AA, BB, pi, model = self.updateParameter(self.DI_all_train, A, B, pi)        
        AA, BB, pi, model = self.updateParameter(self.DI_all_train, AA, BB, pi)
        AA, BB, pi, model = self.updateParameter(self.DI_all_train, AA, BB, pi)
        
        self.model = model
        
    
    def viterbipath(self, DI_dict):
        """
        """
        F = ghmm.Float()
        paths_dict = {}
        
        for d in DI_dict.keys():
            tmp = ghmm.EmissionSequence(F, list(DI_dict[d]))
            paths_dict[d] = self.model.viterbi(tmp)
        
        return paths_dict
    
    
    def BoundaryMask(self, origin_range, mask_str):
        """
            The Boundary Selecting Mask.
        """
        for item in mask_str:
            item_len = len(item[0])
            start_end_flag = (item[1] == item[2])
            
            for i in xrange(origin_range.shape[0] - item_len + 1):
                temp = ''
                for tempi in origin_range['raw_state'][i:i+item_len]:
                    temp += tempi
                if temp == item[0]:
                    if start_end_flag:
                        origin_range['state'][i+item[1]] = 'both'
                    else:
                        if(item[1] >= 0):
                            if origin_range['state'][i+item[1]] == 'end':
                                origin_range['state'][i+item[1]]='both'
                            else:
                                origin_range['state'][i+item[1]]='start'
                        if(item[2] >= 0):
                            if(origin_range['state'][i+item[2]]=='start'):
                                origin_range['state'][i+item[2]]='both'
                            else:
                                origin_range['state'][i+item[2]]='end'
        
        mask_boundary = (origin_range['state'] != 'none')
        
        return mask_boundary
    
    
    def BoundaryCall(self, paths_sub, Gap_sub, DI_len_sub):
        """
            Boundary calling.
        """
        index_range = np.arange(DI_len_sub)
        dtype = [('boundary',np.int),('state','>S5'),('rely',np.float),
                 ('raw_state','>S1')]
        #origin_range:
        #   boundary--site index,every bin corresponding to one chromosome
        #   state--the gap state is '5'
        #   rely--the ghmm rely return
        #   raw_state--the state returned by ghmm
        origin_range = np.zeros((DI_len_sub,), dtype = dtype)
        origin_range['boundary']=index_range
        origin_range['raw_state'] = '5'
        origin_range['state'] = 'none'
        
        for d in paths_sub.keys():
            origin_range['raw_state'][d[0]:d[1]] = paths_sub[d][0]
            origin_range['rely'][d[0]:d[1]] = paths_sub[d][1]
        if self.state_num == 3:
            mask_str = [('220',2,2),('200',1,1),('2221',3,3),('1000',1,1)]
            mask_boundary = self.BoundaryMask(origin_range, mask_str)
        elif self.state_num == 5:
            mask_str = [('40',1,1)]
            mask_boundary = self.BoundaryMask(origin_range, mask_str)
        
        boundary_index = origin_range[mask_boundary]
        boundary_index['boundary'] = boundary_index['boundary'] * self.Res
        
        return boundary_index
        
    
    def modelPredict(self):
        """
             Model Predicting.
        """
        self.modelTrain()
        boundary_index_dict = {}
        
        for chrom in self.DI_all_train.keys():
            DI_dict = self.DI_all_train[chrom]
            paths_sub = self.viterbipath(DI_dict)
            Gap_dict = self.Gap_all[chrom]
            DI_len = len(self.DI_dict[chrom])
            boundary_index = self.BoundaryCall(paths_sub = paths_sub,
                                               Gap_sub= Gap_dict,
                                               DI_len_sub = DI_len)
            
            boundary_index_dict[chrom] = boundary_index
        
        self.boundary_index = boundary_index_dict
    
    
    def Candidate_domains(self):
        """
            Region may be domains contains Gap.
        """
        candidate_type = [('chr','>S2'),('start',np.int),('end',np.int)]
        candidate_domain = {}
        for chrom in self.DI_all_train.keys():
            domain_key = self.DI_all_train[chrom].keys()
            domain_key.sort()
            domain_start = [domain[0] for domain in domain_key]
            domain_end = [domain[1] for domain in domain_key]
            candidate_domain[chrom] = np.zeros((len(domain_start),), dtype = candidate_type)
            
            candidate_domain[chrom]['chr'] = chrom
            candidate_domain[chrom]['start'] = np.array(domain_start) * self.Res
            candidate_domain[chrom]['end'] = np.array(domain_end) * self.Res
            
        self.candidate_domain = candidate_domain
    
    
    def BoundaryFilter(self):
        """
            Filter Boundary
        """
        width = 7
        Boundary_filtered_dic = {}
        
        for chrom in self.boundary_index.keys():
            Gap_sub = self.Gap_all[chrom]
            for boundary_i in xrange(len(self.boundary_index[chrom]['boundary'])):
                boundary = self.boundary_index[chrom]['boundary'][boundary_i]
                boundary_bin = boundary / self.Res           
                if((sum(((boundary_bin-width)<=Gap_sub)&(Gap_sub<=boundary_bin))
                        >=(width-1)/2.0)&
                    (sum(((boundary_bin)<=Gap_sub)&(Gap_sub<=boundary_bin+width))
                        >=(width-1)/2.0)):                
                    self.boundary_index[chrom]['state'][boundary_i]='none'
                elif((sum(((boundary_bin-width)<=Gap_sub)&(Gap_sub<=boundary_bin))
                        >=(width-1)/2.0)&
                    (self.boundary_index[chrom]['state'][boundary_i]!='end')):
                    self.boundary_index[chrom]['state'][boundary_i]='start'
                elif((sum(((boundary_bin-width)<=Gap_sub)&(Gap_sub<=boundary_bin))
                        >=(width-1)/2.0)& 
                    (self.boundary_index[chrom]['state'][boundary_i]=='end')):
                    self.boundary_index[chrom]['state'][boundary_i]='none'
                elif((sum(((boundary_bin)<=Gap_sub)&(Gap_sub<=boundary_bin+width))
                        >=(width-1)/2.0)&
                    (self.boundary_index[chrom]['state'][boundary_i]!='start')):
                   self.boundary_index[chrom]['state'][boundary_i]='end'
                elif((sum(((boundary_bin)<=Gap_sub)&(Gap_sub<=boundary_bin+width))
                        >=(width-1)/2.0)&
                    (self.boundary_index[chrom]['state'][boundary_i]=='start')):
                    self.boundary_index[chrom]['state'][boundary_i]='none'
            mask=self.boundary_index[chrom]['state']!='none'
            Boundary_filtered_dic[chrom]=self.boundary_index[chrom]['boundary'][mask]
        
        self.boundary_filtered = Boundary_filtered_dic
    
    
    def BoundaryToDomain(self):
        """
            Created Domains based on Boundary.
        """
        self.Candidate_domains()
        domain_dict = {}
        domain_type=[('start',np.int),('end',np.int)]
        
        for chrom in self.boundary_index.keys():
            domain_list = [[],[]]
            boundary_index_sub = self.boundary_index[chrom]['boundary']
            state_sub = self.boundary_index[chrom]['state']
            
            for ind in range(len(boundary_index_sub) - 1):
                
                #judge whether the two boundarys be on candidate domain
                start_i = np.where(((self.candidate_domain[chrom]['start']
                                    <=boundary_index_sub[ind])
                                    &(boundary_index_sub[ind]
                                    <=self.candidate_domain[chrom]['end']))==True)[0][0]
                
                end_i = np.where(((self.candidate_domain[chrom]['start']
                                    <=boundary_index_sub[ind + 1])
                                    &(boundary_index_sub[ind + 1]
                                    <=self.candidate_domain[chrom]['end']))==True)[0][0]
                #judge whether the two boundarys can be domain or not
                #the state none、end or start can't be domain
                if((start_i!=end_i)|
                    ((state_sub[ind]=='none')|(state_sub[ind]=='end'))| 
                    ((state_sub[ind+1]=='none')|(state_sub[ind+1]=='start'))):
                        continue
                
                #detect the continuous gap of two neighborhood boundary
                four_gap_num = 0
                three_gap_num = 0
                two_gap_num = 0
                for jnd in xrange(int(boundary_index_sub[ind] / self.Res),
                                      int(boundary_index_sub[ind + 1] / self.Res - 3)):
                    if sum(self.DI_dict[chrom][jnd:jnd + 4] == 0) == 4:
                        four_gap_num += 1
                        break
                    elif sum(self.DI_dict[chrom][jnd:jnd + 3] == 0) == 3:
                        three_gap_num += 1
                        break
                    elif sum(self.DI_dict[chrom][jnd:jnd + 2] == 0) == 2:
                        two_gap_num += 1
                if ((four_gap_num>=1)|(three_gap_num>=2)|(two_gap_num>=3)):
                    continue
                
                #the gap ratio in domain should less than 1/3
                domain_gap = 1.0 / 3.0
                
                if (sum(self.DI_dict[chrom][int(boundary_index_sub[ind] / self.Res):
                    int(boundary_index_sub[ind + 1] / self.Res)] == 0) >
                    ((boundary_index_sub[ind + 1] - boundary_index_sub[ind]) / self.Res 
                                                                            * domain_gap)):
                    continue
                #the dommain should bigger than min TAD size
                if boundary_index_sub[ind + 1] - boundary_index_sub[ind] < self.minTAD:
                    continue
                # the domain should smaller than max TAD size
                if boundary_index_sub[ind + 1] - boundary_index_sub[ind] > self.maxTAD:
                    continue
                
                domain_list[0].append(boundary_index_sub[ind])
                domain_list[1].append(boundary_index_sub[ind + 1])
            
            domain_dict[chrom] = np.zeros((len(domain_list[0]),),dtype = domain_type)
            domain_dict[chrom]['start'] = np.array(domain_list[0])
            domain_dict[chrom]['end'] = np.array(domain_list[1])
        
        self.Domain_dict = domain_dict
    
    
    def Plot_TAD(self, outpdf, length = 4000000):
        """
            Plotting TAD with DI signal 
        """
        pp = PdfPages(outpdf)
        cmap = self.Cmap_Setting()
        self.fig_settings()
        
        for chro in self.chroms:
            Matrix = self.Matrix_Dict[chro]
            Sigs = self.DI_dict[chro]
            TADs = self.Domain_dict[chro]
            N = Matrix.shape[0]
            
            interval = length // self.Res
            
            startHiC = 0
            Len = N // interval
            
            for idx in range(Len):
                endHiC = startHiC + interval
                
                fig = plt.figure(figsize = self.figure_Dict['size'])
                ax = fig.add_axes([self.figure_Dict['Left'], self.figure_Dict['HB'],
                                   self.figure_Dict['width'], self.figure_Dict['HH']])
                
                M = Matrix[startHiC:endHiC,startHiC:endHiC]
                tmp_sigs = Sigs[startHiC:endHiC]
                mask = ((TADs['start'] > startHiC * self.Res) & (TADs['start'] < endHiC * self.Res)) |\
                       ((TADs['end'] > startHiC * self.Res) & (TADs['end'] < endHiC * self.Res))                       
                tmp_tads = TADs[mask]
                
                nonzero = M[np.nonzero(M)]
                if nonzero.size <= 100:
                    plt.close(fig)
                    startHiC = endHiC
                    continue
                
                vmax = np.percentile(nonzero, 95)
                sc = ax.imshow(M, cmap = cmap, aspect = 'auto', interpolation = 'none',
                               extent = (0, interval, 0, interval), vmax = vmax, origin = 'lower')
                
                for tad in tmp_tads:
                    s_t = tad['start'] // self.Res - startHiC
                    s_e = tad['end'] // self.Res - startHiC
                    ax.plot([s_t,s_e],[s_t,s_t], color = 'black', lw = 1.0)
                    ax.plot([s_t,s_e],[s_e,s_e], color = 'black', lw = 1.0)
                    ax.plot([s_t,s_t],[s_t,s_e], color = 'black', lw = 1.0)
                    ax.plot([s_e,s_e],[s_t,s_e], color = 'black', lw = 1.0)
                    
                ax.set_xlim(0,(endHiC - startHiC))
                ax.set_ylim(0,(endHiC - startHiC))
                
                ticks = list(np.linspace(0,interval,5).astype(int))
                pos = [((startHiC + t) * self.Res) for t in ticks]
                labels = [self.properU(p) for p in pos]
                ax.set_xticks(ticks)
                ax.set_xticklabels(labels)
                ax.set_yticks(ticks)
                ax.set_yticklabels(labels)
                if self.Allelic == False:
                    ax.set_xlabel('Chr'+chro, size = 14)
                else:
                    ax.set_xlabel('Chr'+chro[1:], size = 14)
                
                ax = fig.add_axes([self.figure_Dict['Left']+self.figure_Dict['width']+0.02,
                                   self.figure_Dict['HB'],0.01,self.figure_Dict['HH']])
                
                fig.colorbar(sc, cax = ax)
                
                index, sigs = self.Sig_update(tmp_sigs)
                ax2 = fig.add_axes([self.figure_Dict['Left'],self.figure_Dict['SB'],
                                    self.figure_Dict['width'],self.figure_Dict['HB']])
                
                for spine in ['right','top','left']:
                    ax2.spines[spine].set_visible(False)
                
                ax2.fill_between(index, sigs, where = sigs <= 0, color = '#7093DB')
                ax2.fill_between(index, sigs, where = sigs >= 0, color = '#E47833')
                ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                                right = False, labelbottom =False, labeltop = False,
                                labelleft = False, labelright = False)
                ax2.set_xlim(0,len(tmp_sigs))
                ax2.set_ylabel('DI', size = 12)
                
                pp.savefig(fig)
                plt.close(fig)
                startHiC = endHiC
        
        pp.close()
        
        
                
    def run_TADs(self, OutPath, **kwargs):
        """
            Run TAD calling.
        Parameters
        ----------
        OutPath : path
            Out put Path
        
        Kwargs : dict
            
            TAD Calling parameters
            
            minTAD : int
                the minimum TAD size
            
            maxTAD : int
                the maxinum TAD size
            
            state_num : int
                number of state in ghmm
            
            window : int
                window size for DI algorithm
            
            test_type : str
                test type for DI algorithm. (ttest, chitest) 
                
            plot : bool
                whether plot
                
        """
        
        minTAD = kwargs.get('minTAD', 200000)
        maxTAD = kwargs.get('maxTAD', 4000000)
        state_num = kwargs.get('state_num', 3)
        window = kwargs.get('window', 600000)
        test_type = kwargs.get('test_type','ttest')
        PLOT = kwargs.get('plot',True)
        
        
        print " Parameters Initialization ..."
        self.TAD_parameter_init(minTAD = minTAD,
                                maxTAD = maxTAD,
                                state_num = state_num,
                                window = window,
                                test_type = test_type)
        
        print " DI calculting and ghmm model building ..."
        self.modelPredict()
        print " Boundary filtering ..."
        self.BoundaryFilter()
        print " Domain preparing ..."
        self.BoundaryToDomain()
        
        print " Output the results in %s"  % OutPath
        
        if os.path.exists(OutPath):
            pass
        else:
            os.mkdir(OutPath)
        
        res = self.properU(self.Res)
        prefix = os.path.split(OutPath.rstrip('/'))[-1]
        
        #Output DI
        DI_out = open(os.path.join(OutPath, prefix+'_DI_'+res+'.txt'),'w')
        if self.Allelic == False:
            for chro, DI_sub in self.DI_dict.items():
                for v in DI_sub:
                    DI_out.writelines(chro+'\t'+str(v)+'\n')
        elif self.Allelic == 'Maternal':
            for chro, DI_sub in self.DI_dict.items():
                for v in DI_sub:
                    DI_out.writelines(chro[1:]+'\t'+str(v)+'\n')
        elif self.Allelic == 'Paternal':
            for chro, DI_sub in self.DI_dict.items():
                for v in DI_sub:
                    DI_out.writelines(chro[1:]+'\t'+str(v)+'\n')
        else:
            raise Exception ('Unkonwn key word %s, Only Maternal, Paternal, False allowed' % self.Allelic)
    
        DI_out.close()
        
        #Output All Boundary.
        All_B_out = open(os.path.join(OutPath, prefix+'_All_Boundary_'+res+'.txt'),'w')
        if self.Allelic == False:            
            for chro, B_sub in self.boundary_index.items():
                for b in B_sub['boundary']:
                    All_B_out.writelines(chro+'\t'+str(b)+'\n')
        else: 
            for chro, B_sub in self.boundary_index.items():
                for b in B_sub['boundary']:
                    All_B_out.writelines(chro[1:]+'\t'+str(b)+'\n')
            
        All_B_out.close()
        
        #Output filtered Boundary.
        filtered_B_out = open(os.path.join(OutPath, prefix+'_Filtered_Boundary_'+res+'.txt'),'w')
        if self.Allelic == False:
            for chro, B_sub in self.boundary_filtered.items():
                for b in B_sub:
                    filtered_B_out.writelines(chro+'\t'+str(b)+'\n')
        else:
            for chro, B_sub in self.boundary_filtered.items():
                for b in B_sub:
                    filtered_B_out.writelines(chro[1:]+'\t'+str(b)+'\n')
        
        filtered_B_out.close()
        
        #Output Domains
        Domains_out = open(os.path.join(OutPath, prefix+'_Domain_'+res+'.txt'),'w')
        if self.Allelic == False:
            for chro in self.Domain_dict.keys():
                for i in range(len(self.Domain_dict[chro]['start'])):
                    Domains_out.writelines(chro+'\t'+str(self.Domain_dict[chro]['start'][i])+
                    '\t'+str(self.Domain_dict[chro]['end'][i])+'\n')
        else:
            for chro in self.Domain_dict.keys():
                for i in range(len(self.Domain_dict[chro]['start'])):
                    Domains_out.writelines(chro[1:]+'\t'+str(self.Domain_dict[chro]['start'][i])+
                    '\t'+str(self.Domain_dict[chro]['end'][i])+'\n')
        
        Domains_out.close()
        
        #Plotting
        if PLOT:
            print " Plotting ..."
        
            Pdf_out = os.path.join(OutPath, prefix+'_TADs_Plot_'+res+'.pdf')
            self.Plot_TAD(Pdf_out)
        
        print " Done!"
        
#==============================Loop Analysis module=================================


 
    def Peaks_Parameter(self):
        """
        Get Loop Calling Parameters
        
        pw : int
            Width of the interaction region surrounding the peak.According to experience,
            We suggest :
                1 for 20kb resolution.
                2 for 10kb resolution.
                4 for 5kb  resolution.
    
        ww : int
            The size of the donut sampled.
            We suggest :
                3 for 20kb resolution.
                5 for 10kb resolution.
                7 for 5kb  resolution.
        maxww : int
            Maximum donut size. 
    
        sig : float
            Significant Level.
    
        maxapart : int
            Maximum genomic distance between two loci.
            
        
        """
        if self.Res >= 20000:
            self.pw = 1
            self.ww = 3

        elif self.Res >= 10000:
            self.pw = 2
            self.ww = 5

        else:
            self.pw = 4
            self.ww = 7
        
        self.maxww = 20
        self.maxapart = 2000000
        self.sig = 0.05
    
    def lambdachunk(self, E):
        
        numbin = np.int(np.ceil(np.log(E.max()) / np.log(2) * 3 + 1))
        Pool = []
        for i in xrange(1, numbin + 1):
            if i == 1:
                lv = 0; rv = 1
            else:
                lv = np.power(2, ((i - 2)/3.))
                rv = np.power(2, ((i - 1)/3.))
            idx = np.where((E > lv) & (E < rv))[0]
            Pool.append((lv, rv, idx))
    
        return Pool
        
    def pcaller(self, M, cM, biases, IR, chromLen, Diags, cDiags, num, Allelic = False, Gap = None):
        """
            HICCUPS algorithm
        """
        pw = self.pw
        ww = self.ww
        sig = self.sig
        maxww = self.maxww
        maxapart = self.maxapart
        res = self.Res
        
        extDiags = {}
        for w in range(ww, maxww + 1):
            temp = []
            for i in xrange(num):
                OneDArray = Diags[i]
                extODA = np.zeros(chromLen - i + w*2)
                extODA[w:-w] = OneDArray
                temp.append(extODA)
            extDiags[w] = temp
    
        x = np.arange(ww, num)
        predictE = IR.predict(x)
        predictE[predictE < 0] = 0
        EDiags = []
        for i in xrange(x.size):
            OneDArray = np.ones(chromLen - x[i]) * predictE[i]
            EDiags.append(OneDArray)
    
        EM = sparse.diags(EDiags, x, format = 'csr')
    
        extCDiags = {}
        extEDiags = {}
        for w in range(ww, maxww + 1):
            tempC = []
            tempE = []
            for i in xrange(x.size):
                extODA_E = np.zeros(chromLen - x[i] + w*2)
                extODA_E[w:-w] = EDiags[i]
                tempE.append(extODA_E)
                extODA_C = np.zeros(chromLen - x[i] + w*2)
                extODA_C[w:-w] = cDiags[i]
                tempC.append(extODA_C)
            extCDiags[w] = tempC
            extEDiags[w] = tempE
    
        ps = 2 * pw + 1 # Peak Size
    
        Pool_Diags = {}
        Pool_EDiags = {}
        Pool_cDiags = {}
        Offsets_Diags = {}
        Offsets_EDiags = {}
    
        for w in range(ww, maxww + 1):
            ws = 2 * w + 1 # Window size
            ss = range(ws)
            Pool_Diags[w] = {}
            Pool_EDiags[w] = {}
            Pool_cDiags[w] = {}
            Offsets_Diags[w] = {}
            Offsets_EDiags[w] = {}
            for i in ss:
                for j in ss:
                    Pool_Diags[w][(i,j)] = []
                    Pool_EDiags[w][(i,j)] = []
                    Pool_cDiags[w][(i,j)] = []
                    Offsets_Diags[w][(i,j)] = np.arange(num) + (i - j)
                    Offsets_EDiags[w][(i,j)] = x + (i - j)
                    for oi in np.arange(num):
                        if Offsets_Diags[w][(i,j)][oi] >= 0:
                            starti = i
                            endi = i + chromLen - Offsets_Diags[w][(i,j)][oi]
                        else:
                            starti = i - Offsets_Diags[w][(i,j)][oi]
                            endi = starti + chromLen + Offsets_Diags[w][(i,j)][oi]
                        Pool_Diags[w][(i,j)].append(extDiags[w][oi][starti:endi])
                    for oi in xrange(x.size):
                        if Offsets_EDiags[w][(i,j)][oi] >= 0:
                            starti = i
                            endi = i + chromLen - Offsets_EDiags[w][(i,j)][oi]
                        else:
                            starti = i - Offsets_EDiags[w][(i,j)][oi]
                            endi = starti + chromLen + Offsets_EDiags[w][(i,j)][oi]
                        Pool_EDiags[w][(i,j)].append(extEDiags[w][oi][starti:endi])
                        Pool_cDiags[w][(i,j)].append(extCDiags[w][oi][starti:endi])
                
        ## Peak Calling ...    
        xi, yi = M.nonzero()
        Mask = ((yi - xi) >= ww) & ((yi - xi) <= (maxapart // res))
        xi = xi[Mask]
        yi = yi[Mask]
        if Allelic != False:
            H = M.toarray()
            Non_Gap = np.ones(xi.shape, dtype = bool)
            N = len(xi)
            print "Allelic HiC Loops Calling ..."
            print " Remove Gap and Blanking area ..."
            for i in range(N):
                #Gap region                
                if xi[i] in Gap and yi[i] in Gap:
                    Non_Gap[i]  = False
                
                #Blanking region
                try:
                    left = H[xi[i] - 1][yi[i]]
                except:
                    left = 1
                try:
                    right = H[xi[i] - 1][yi[i]]
                except:
                    right = 1
                try:
                    top = H[xi[i]][yi[i] + 1]
                except:
                    top = 1
                try:
                    bottom = H[xi[i]][yi[i] - 1]
                except:
                    bottom = 1
                if left * right * top * bottom == 0:
                    Non_Gap[i] = False
            xi = xi[Non_Gap]
            yi = yi[Non_Gap]
        else:
            print "Non-Allelic HiC Loops Calling ..."
                
        flocals = ['K', 'Y']
        bSV = {}; bEV = {}
        for fl in flocals:
            bSV[fl] = np.zeros(xi.size)
            bEV[fl] = np.zeros(xi.size)
    
        #logger.info('Observed Contact Number: %d', xi.size)
        print 'Observed Contact Number: %d' % xi.size
    
        RefIdx = np.arange(xi.size)
        RefMask = np.ones_like(xi, dtype = bool)
    
        iniNum = xi.size
    
        #logger.info('Two local neighborhoods, two expected matrices ...')
        print 'Two local neighborhoods, two expected matrices...'
        for w in range(ww, maxww + 1):
            ws = 2 * w + 1
            bS = {}; bE = {}
            for fl in flocals:
                bS[fl] = sparse.csr_matrix((chromLen, chromLen))
                bE[fl] = sparse.csr_matrix((chromLen, chromLen))
            Reads = sparse.csr_matrix((chromLen, chromLen))
            #logger.info('    Current window width: %s' % w)
            print '    Current window width: %s' % w
            P1 = set([(i,j) for i in range(w-pw, ps+w-pw) for j in range(w-pw, ps+w-pw)]) # Center Peak Region
            P_1 = set([(i,j) for i in range(w+1, ws) for j in range(w)])
            P_2 = set([(i,j) for i in range(w+1, ps+w-pw) for j in range(w-pw, w)])
            P2 = P_1 - P_2 # Lower-left Region
        
            for key in Pool_Diags[w]:
                if (key[0] != w) and (key[1] != w) and (key not in P1) and (key not in P2):
                    bS['K'] = bS['K'] + sparse.diags(Pool_cDiags[w][key], Offsets_EDiags[w][key], format = 'csr')
                    bE['K'] = bE['K'] + sparse.diags(Pool_EDiags[w][key], Offsets_EDiags[w][key], format = 'csr')
                if key in P2:
                    bS['K'] = bS['K'] + sparse.diags(Pool_cDiags[w][key], Offsets_EDiags[w][key], format = 'csr')
                    bE['K'] = bE['K'] + sparse.diags(Pool_EDiags[w][key], Offsets_EDiags[w][key], format = 'csr')
                    bS['Y'] = bS['Y'] + sparse.diags(Pool_cDiags[w][key], Offsets_EDiags[w][key], format = 'csr')
                    bE['Y'] = bE['Y'] + sparse.diags(Pool_EDiags[w][key], Offsets_EDiags[w][key], format = 'csr')
                    Reads = Reads + sparse.diags(Pool_Diags[w][key], Offsets_Diags[w][key], format = 'csr')
                
            Txi = xi[RefIdx]
            Tyi = yi[RefIdx]
            RNums = np.array(Reads[Txi, Tyi]).ravel()
            EIdx = RefIdx[RNums >= 16]
            #logger.info('    Valid Contact Number: %d', EIdx.size)
            print '    Valid Contact Number: %d' % EIdx.size
            Valid_Ratio = EIdx.size/float(iniNum)
            #logger.info('    Valid Contact Ratio: %.3f', Valid_Ratio)
            print '    Valid Contact Ratio: %.3f' % Valid_Ratio
            Exi = xi[EIdx]
            Eyi = yi[EIdx]
            for fl in flocals:
                bSV[fl][EIdx] = np.array(bS[fl][Exi, Eyi]).ravel()
                bEV[fl][EIdx] = np.array(bE[fl][Exi, Eyi]).ravel()
                
            RefIdx = RefIdx[RNums < 16]
            
            iniNum = RefIdx.size
        
            if Valid_Ratio < 0.1:
                #logger.info('    Ratio of valid contact is too small, break the loop ...')
                print '    Ratio of valid contact is too small, break the loop ...'
                break
        
            #logger.info('    Continue ...')
            print '    Continue ...'
            #logger.info('    %d Contacts will get into next loop ...', RefIdx.size)
            print '    %d Contacts will get into next loop ...' % RefIdx.size
        RefMask[RefIdx] = False
    
        Mask = (bEV['K'] != 0) & (bEV['Y'] != 0) & RefMask
        xi = xi[Mask]
        yi = yi[Mask]
        bRV = {}
        for fl in flocals:
            bRV[fl] = bSV[fl][Mask] / bEV[fl][Mask]
    
        bR = {}
        for fl in flocals:
            bR[fl] = sparse.csr_matrix((chromLen, chromLen))
            bR[fl][xi, yi] = bRV[fl]
    
        ## Corrected Expected Matrix
        cEM = {}
        for fl in flocals:
            cEM[fl] = EM.multiply(bR[fl])
    
        #logger.info('Poisson Models and Benjamini-Hochberg Correcting for lambda chunks ...')
        print 'Poisson Models and Benjamini-Hochberg Correcting for lambda chunks ...'
        Description = {'K': 'Donut backgrounds', 'Y': 'Lower-left backgrounds'}
        xpos = {}; ypos = {}; Ovalues = {}; Evalues = {}
        Fold = {}; pvalues = {}; qvalues = {}
        gaps = set(np.where(np.array(M.sum(axis=1)).ravel() == 0)[0])
        for fl in flocals:
            #logger.info('    %s ...', Description[fl])
            print '    %s ...' % Description[fl]
            xi, yi = cEM[fl].nonzero()
            Evalues[fl] = np.array(cEM[fl][xi, yi]).ravel() * biases[xi] * biases[yi]
            Mask = (Evalues[fl] > 0)
            Evalues[fl] = Evalues[fl][Mask]
            xi = xi[Mask]
            yi = yi[Mask]
            Ovalues[fl] = np.array(M[xi, yi]).ravel()
            Fold[fl] =  Ovalues[fl] / Evalues[fl]
            #logger.info('    Valid contact number: %d', xi.size)
            print '    Valid contact number: %d', xi.size
        
            pvalue = np.ones(xi.size)
            qvalue = np.ones(xi.size)
        
            #logger.info('    Lambda chunking ...')
            print '    Lambda chunking ...'
            chunks = self.lambdachunk(Evalues[fl])
            #logger.info('    Number of chunks: %d', len(chunks))
            print '    Number of chunks: %d', len(chunks)
            for chunk in chunks:
                #logger.debug('        lv: %.4g, rv: %.4g, Num: %d', chunk[0], chunk[1], chunk[2].size)
                print '        lv: %.4g, rv: %.4g, Num: %d' % (chunk[0], chunk[1], chunk[2].size)
                if chunk[2].size > 0:
                    Poiss = poisson(chunk[1])
                    #logger.debug('        Assign P values ...')
                    print '        Assign P values ...'
                    chunkP = 1 - Poiss.cdf(Ovalues[fl][chunk[2]])
                    pvalue[chunk[2]] = chunkP
                    #logger.debug('        Multiple testing ...')
                    print '        Multiple testing ...'
                    cResults = multipletests(chunkP, alpha = sig, method = 'fdr_bh')
                    cP = cResults[1] # Corrected Pvalue
                    qvalue[chunk[2]] = cP
                else:
                    #logger.debug('        Skipping ...')
                    print '        Skipping ...'
        
            reject = qvalue <= sig
            qvalue = qvalue[reject]
            pvalue = pvalue[reject]
            Ovalues[fl] = Ovalues[fl][reject]
            Evalues[fl] = Evalues[fl][reject]
            Fold[fl] = Fold[fl][reject]
            xi = xi[reject]
            yi = yi[reject]
        
            #logger.info('    Remove Gap Effects ...')
            print '    Remove Gap Effects ...'
        
            if len(gaps) > 0:
                fIdx = []
                for i in xrange(xi.size):
                    lower = (xi[i] - 5) if (xi[i] > 5) else 0
                    upper = (xi[i] + 5) if ((xi[i] + 5) < chromLen) else (chromLen - 1)
                    cregion_1 = range(lower, upper)
                    lower = (yi[i] - 5) if (yi[i] > 5) else 0
                    upper = (yi[i] + 5) if ((yi[i] + 5) < chromLen) else (chromLen - 1)
                    cregion_2 = range(lower, upper)
                    cregion = set(cregion_1) | set(cregion_2)
                    intersect = cregion & gaps
                    if len(intersect) == 0:
                        fIdx.append(i)
        
                xi = xi[fIdx]
                yi = yi[fIdx]
                Ovalues[fl] = Ovalues[fl][fIdx]
                pvalue = pvalue[fIdx]
                qvalue = qvalue[fIdx]
                Fold[fl] = Fold[fl][fIdx]
                Evalues[fl] = Evalues[fl][fIdx]
        
            xpos[fl] = xi
            ypos[fl] = yi
            pvalues[fl] = pvalue
            qvalues[fl] = qvalue
    
        #logger.info('Combine two local filters ...')
        print 'Combine two local filters ...'
    
        preDonuts = dict(zip(zip(xpos['K']*res, ypos['K']*res), zip(Ovalues['K'], Fold['K'], pvalues['K'], qvalues['K'])))
        preLL = dict(zip(zip(xpos['Y']*res, ypos['Y']*res), zip(Ovalues['Y'], Fold['Y'], pvalues['Y'], qvalues['Y'])))
    
        commonPos = set(preDonuts.keys()) & set(preLL.keys())
        Donuts = {}; LL = {}
        for pos in commonPos:
            Donuts[pos] = preDonuts[pos]
            LL[pos] = preLL[pos]
    
        return Donuts, LL        
    
    def bias_handle(self, bias):
        
        bias = bias.reshape(bias.shape[0],)
        bias[np.isnan(bias)] = 1.0
        return bias
        
    def CallPeaks(self, outfil, Allelic = False):
        """
            Identify chromotin loops from Hi-C data by HICCUPS algorithm.
        
        Parameters
        ----------
        outfil : str
            Output file
        
        Allelic : 
                    False for Traditional Loop calling.
                    Maternal for Haplotype-resolved Loop calling (in Maternal)
                    Paternal for Haplotype-resolved Loop calling (in Paternal)
        
        """
        print "===============Calling Loops==============="
        self.cooler = cooler.Cooler(self.cooler_fil)
        self.Matrix = {}
        #Chromosome List
        if Allelic == False:
            chroms = self.cooler.chromnames
        elif Allelic == 'Maternal':
            chroms = [i for i in self.cooler.chromnames if i.startswith('M')]
        elif Allelic == 'Paternal':
            chroms = [i for i in self.cooler.chromnames if i.startswith('P')]
        else:
            raise Exception ('Unkonwn key word %s, Only Maternal, Paternal, False allowed' % Allelic)
        self.chroms = chroms
        #Gap List
        if Allelic == False:
            pass
        else:
            if self.Gap_file == None:
                raise Exception ('Gap file needed for haplotype-resolved loop calling ...')
            Gap = np.load(self.Gap_file, allow_pickle = True)
            Gap = Gap[str(self.Res)][()]
            gap = {}
            for chro in chroms:
                gap[chro] = Gap[chro]

            
        OF = open(os.path.join(outfil),'w')
        head = '\t'.join(['chromLabel', 'loc_1', 'loc_2', 'IF', 'D-Enrichment', 'D-pvalue', 'D-qvalue',
                          'LL-Enrichment', 'LL-pvalue', 'LL-qvalue']) + '\n'
    
        OF.write(head)        
        #---------------Initialize Parameters----------------
        self.Peaks_Parameter()
        
        for chro in self.chroms:
            print "Chromsome %s ..." % chro
            if Allelic == False:
                H = self.cooler.matrix(balance = False).fetch(chro)
                self.Matrix[chro] = H
                cH = self.cooler.matrix(balance = True).fetch(chro)
                cH = np.nan_to_num(cH)
                tmp = self.cooler.bins().fetch(chro)['weight'].values
                mask = np.logical_not((tmp == 0)) | np.isnan(tmp)
                biases = np.zeros_like(tmp)
                biases[mask] = 1 / tmp[mask]
            else:
                H = self.cooler.matrix(balance = False).fetch(chro)
                self.Matrix[chro] = H
                cH = copy.deepcopy(H)
                gap_sub = gap[chro]
                biases = np.ones(H.shape[0],)
            
            H = H - np.diag(H.diagonal())
            print "Customize sparse Matrix"
            chromLen = H.shape[0]
            num = self.maxapart // self.Res + self.maxww + 1
            Diags = [np.diagonal(H,i) for i in np.arange(num)]
            M = sparse.diags(Diags, np.arange(num), format = 'csr')
            x = np.arange(self.ww, num)
            y = []
            cDiags = []
            for i in x:
                diag = np.diagonal(cH,i)
                y.append(diag.mean())
                cDiags.append(diag)
            cM = sparse.diags(cDiags, x, format = 'csr')
            IR = isotonic.IsotonicRegression(increasing = 'auto')
            IR.fit(x, y)
            
            del H, cH
            if Allelic == False:
                Donuts, LL = self.pcaller(M = M, cM = cM, biases = biases, IR = IR, 
                                          chromLen = chromLen, Diags = Diags,
                                          cDiags = cDiags, num = num, Allelic = False)
            else:
                Donuts, LL = self.pcaller(M = M, cM = cM, biases = biases, IR = IR, 
                                          chromLen = chromLen, Diags = Diags, Gap = gap_sub,
                                          cDiags = cDiags, num = num, Allelic = True)
            
            for i in Donuts:
                lineFormat = '%s\t%d\t%d\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\n'
                if self.Allelic == False:
                    contents = (chro,) + i + Donuts[i] + LL[i][1:]
                else:
                    contents = (chro[1:],)+ i + Donuts[i] + LL[i][1:]
                line = lineFormat % contents
                OF.write(line)
                    
        OF.flush()
        OF.close()
        
        print "HICCUPS Done!"
        
    
    def Loop_Selecting(self, input_fil, output_fil):
        """
        """
        print "=================================="
        print "Loop Selecting ..."
        
        f = open(input_fil,'r')
        head = ['chromLabel', 'loc_1', 'loc_2', 'IF', 'D-Enrichment', 
                'D-pvalue', 'LL-Enrichment', 'LL-pvalue', 'LL-qvalue']
        o = open(output_fil, 'w')
        o.writelines('\t'.join(head)+'\n')
        
        for line in islice(f, 1, None):
            l = line.strip().split()
            chro = l[0]
            bin1 = int(l[1]) // 40000
            bin2 = int(l[2]) // 40000
            M = self.Matrix[chro]
            IF = M[bin1][bin2]
            
            distance = M.diagonal(bin1 - bin2).copy()
            distance.sort()
            
            index = bisect.bisect_left(distance, IF)
            ratio = index / len(distance)
            
            if ratio < self.ratio or IF < self.LoopStrength:
                continue
            else:
                o.writelines(line)
        f.close()
        o.close()
                
    
    def center_sites(self, lists):
        sum_x = 0
        sum_y = 0
        for x in lists:
            sum_x += x[1]
            sum_y += x[2]
        n = len(lists)
        return [float(sum_x)/n, float(sum_y) / n]
    
    def distance(self, sites_1,sites_2):
        return math.sqrt((sites_1[0]-sites_2[1])**2 + (sites_1[1]-sites_2[2])**2)
    
    def peakcluster(self, loop_sites, dis, chrom):
        classes = []
        for i in chrom:
            c_loops = sorted(loop_sites[loop_sites['chr'] == i], key = lambda x:x[1])
            while True:
                c_class = []
                c_class.append(c_loops[0])
                c_loops.remove(c_loops[0])
                center = self.center_sites(c_class)
                for loop in c_loops:
                    if self.distance(center, loop) <= dis:
                        c_class.append(loop)
                        center = self.center_sites(c_class)
                        c_loops.remove(loop)
                classes.append(c_class)
                if len(c_loops) == 0:
                    break
        
        return classes
    
    def filter_initial(self,cls):
        loop = []
        for outs in cls:
            lens = len(outs)
            if lens >= 2:
                outs = sorted(outs, key = lambda x:x[3])
                loop.append((outs[0][0], outs[0][1], outs[0][2], outs[0][3], lens))
            else:
                loop.append((outs[0][0], outs[0][1], outs[0][2], outs[0][3], lens))
        return np.array(loop, dtype = [ ('chr', '<S4'), ('S1', '<i4'), ('E1', '<i4'), ('Q', '<f8'), ('sums', '<f8')])
        
    def filter_next(self,cls):
        loop = []
        for outs in cls:
            lens = len(outs)
            sums = 0
            for x in outs:
                sums += x[4]
            if  lens >= 2 :
                outs = sorted(outs, key=lambda x:x[3])
                loop.append((outs[0][0], outs[0][1], outs[0][2], outs[0][3], sums))
            else :
                loop.append((outs[0][0], outs[0][1], outs[0][2], outs[0][3], sums))
        return np.array(loop, dtype = [ ('chr', '<S4'), ('S1', '<i4'), ('E1', '<i4'), ('Q', '<f8'), ('sums', '<f8')])
        
    def LoopCluster(self, rawfil, weight_q_value = 0.0001):
        """
            Use the Cluster algorithm for final loop select.
            
        Parameters
        ----------
        rawfil : fil
            HICCUPS results file.
        
        weight_q_value : float
            threshold for cluster.
            
        """
        print "=================================="
        print "Loop Clustering ..."
        loopType = np.dtype({'names':['chr','S1','E1','Q-value'],
                             'formats':['S4', np.int, np.int,np.float]})
        
        loops = np.loadtxt(rawfil, dtype = loopType, usecols = [0,1,2,9], skiprows = 1)
        
        initial_distance = self.Res * math.sqrt(2) + 1000
        chroms = list(set(loops['chr']))
        
        cluster_initial = self.peakcluster(loops, initial_distance, chroms)
        loop_inital = self.filter_initial(cluster_initial)
        
        while True:
            cluster_loop = self.peakcluster(loop_inital,initial_distance * 2, chroms)
            loop_final = self.filter_next(cluster_loop)
            
            if len(loop_inital) == len(loop_final):
                break
            loop_inital = loop_final
        
        path = os.path.split(rawfil)[0]
        if path == '':
            path = './'
        fil = os.path.split(rawfil)[1]
        
        cluster_fil = os.path.join(path,'Cluster_'+fil)
        out_fil = open(cluster_fil,'w')
        out_fil.writelines('chr\tstart\tend\tIF\tweight_Q-value\taggregateNum\n')
        
        if self.Allelic == False:          
            for outs in loop_final:
                q_value = outs[3] / (10 ** outs[4])
                x, y = outs[1] // self.Res, outs[2] // self.Res
                if q_value < weight_q_value:
                    lst = [outs[0],outs[1],outs[2],self.Matrix[outs[0]][x][y],q_value,outs[4]]
                    strs = '\t'.join(map(str,lst)) + '\n'
                    out_fil.write(strs)
        else:
            weighted_Loops = []                       
            for outs in loop_final:
                if self.Allelic == 'Maternal':
                    M = self.Matrix['M'+outs[0]]
                else:
                    M = self.Matrix['P'+outs[0]]
                x, y = outs[1] // self.Res, outs[2] // self.Res
                q_value = outs[3] / (10 ** outs[4])
                if q_value< weight_q_value:
                    lst = (outs[0],outs[1],outs[2],M[x][y],q_value,outs[4])
                    weighted_Loops.append(lst)
            
            weighted_type = np.dtype({'names':['chr','start','end','IF','w_q','level'],
                                      'formats':['S4',np.int,np.int,np.float,np.float64,np.int]})
            weighted_Loops = np.array(weighted_Loops,dtype = weighted_type)
            
            Loop_threshould = {}
            chros = set(weighted_Loops['chr'])
            
            #Transform zero ->  10 ** -20 for compute
            weighted_Loops['w_q'][weighted_Loops['w_q'] == 0] = 10 ** (-20)
            
            for chro in chros:
                tmp_loops = weighted_Loops[weighted_Loops['chr'] == chro]
                tmp_threshould = np.percentile(tmp_loops['IF'] * -np.log10(tmp_loops['w_q']),15)
                Loop_threshould[chro] = tmp_threshould
                
            for outs in weighted_Loops:
                if (outs['IF'] * -np.log10(outs['w_q'])) < Loop_threshould[outs['chr']]:
                    continue
                else:
                    lst = tuple(outs)
                    strs = '\t'.join(map(str,lst)) + '\n'
                    out_fil.write(strs)
                
        out_fil.close()
        print "Done !"
        self.out_fil = cluster_fil
    
    
    
    def Loading_Loops(self):
        """
        """
        dtype = np.dtype({'names':['chr', 'start', 'end'],
                          'formats':['S4', np.int, np.int]})
        
        self.loops = np.loadtxt(self.out_fil, skiprows = 1, usecols = (0,1,2), dtype = dtype)
        
    
    
    
    
    def Plot_Loops(self, outpdf, length = 4000000):
        """
            Plotting Loop with HeatMap
        """
        pp = PdfPages(outpdf)
        cmap = self.Cmap_Setting()
        self.fig_settings()
        self.Loading_Loops()
        if self.Allelic == False:
            self.Matrix = {}
            for chro in self.chroms:
                M = self.cooler.matrix(balance = True).fetch(chro)
                M = np.nan_to_num(M)
                self.Matrix[chro] = M
        else:
            pass
        
        for chro in self.chroms:
            Matrix = self.Matrix[chro]
            N = Matrix.shape[0]
            
            interval = length // self.Res
            
            startHiC = 0
            Len = N // interval
            if self.Allelic == False:
                loops = self.loops[self.loops['chr'] == chro]
            else:
                loops = self.loops[self.loops['chr'] == chro[1:]]
            
            for idx in range(Len):
                endHiC = startHiC + interval
                
                fig = plt.figure(figsize = self.figure_Dict['size'])
                ax = fig.add_axes([self.figure_Dict['Left'], self.figure_Dict['HB'],
                                   self.figure_Dict['width'], self.figure_Dict['HH']])
                M = Matrix[startHiC:endHiC, startHiC:endHiC]
                
                mask = (loops['start'] >= startHiC * self.Res) & (loops['end'] <= endHiC * self.Res)
                
                tmp_loops = loops[mask]
                
                nonzero = M[np.nonzero(M)]
                if nonzero.size <= 100 or tmp_loops.size == 0:
                    plt.close(fig)
                    startHiC = endHiC
                    continue
                
                vmax = np.percentile(nonzero, 95)
                sc = ax.imshow(M, cmap = cmap, aspect = 'auto', interpolation = 'none',
                               extent = (0, M.shape[0], 0, M.shape[0]), vmax = vmax, origin = 'lower')
                for lp in tmp_loops:
                    s_t = lp['start'] // self.Res - startHiC
                    s_e = lp['end'] // self.Res - startHiC
                    ax.scatter(s_t + 0.5, s_e + 0.5, color = '', edgecolors = 'b', s = 10)
                
                
                ticks = list(np.linspace(0,interval,5).astype(int))
                pos = [((startHiC + t) * self.Res) for t in ticks]
                labels = [self.properU(p) for p in pos]
                ax.set_xticks(ticks)
                ax.set_xticklabels(labels)
                ax.set_yticks(ticks)
                ax.set_yticklabels(labels)
                ax.set_xlim(0,(endHiC - startHiC))
                ax.set_ylim(0,(endHiC - startHiC))
                if self.Allelic == False:
                    ax.set_xlabel('Chr'+chro, size = 14)
                else:
                    ax.set_xlabel('Chr'+chro[1:], size = 14)
                
                ax = fig.add_axes([self.figure_Dict['Left'] + self.figure_Dict['width'] + 0.02, 
                                   self.figure_Dict['HB'], 0.01, self.figure_Dict['HH']])
                fig.colorbar(sc,cax = ax)
                
                pp.savefig(fig)
                plt.close(fig)
                startHiC = endHiC
        pp.close()
                
    
    def run_Loops(self, OutPath, plot = False):
        """
            run Loop calling.
        
        Parameters
        ----------
        OutPath : path
            Output Path
            
        """
        if os.path.exists(OutPath):
            pass
        else:
            os.mkdir(OutPath)
        
        res = self.properU(self.Res)
        prefix = os.path.split(OutPath.rstrip('/'))[-1]
        outfil = os.path.join(OutPath,prefix+'_Loops_'+res+'.txt')
        if self.Allelic == False:
            select_fil = os.path.join(OutPath, 'Selected_'+prefix+'_Loops_'+res+'.txt')
            self.CallPeaks(outfil,Allelic=self.Allelic)
            self.Loop_Selecting(input_fil = outfil, output_fil = select_fil)
            self.LoopCluster(rawfil = select_fil)
        else:
            self.CallPeaks(outfil,Allelic=self.Allelic)
            self.LoopCluster(rawfil = outfil)
        
        if plot:
            print "Plotting..."
            
            Pdf_out = os.path.join(OutPath, prefix+'_Loops_Plot_'+res+'.pdf')
            self.Plot_Loops(Pdf_out)
        
        print 'Done !'
        
