# -*- coding: utf-8 -*-
"""
Created on Mon May 13 21:24:03 2019

@author: han-luo
"""

from __future__ import division
from scipy import sparse
import cooler
from cooler.util import binnify
from cooler.reduce import CoolerMerger
from cooler.api import Cooler
from cooler import create_cooler
import os, tempfile, logging, math, copy
import pandas as pd
import numpy as np
import subprocess, gc, time

log = logging.getLogger(__name__)

def create_from_unordered(cool_uri, bins, chunks, columns=None, dtypes=None, mergebuf=int(20e6),
                         delete_temp=True, temp_dir=None, **kwargs):
    """
    Create a Cooler in two passes via an external sort mechanism. In the first 
    pass, a sequence of data chunks are processed and sorted in memory and saved
    to temporary Coolers. In the second pass, the temporary Coolers are merged 
    into the output. This way the individual chunks do not need to be provided
    in any particular order.
    
    Parameters
    ----------
    cool_uri : str
        Path to Cooler file or URI to Cooler group. If the file does not exist,
        it will be created.
    bins : DataFrame
        Segmentation of the chromosomes into genomic bins. May contain 
        additional columns.
    chunks : iterable of DataFrames
        Sequence of chunks that get processed and written to separate Coolers 
        and then subsequently merged.
    columns : sequence of str, optional
        Specify here the names of any additional value columns from the input 
        besides 'count' to store in the Cooler. The standard columns ['bin1_id', 
        'bin2_id', 'count'] can be provided, but are already assumed and don't 
        need to be given explicitly. Additional value columns provided here will 
        be stored as np.float64 unless otherwised specified using `dtype`.
    dtypes : dict, optional
        Dictionary mapping column names to dtypes. Can be used to override the
        default dtypes of ``bin1_id``, ``bin2_id`` or ``count`` or assign
        dtypes to custom value columns. Non-standard value columns given in
        ``dtypes`` must also be provided in the ``columns`` argument or they
        will be ignored.
    assembly : str, optional
        Name of genome assembly.
    mode : {'w' , 'a'}, optional [default: 'w']
        Write mode for the output file. 'a': if the output file exists, append
        the new cooler to it. 'w': if the output file exists, it will be
        truncated. Default is 'w'.
    metadata : dict, optional
        Experiment metadata to store in the file. Must be JSON compatible.
    mergebuf : int, optional
        Maximum number of records to buffer in memory at any give time during 
        the merge step.
    delete_temp : bool, optional
        Whether to delete temporary files when finished. 
        Useful for debugging. Default is False.
    temp_dir : str, optional
        Create temporary files in this directory.
    See also
    --------
    sanitize_records
    sanitize_pixels
    """
    bins = bins.copy()
    bins['chrom'] = bins['chrom'].astype(object)

    tf = tempfile.NamedTemporaryFile(
                suffix='.multi.cool', 
                delete=delete_temp,
                dir=temp_dir)
        
    uris = []
    for i, chunk in enumerate(chunks):
        uri = tf.name + '::' + str(i)
        uris.append(uri)
        #log.info('Writing chunk {}: {}'.format(i, uri))
        create_cooler(uri, bins, chunk, columns=columns, mode='a', boundscheck=False,
                      triucheck=False, dupcheck=False, ensure_sorted=False, ordered=True,
                      dtypes=dtypes)
        
    chunks = CoolerMerger([Cooler(uri) for uri in uris], mergebuf)

    #log.info('Merging into {}'.format(cool_uri))
    create_cooler(cool_uri, bins, chunks, columns=columns, dtypes=dtypes, ordered=True,
                  **kwargs)



class NPZ2Cooler(object):
    """
    Load bin-level Hi-C data of NPZ format, and save it into cooler files.
    
    This class was created with the help of xiaoTao Wang's original package which 
    support an effective  way for searching chromatin loops.
    More details : (https://github.com/XiaoTaoWang/HiCPeaks)
    
    Thanks a lot, my old brother :)
    
    Parameters
    ----------
    datasets : dict, {resolution(int):NPZ_file}
        *resolution* should be in base-pair unit. *data_path* indicates the
        absolute data path.
        
    outfil : str
        Path of the output Cooler file.
    
    assembly : str
        Genome assembly name. (Default: hg38)
        
    chromsizes_file : str
        Path to the file containing chromosome size information.
    
    chroms : list
        List of chromosome labels. Only Hi-C data within the specified chromosomes
        will be included. Specially, '#' stands for chromosomes with numerical
        labels. If an empty list is provided, all chromosome data will be loaded.
        (Default: ['#', 'X'])
        
    onlyIntra : bool
        If specified, only include intra-chromosomal data.
    
    dtype : {'int', float}
        The desired data type for your contact matrices.
    """    

    def __init__(self, datasets, outfil, genomeSizes_file, chroms = ['#', 'X'],
                 onlyIntra = True, dtype = 'int'):
        
        self.outfil = os.path.abspath(os.path.expanduser(outfil))
        self.chroms = set(chroms)
        self.onlyIntra = onlyIntra
        self.genomeSize_file = genomeSizes_file
        data = datasets
        
        chromsizes = self.readChromSize()
        chromlist = Sort_Chromosomes(chromsizes.keys())
        lengths = [chromsizes[i] for i in chromlist]
        self.chromsizes = pd.Series(data = lengths, index = chromlist)
        
        ## Store the data index into Map, for cooler creating...
        self.Map = {}
        for res in data:
            self.Map[res] = {}
            lib = data[res]
            for i in lib.keys():
                if (not '_' in i) and ((not self.chroms) or (i.isdigit() and '#' in self.chroms) or (i in self.chroms)):
                    c1 = c2 = i
                    self.Map[res][(c1, c2)] = lib
                else:
                    tmp = i.split('_')
                    if len(tmp) != 2:
                        continue
                    else:
                        c1, c2 = tmp
                        check1 = ((not self.chroms) or (c1.isdigit() and '#' in self.chroms) or (c1 in self.chroms))
                        check2 = ((not self.chroms) or (c2.isdigit() and '#' in self.chroms) or (c2 in self.chroms))
                        if check1 and check2:
                            self.Map[res][c1, c2] = lib
        
        print 'Extract and save data into cooler format for each resolution ...'
        log.log(21, '    Extract and save data into cooler format for each resolution ...')
        for res in self.Map:
            print '    Current resolution : {}bp'.format(res)
            log.log(21, '    Current resolution : {}bp'.format(res))
            byres = self.Map[res]
            subset = []
            for c1, c2 in byres:
                subset.extend([c1, c2])
            subset = set(subset)
            Bool = [(i in subset) for i in self.chromsizes.index]
            chromsizes = self.chromsizes[Bool]
            bin_cumnums = self.binCount(chromsizes, res)
            print '    Generate bin table ...'
            log.log(21, '    Generate bin table ...')
            bintable = binnify(chromsizes, res)
            pixels = self._generator(byres, bin_cumnums)
            
            if os.path.exists(self.outfil):
                mode = 'a'
            else:
                mode = 'w'
            
            if dtype == 'int':
                dtypes = {'count': np.int32}
            else:
                dtypes = {'count': np.float64}
            
            cooler_uri = '{}::{}'.format(self.outfil, res)
            
            print "    save to cooler ..."
            log.log(21, '    save to cooler ...')
            if self.onlyIntra:
                create_cooler(cooler_uri, bintable, pixels, mode=mode,boundscheck=False,
                              triucheck=False, dupcheck=False, ensure_sorted=False,
                              ordered=True,metadata={'onlyIntra':str(self.onlyIntra)}, dtypes=dtypes)
            else:
                create_from_unordered(cooler_uri, bintable, pixels, mode=mode, delete_temp=True,
                                      metadata={'onlyIntra':str(self.onlyIntra)},boundscheck=False,
                                      triucheck=False, dupcheck=False,ensure_sorted=False, dtypes=dtypes)
            
            print '    Done!'
            log.log(21, '    Done!')

    def readChromSize(self):
        
        chromsize = {}
        with open(self.genomeSize_file,'r') as source:
            for line in source:
                parse = line.strip().split()
                c, s = parse[0].lstrip('chr'), parse[1]
                check = ((not self.chroms) or (c.isdigit() and ('#' in self.chroms)) or (c in self.chroms))
                if check:
                    chromsize[c] = int(s)
        
        return chromsize
    

    def binCount(self, chromsizes, res):
        
        def _each(chrom):
            clen = chromsizes[chrom]
            n_bins = int(np.ceil(clen / res))
            return n_bins
        
        data = [_each(c) for c in chromsizes.index]
        n_bins = pd.Series(data, index=chromsizes.index)
        cum_n = n_bins.cumsum()
        
        return cum_n

    def _generator(self, byres, bin_cumnums):
        
        S_dtype = np.dtype({'names':['bin1', 'bin2', 'IF'],
                                    'formats':[np.int, np.int, np.float]})
                                    
        for i in range(self.chromsizes.size):
            for j in range(i, self.chromsizes.size):
                c1, c2 = self.chromsizes.index[i], self.chromsizes.index[j]
                if self.onlyIntra:
                    if c1!=c2:
                        continue
                if (c1,c2) in byres:
                    ci, cj = i, j
                else:
                    if (c2,c1) in byres:
                        c1, c2 = c2, c1
                        ci, cj = j, i
                    else:
                        continue
                
                if type(byres[(c1,c2)])==str:
                    data = np.loadtxt(byres[(c1,c2)], dtype=S_dtype)
                else:
                    # Make it compatible with TADLib and old version of runHiC
                    if c1!=c2:
                        data = byres[(c1,c2)][(c1+'_'+c2)]
                    else:
                        if c1 in byres[(c1,c2)].keys():
                            data = byres[(c1,c2)][c1]
                        else:
                            data = byres[(c1,c2)][(c1+'_'+c2)]

                x, y = data['bin1'], data['bin2']
                # Fast guarantee triu matrix
                if ci > cj:
                    x, y = y, x
                    ci, cj = cj, ci
                
                xLen = x.max() + 1
                yLen = y.max() + 1
                if ci != cj:
                    tmp = sparse.csr_matrix((data['IF'], (x,y)), shape=(xLen, yLen))
                else:
                    Len = max(xLen, yLen)
                    tmp = sparse.csr_matrix((data['IF'], (x,y)), shape=(Len, Len))
                    tmp = sparse.lil_matrix(tmp)
                    tmp[y,x] = tmp[x,y]
                    tmp = sparse.triu(tmp)
                
                x, y = tmp.nonzero()
                if ci > 0:
                    x = x + bin_cumnums[ci-1]
                if cj > 0:
                    y = y + bin_cumnums[cj-1]
                
                data = tmp.data

                current = pd.DataFrame({'bin1_id':x, 'bin2_id':y, 'count':data},
                                       columns=['bin1_id', 'bin2_id', 'count'])

                yield current



def Merge_beds(bed_lst):
    """
        Merging the  beds To Generate the Matrix.
    """
    Merge_cmd = ['cat'] + list(bed_lst)
    
    return Merge_cmd


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


def Load_Genome(genomeSize, chroms):
    """
        Loading the Genome Size
        Return the genomeSize dict.
    """
    
    genome = {}
    with open(genomeSize,'r') as f:
        for line in f:
            line = line.strip().split()
            c = line[0].lstrip('chr')
            check = ((not chroms) or (c.isdigit() and ('#' in chroms)) or (c in chroms))
            if check:
                genome[c] = int(line[1])
            else:
                continue
    
    return genome


def Load_HaplotypeGenome(genomeSize, chroms):
    """
        Loading the Genome Size
        Return the genomeSize dict
    """
    genome = {}
    with open(genomeSize,'r') as f:
        for line in f:
            line = line.strip().split()
            c = line[0].lstrip('chr')
            check = ((not chroms) or (c.isdigit() and ('#' in chroms)) or (c in chroms))
            if check:
                genome['M'+c] = int(line[1])
                genome['P'+c] = int(line[1])
            else:
                continue
            
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
    

def Get_Chro_Bins(genomeSize,Resolution, chroms):
    """
    """
    ChromBins = Load_Genome(genomeSize, chroms)
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


def Get_Chro_Bins_Haplotypes(genomeSize,Resolution, chroms):
    """
    """
    Haplotype_Bins = {}
    ChromBins = Load_Genome(genomeSize, chroms)
    Chro = Sort_Chromosomes(ChromBins)
    
    for c, l in ChromBins.items():
        ChromBins[c] = l // Resolution
    
    Haplotype_Chro = []
    for i in Chro:
        Haplotype_Chro.append('M'+i)
    for i in Chro:
        Haplotype_Chro.append('P'+i)
    
    for index, chro in enumerate(Haplotype_Chro):
        if index == 0:
            Haplotype_Bins[chro] = (0, ChromBins[chro[1:]])
        else:
            Haplotype_Bins[chro] = (Haplotype_Bins[Haplotype_Chro[index - 1]][1] +1, 
                                    ChromBins[chro[1:]] + Haplotype_Bins[Haplotype_Chro[index - 1]][1]+1)
    
    Haplotype_Bins_num = Haplotype_Bins[Haplotype_Chro[-1]][1]
    
    return Haplotype_Bins, Haplotype_Bins_num + 1
    

def WholeMatrixToSparseDict(Bins, Matrix):
    """
    """
    S_dtype = np.dtype({'names':['bin1', 'bin2', 'IF'],
                        'formats':[np.int, np.int, np.float]})
    chroms = Sort_Chromosomes(Bins.keys())
    
    W_Lib = {}
    for i in range(len(chroms)): 
        if i == len(chroms) - 1:
            chro = chroms[i]    
            start = Bins[chro][0]        
            end = Bins[chro][1] + 1                     
            W_Lib[chro] = Matrix[start:end, start:end]
        else:              
            for j in range(i, len(chroms)):
                chro1 = chroms[i]
                chro2 = chroms[j]           
                if chro1 == chro2:
                    start = Bins[chro1][0]              
                    end = Bins[chro1][1] + 1
                    W_Lib[chro1] = Matrix[start:end, start:end]
                else:         
                    key = chro1+'_'+chro2
                    start1 = Bins[chro1][0]
                    end1 = Bins[chro1][1] + 1
                    start2 = Bins[chro2][0]
                    end2 = Bins[chro2][1] + 1
                    W_Lib[key] = Matrix[start1:end1, start2:end2]
    
    for key, MM in W_Lib.items():
        if '_' not in key:
            Triu = np.triu(MM)
            x, y = np.nonzero(Triu)
            values = Triu[x, y]
            tmp = np.zeros(values.size, dtype = S_dtype)
            tmp['bin1'] = x
            tmp['bin2'] = y
            tmp['IF'] = values
            W_Lib[key] = tmp
        else:
            x, y = np.nonzero(MM)
            values = MM[x, y]
            tmp = np.zeros(values.size, dtype = S_dtype)
            tmp['bin1'] = x
            tmp['bin2'] = y
            tmp['IF'] = values
            W_Lib[key] = tmp
    
    return W_Lib

def IntraMatrixToSparseDict(Dict):
    """
    """
    S_dtype = np.dtype({'names':['bin1', 'bin2', 'IF'],
                        'formats':[np.int, np.int, np.float]})
    Sparse_Dict = {}
    for chro in Dict.keys():
        M = np.triu(Dict[chro])
        x, y = np.nonzero(M)
        values = M[x, y]
        tmp = np.zeros(values.size, dtype = S_dtype)
        tmp['bin1'] = x
        tmp['bin2'] = y
        tmp['IF'] = values
        Sparse_Dict[chro] = tmp
    
    return Sparse_Dict



def TraditionalMatrixBuilding(bed_IO, genomeSize, wholeRes, localRes, chroms):
    """
        Tranditional Matrix construction.
        
        bed_IO : input instream
        
        genomeSize : genomeSize file
        
        wholeRes : list or tuple
                    Resolution for whole-genome contact matrix construction
        
        LocalRes : list or tuple
                    Resolution for intra-chromosome contact matrix construction
        
        chroms : list or tuple
                    List of chromosome labels. Only Hi-C data within the specified chromosomes
                    will be included. Specially, '#' stands for chromosomes with numerical
                    labels. If an empty list is provided, all chromosome data will be loaded.
                    (Default: ['#', 'X'])
        
    """
    Bins_Pos = {}; Bins_Sum = {}
    genome = Load_Genome(genomeSize, chroms)
    Whole_Lib = {}
    Local_Lib = {}
    for res in wholeRes:
        Bins, Sum = Get_Chro_Bins(genomeSize, res, chroms)
        Bins_Pos[res] = Bins
        Bins_Sum[res] = Sum
        Whole_Lib[res] = {}
        Whole_Lib[res]['Bins'] = Bins
        Whole_Lib[res]['Matrix'] = np.zeros((Sum,Sum), dtype = np.int)
    
    for res in localRes:
        Local_Lib[res] = {}
        for c,l in genome.items():
            bin_size = (l // res) + 1
            Local_Lib[res][c] = np.zeros((bin_size,bin_size), dtype = np.int)
    
    count = 0
    for line in bed_IO:
        count += 1
        if count % 20000000 == 0:
            print "    %d pairs completed" % count
            log.log(21, "    %d pairs completed" % count)
        line = line.strip().split()
        
        c1 = line[1].lstrip('chr')
        c2 = line[8].lstrip('chr')
        check1 = ((not chroms) or (c1.isdigit() and ('#' in chroms)) or (c1 in chroms))
        check2 = ((not chroms) or (c2.isdigit() and ('#' in chroms)) or (c2 in chroms))

        if check1 and check2:
            ##WholeGenome Matrix Building
            for res in wholeRes:
                s1 = int(line[6]) // res
                bin1 = s1 + Bins_Pos[res][c1][0]
                s2 = int(line[13]) // res
                bin2 = s2 + Bins_Pos[res][c2][0]
                
                if bin1 != bin2:
                    Whole_Lib[res]['Matrix'][bin1][bin2] += 1
                    Whole_Lib[res]['Matrix'][bin2][bin1] += 1
                else:
                    Whole_Lib[res]['Matrix'][bin1][bin2] += 1
            
            ##Intra-Chromosome Matrix Building
            if c1 == c2:
                for res in localRes:
                    pos1 = int(line[6]) // res
                    pos2 = int(line[13]) // res
                    if pos1 != pos2:
                        Local_Lib[res][c1][pos1][pos2] += 1
                        Local_Lib[res][c1][pos2][pos1] += 1
                    else:
                        Local_Lib[res][c1][pos1][pos2] += 1
    
    for res, lib in Whole_Lib.items():
        Matrix = lib['Matrix']
        Bins = lib['Bins']
        Whole_Lib[res] = WholeMatrixToSparseDict(Bins, Matrix)
    
    for res, lib in Local_Lib.items():
        Local_Lib[res] = IntraMatrixToSparseDict(lib)
    
    return Whole_Lib, Local_Lib



def TraditionalMatrixConstruction(OutPath, RepPath, genomeSize, wholeRes, 
                                   localRes, chroms = ['#', 'X'], balance = True):
                                       
    """
    """
    print "Building Replicate Matrix respectively"
    log.log(21, "Building Replicate Matrix respectively")
    
    CoolerPath = os.path.join(OutPath,'Cooler')
    if not os.path.exists(CoolerPath):
        os.mkdir(CoolerPath)
    
    Replicate_Cooler = []
    
    for rep_p in RepPath:
        files = [i for i in os.listdir(rep_p) if '_Valid.bed' in i]
        prefix = files[0].split('Valid')[0]
        files = [os.path.join(rep_p,fil) for fil in files]
        
        Merge_cmd = Merge_beds(files)
        pread = subprocess.Popen(Merge_cmd, stdout=subprocess.PIPE, bufsize=-1)        
        instream = pread.stdout
        
        Whole_Lib, Local_Lib = TraditionalMatrixBuilding(bed_IO = instream,
                                                          genomeSize = genomeSize,
                                                          wholeRes = wholeRes,
                                                          localRes = localRes,
                                                          chroms = chroms)
        
#        if saveNPZ:
#            ## NPZ Path, Save the data into Numpy NPZ files
#            NPZPath = os.path.join(OutPath,'NPZ')
#            if not os.path.exists(NPZPath):
#                os.mkdir(NPZPath)
#            outfil_whole = os.path.join(NPZPath,prefix+'Whole_Genome_Matrix.npz')
#            outfil_local = os.path.join(NPZPath,prefix+'Intra_Chromosome_Matrix.npz')
#            
#            np.savez_compressed(outfil_whole, **Whole_Lib)
#            np.savez_compressed(outfil_local, **Local_Lib)
        
        
        ## Cooler Path, Save the data into cooler files
        
        instream.close()
        Multi_cooler = os.path.join(CoolerPath,prefix+'Multi.cool')
        
        NPZ2Cooler(datasets = Whole_Lib, 
                   outfil = Multi_cooler,
                   genomeSizes_file = genomeSize,
                   chroms = chroms,
                   onlyIntra = False,
                   dtype = 'int')
        
        NPZ2Cooler(datasets = Local_Lib,
                   outfil = Multi_cooler,
                   genomeSizes_file = genomeSize,
                   chroms = chroms,
                   onlyIntra = True,
                   dtype = 'int')
        
        Replicate_Cooler.append(Multi_cooler)
        
        print "    %s finished" % Multi_cooler
        log.log(21, "    %s finished" % Multi_cooler)
        
    
    print "Merging the replicates ..."
    log.log(21, "Merging the replicates ...")
    
    Merged_cooler = os.path.join(CoolerPath,'Merged_Multi.cool')
    mergebuf = int(20e6)
    
    for res in wholeRes+localRes:
        merge_cooler = '{}::{}'.format(Merged_cooler, res)
        replicate_cooler = ['{}::{}'.format(i,res) for i in Replicate_Cooler]
        cooler.merge_coolers(output_uri = merge_cooler,
                             input_uris = replicate_cooler,
                             mergebuf = mergebuf)
    print "    %s finished" % Merged_cooler
    log.log(21, "    %s finished" % Merged_cooler)
    
    ## ICE balance with cooler
    if balance:
        
        print "    Balancing start ..."
        log.log(21, "    Balancing start ...")
        coolers = Replicate_Cooler + [Merged_cooler]
        
        for res in wholeRes:
            for cooler_file in coolers:
                cmd = 'cooler balance --ignore-diags 1 --force {}::{}'.format(cooler_file, res)
                subprocess.call(cmd, shell = True)
        
        for res in localRes:
            for cooler_file in coolers:
                cmd = 'cooler balance --ignore-diags 1 --cis-only --force {}::{}'.format(cooler_file, res)
                subprocess.call(cmd, shell = True)
    
    print "All Done!"
    log.log(21, "    All Done!")
    


def GetNeighborhoodIndex(L):
    """
    """
    center = L+1 # center index
    i_index=[]
    j_index=[]
    for i in range(L*2+1):
        for j in range(L*2+1):
            if math.sqrt((i-center)**2 + (j-center)**2) < math.sqrt(L):
                i_index.append(i)
                j_index.append(j)
    return i_index, j_index
    
def GetNeighborhoodContacts(M, i_index, j_index):
    """
        Subset Contacts for imputation.
    """
    return M[i_index,j_index]
    


def Gap_definedLowRes(Matrix):
    """
    """
    threshold = 0.1
    
    gap = []
    for index,i in enumerate(Matrix):
        tmp = 1 - (i == 0).sum() / float(len(i))
        if tmp < threshold:
            gap.append(index)
                        
    return np.array(gap)


def Non_Gap_DefinedLowRes(N,gap):
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
    

def Trans2symmetryLowRes(Matrix):
    """
    """
    upper_M = np.triu(Matrix) + np.tril(Matrix,-1).T
    Sym_M = np.triu(upper_M,1).T + upper_M
        
    return Sym_M
        


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


def TraditionalMatrixInAllelic(bed_IO, genomeSize, wholeRes, localRes, chroms):
    """
         Build Tranditional Matrix in Haplotype-resolved Hi-C Pipeline.
    """

    Bins_Pos = {}; Bins_Sum = {}
    genome = Load_Genome(genomeSize, chroms)
    Whole_Lib = {}
    Local_Lib = {}
    for res in wholeRes:
        Bins, Sum = Get_Chro_Bins(genomeSize, res, chroms)
        Bins_Pos[res] = Bins
        Bins_Sum[res] = Sum
        Whole_Lib[res] = {}
        Whole_Lib[res]['Bins'] = Bins
        Whole_Lib[res]['Matrix'] = np.zeros((Sum,Sum), dtype = np.int)
    
    for res in localRes:
        Local_Lib[res] = {}
        for c,l in genome.items():
            bin_size = (l // res) + 1
            Local_Lib[res][c] = np.zeros((bin_size,bin_size), dtype = np.int)
    
    count = 0
    for line in bed_IO:
        count += 1
        if count % 20000000 == 0:
            print "    %d pairs completed" % count
            log.log(21, "    %d pairs completed" % count)
        line = line.strip().split()
        
        c1 = line[0].lstrip('chr')
        c2 = line[2].lstrip('chr')
        check1 = ((not chroms) or (c1.isdigit() and ('#' in chroms)) or (c1 in chroms))
        check2 = ((not chroms) or (c2.isdigit() and ('#' in chroms)) or (c2 in chroms))

        if check1 and check2:
            ##WholeGenome Matrix Building
            for res in wholeRes:
                s1 = int(line[1]) // res
                bin1 = s1 + Bins_Pos[res][c1][0]
                s2 = int(line[3]) // res
                bin2 = s2 + Bins_Pos[res][c2][0]
                
                if bin1 != bin2:
                    Whole_Lib[res]['Matrix'][bin1][bin2] += 1
                    Whole_Lib[res]['Matrix'][bin2][bin1] += 1
                else:
                    Whole_Lib[res]['Matrix'][bin1][bin2] += 1
            
            ##Intra-Chromosome Matrix Building
            if c1 == c2:
                for res in localRes:
                    pos1 = int(line[1]) // res
                    pos2 = int(line[3]) // res
                    if pos1 != pos2:
                        Local_Lib[res][c1][pos1][pos2] += 1
                        Local_Lib[res][c1][pos2][pos1] += 1
                    else:
                        Local_Lib[res][c1][pos1][pos2] += 1
    
    return Whole_Lib, Local_Lib


def GenomeWideMatrixCorrection(Bins_Pos, Hap_Bins_Pos, T_M, H_M):
    """
    """
    Beta = {}
    for chro in Bins_Pos.keys():
        if chro not in Beta.keys():
            Beta[chro] = []
            Tra_M = T_M[Bins_Pos[chro][0]:Bins_Pos[chro][1]+1,
                        Bins_Pos[chro][0]:Bins_Pos[chro][1]+1]
            M_M = H_M[Hap_Bins_Pos['M'+chro][0]:Hap_Bins_Pos['M'+chro][1]+1,
                      Hap_Bins_Pos['M'+chro][0]:Hap_Bins_Pos['M'+chro][1]+1]
            
            P_P = H_M[Hap_Bins_Pos['P'+chro][0]:Hap_Bins_Pos['P'+chro][1]+1,
                      Hap_Bins_Pos['P'+chro][0]:Hap_Bins_Pos['P'+chro][1]+1]
            
            Gap = Gap_definedLowRes(Tra_M)
            N = Tra_M.shape[0]
            Non_Gap_M = Non_Gap_DefinedLowRes(N, Gap)
            
            alpha = []
            
            for i in range(N):
                alpha.append((M_M[i].sum()+P_P[i].sum()) / (Tra_M[i].sum() + 1))
            
            alpha = np.array(alpha)
            alpha /= np.max(alpha[Non_Gap_M])
            alpha[alpha == 0] = 1
            threshold = np.percentile(alpha[Non_Gap_M], 20)
            alpha[alpha < threshold] = threshold
            Beta[chro] = alpha
    
    chros = Sort_Chromosomes(Beta.keys())
    Alpha = []
    for i in chros:
        Alpha.extend(Beta[i])
    Alpha += Alpha
    Alpha = np.array(Alpha)
    
    S_H_M = H_M / Alpha[:, None]
    Sym_S_H_M = Trans2symmetryLowRes(S_H_M)
    Cor_S_H_M = Correct_VC(Sym_S_H_M, 2/3)
    R_F = H_M.mean() / Cor_S_H_M.mean()
    Nor_S_H_M = R_F * Cor_S_H_M
    
    return Nor_S_H_M


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
    if threshold > 0.2:
        threshold = 0.2
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
        for j in gap:
            if i == j:
                New_Matrix[i][j] = Matrix[i][j]
            else:
                value = max(Matrix[i][j], Matrix[j][i])
                New_Matrix[i][j] = value
                New_Matrix[j][i] = value
        for j in Non_Gap:
            New_Matrix[i][j] = Matrix[j][i]
            New_Matrix[j][i] = Matrix[j][i]

    # NonGap <-> NonGap    
    for i in Non_Gap:
        for j in Non_Gap:
            if i == j:
                New_Matrix[i][j] = Matrix[i][j]
            else:
                value = (Matrix[i][j] + Matrix[j][i]) / 2.0
                New_Matrix[i][j] = value
                New_Matrix[j][i] = value
    
    return New_Matrix
        



def TwoStepCorrection(TM, MM, PM):
    """
    """
    alpha = []
    N = TM.shape[0]
    Gap_M = Gap_defined(MM)
    Gap_P = Gap_defined(PM)
    Non_Gap_M = Non_Gap_Defined(N, Gap_M)
    Non_Gap_P = Non_Gap_Defined(N, Gap_P)
    
    for i in range(N):
        alpha.append((MM[i].sum() + PM[i].sum()) / (TM[i].sum() + 1))
    
    alpha = np.array(alpha)
    # remove Gap effect
    Non_Gap_index = list(set(Non_Gap_M) | set(Non_Gap_P))
    
    # SNP Correcting
    alpha /= np.max(alpha[Non_Gap_index])
    alpha[alpha == 0] = 1
    threshold = np.percentile(alpha[Non_Gap_index], 20)
    alpha[alpha < threshold] = threshold
    
    S_MM = MM / alpha[:, None]
    S_PM = PM / alpha[:, None]
    
    Sym_MM = Trans2symmetry(S_MM, Gap_M)
    Sym_PM = Trans2symmetry(S_PM, Gap_P)
    
    # HiC Correcting
    Cor_MM = Correct_VC(Sym_MM, 2/3)
    Cor_PM = Correct_VC(Sym_PM, 2/3)
    
    M_RF = MM.mean() / Cor_MM.mean()
    P_RF = PM.mean() / Cor_PM.mean()
    
    Nor_MM = M_RF * Cor_MM
    Nor_PM = P_RF * Cor_PM
    
    return Nor_MM, Nor_PM, Gap_M, Gap_P
    
    
def IntraChromMatrixCorrection(Tra_Lib, Hap_Lib):
    """
    """
    Nor_Lib = {}
    Gap_Lib = {}
    for chro in Tra_Lib.keys():
        TM = Tra_Lib[chro]
        MM = Hap_Lib['M'+chro]
        PM = Hap_Lib['P'+chro]
        Nor_MM, Nor_PM, Gap_M, Gap_P = TwoStepCorrection(TM, MM, PM)
        Nor_Lib['M'+chro] = Nor_MM
        Nor_Lib['P'+chro] = Nor_PM
        Gap_Lib['M'+chro] = Gap_M
        Gap_Lib['P'+chro] = Gap_P
    
    return Nor_Lib, Gap_Lib


def HaplotypeMatrixBuilding(OutPath, BedPath, genomeSize, wholeRes, localRes, 
                            Imputation_region, Imputation_min, Imputation_ratio, chroms):
    """
        Build and Imputation All Haplotype Matrix. Include:
            Maternal -- Maternal
            Paternal -- Paternal
            Maternal -- Paternal
            Paternal -- Maternal
        
        Whole_Res : 
            Lower resolution
            We suggest the value may bigger that 500K.
                
    """
    
    DataSets = {}    
    
    files = [i for i in os.listdir(BedPath) if 'Bi_Allelic.bed' in i or \
                'M_M.bed' in i or 'M_P.bed' in i or 'P_P.bed' in i or 'P_M.bed' in i]
    
    files.sort()
    prefix = files[0].split('Valid')[0]
    
    print "Matrix Construction for %s " % prefix
    log.log(21, "Matrix Construction for %s " % prefix)
    
    if len(files) != 5:
        mark, tmp = Check_Bed(files)
        if mark == True:
            pass
        else:
            raise Exception ('Missing file %s.bed in %s' % (tmp, BedPath))
    
    files = [os.path.join(BedPath, fil) for fil in files]
    
    print "    Tranditional Matrix Building ..."
    log.log(21, "    Tranditional Matrix Building ...")
    Trans_fil = files
    Trans_cmd = Merge_beds(Trans_fil)
    pread = subprocess.Popen(Trans_cmd, stdout=subprocess.PIPE,bufsize=-1)
    instream = pread.stdout
    
    Whole_Lib, Local_Lib = TraditionalMatrixInAllelic(bed_IO = instream,
                                                       genomeSize = genomeSize,
                                                       wholeRes = wholeRes,
                                                       localRes = localRes,
                                                       chroms = chroms)
    instream.close()
    
    DataSets['Tradition_Whole'] = Whole_Lib
    DataSets['Tradition_Local'] = Local_Lib
      
    ## Chromosome Bins prepare for Haplotype-resolved Matrix.
    Hap_Bins_Pos = {}
    Hap_Bins_Sum = {}
    for res in wholeRes:
        Bins, Sum = Get_Chro_Bins_Haplotypes(genomeSize, res, chroms)
        Hap_Bins_Pos[res] = Bins
        Hap_Bins_Sum[res] = Sum
        
    ## UnImputated Matrix Building
    
    print "    UnImputated Matrix Building ..."
    log.log(21, "    UnImputated Matrix Building ...")
    
    UnImputated_Whole_Lib = {}
    UnImputated_Local_Lib = {}
    
    for res in wholeRes:
        UnImputated_Whole_Lib[res] = {}
        UnImputated_Whole_Lib[res]['Bins'] = Hap_Bins_Pos[res]
        UnImputated_Whole_Lib[res]['Matrix'] = np.zeros((Hap_Bins_Sum[res], Hap_Bins_Sum[res]), dtype = np.int)
    
    Hap_Genome = Load_HaplotypeGenome(genomeSize, chroms)
    
    for res in localRes:
        UnImputated_Local_Lib[res] = {}
        for chro, l in Hap_Genome.items():
            binSize = (l // res) + 1
            UnImputated_Local_Lib[res][chro] = np.zeros((binSize,binSize), dtype = np.int)
    
    
    ## =================Save Maternal-Maternal contact data ==================
    M_M_fil = [i for i in files if 'M_M.bed' in i]
    M_M_cmd = Merge_beds(M_M_fil)
    pread = subprocess.Popen(M_M_cmd, stdout=subprocess.PIPE, bufsize=-1)
    instream = pread.stdout
    for line in instream:
        line = line.strip().split()
        if line[-1] != 'Both':
            continue
        else:
            c1 = line[0].lstrip('chr')
            c2 = line[2].lstrip('chr')
            check1 = ((not chroms) or (c1.isdigit() and ('#' in chroms)) or (c1 in chroms))
            check2 = ((not chroms) or (c2.isdigit() and ('#' in chroms)) or (c2 in chroms))
            if not (check1 and check2):
                continue
            c1 = 'M'+c1
            c2 = 'M'+c2
            for res in wholeRes:
                bin1 = int(line[1]) // res + Hap_Bins_Pos[res][c1][0]
                bin2 = int(line[3]) // res + Hap_Bins_Pos[res][c2][0]
                if bin1 == bin2:
                    UnImputated_Whole_Lib[res]['Matrix'][bin1][bin2] += 1
                else:
                    UnImputated_Whole_Lib[res]['Matrix'][bin1][bin2] += 1
                    UnImputated_Whole_Lib[res]['Matrix'][bin2][bin1] += 1
            
            if c1 == c2:
                for res in localRes:
                    bin1 = int(line[1]) // res 
                    bin2 = int(line[3]) // res 
                    if bin1 == bin2:
                        UnImputated_Local_Lib[res][c1][bin1][bin2] += 1
                    else:
                        UnImputated_Local_Lib[res][c1][bin1][bin2] += 1
                        UnImputated_Local_Lib[res][c1][bin2][bin1] += 1
    instream.close()
    
    ## =================Save Paternal-Paternal contact data ==================
    P_P_fil = [i for i in files if 'P_P.bed' in i]
    P_P_cmd = Merge_beds(P_P_fil)
    pread = subprocess.Popen(P_P_cmd, stdout=subprocess.PIPE, bufsize=-1)
    instream = pread.stdout
    for line in instream:
        line = line.strip().split()
        if line[-1] != 'Both':
            continue
        else:
            c1 = line[0].lstrip('chr')
            c2 = line[2].lstrip('chr')
            check1 = ((not chroms) or (c1.isdigit() and ('#' in chroms)) or (c1 in chroms))
            check2 = ((not chroms) or (c2.isdigit() and ('#' in chroms)) or (c2 in chroms))
            if not (check1 and check2):
                continue
            c1 = 'P'+c1
            c2 = 'P'+c2
            for res in wholeRes:
                bin1 = int(line[1]) // res + Hap_Bins_Pos[res][c1][0]
                bin2 = int(line[3]) // res + Hap_Bins_Pos[res][c2][0]
                if bin1 == bin2:
                    UnImputated_Whole_Lib[res]['Matrix'][bin1][bin2] += 1
                else:
                    UnImputated_Whole_Lib[res]['Matrix'][bin1][bin2] += 1
                    UnImputated_Whole_Lib[res]['Matrix'][bin2][bin1] += 1
            
            if c1 == c2:
                for res in localRes:
                    bin1 = int(line[1]) // res 
                    bin2 = int(line[3]) // res 
                    if bin1 == bin2:
                        UnImputated_Local_Lib[res][c1][bin1][bin2] += 1
                    else:
                        UnImputated_Local_Lib[res][c1][bin1][bin2] += 1
                        UnImputated_Local_Lib[res][c1][bin2][bin1] += 1
    instream.close()    
    
    ## =================Save Maternal- Paternal contact data ==================
    M_P_fil = [i for i in files if 'M_P.bed' in i]
    M_P_cmd = Merge_beds(M_P_fil)
    pread = subprocess.Popen(M_P_cmd, stdout=subprocess.PIPE, bufsize=-1)
    instream = pread.stdout
    for line in instream:
        line = line.strip().split()
        c1 = line[0].lstrip('chr')
        c2 = line[2].lstrip('chr')
        check1 = ((not chroms) or (c1.isdigit() and ('#' in chroms)) or (c1 in chroms))
        check2 = ((not chroms) or (c2.isdigit() and ('#' in chroms)) or (c2 in chroms))
        if not (check1 and check2):
            continue
        c1 = 'M'+c1
        c2 = 'P'+c2
        for res in wholeRes:
            bin1 = int(line[1]) // res + Hap_Bins_Pos[res][c1][0]
            bin2 = int(line[3]) // res + Hap_Bins_Pos[res][c2][0]
            UnImputated_Whole_Lib[res]['Matrix'][bin1][bin2] += 1
            UnImputated_Whole_Lib[res]['Matrix'][bin2][bin1] += 1
    
    instream.close()
    
    P_M_fil = [i for i in files if 'P_M.bed' in i]
    P_M_cmd = Merge_beds(P_M_fil)
    pread = subprocess.Popen(P_M_cmd, stdout=subprocess.PIPE, bufsize=-1)
    instream = pread.stdout
    for line in instream:
        line = line.strip().split()
        c1 = line[0].lstrip('chr')
        c2 = line[2].lstrip('chr')
        check1 = ((not chroms) or (c1.isdigit() and ('#' in chroms)) or (c1 in chroms))
        check2 = ((not chroms) or (c2.isdigit() and ('#' in chroms)) or (c2 in chroms))
        if not (check1 and check2):
            continue
        c1 = 'P'+c1
        c2 = 'M'+c2
        for res in wholeRes:
            bin1 = int(line[1]) // res + Hap_Bins_Pos[res][c1][0]
            bin2 = int(line[3]) // res + Hap_Bins_Pos[res][c2][0]
            UnImputated_Whole_Lib[res]['Matrix'][bin1][bin2] += 1
            UnImputated_Whole_Lib[res]['Matrix'][bin2][bin1] += 1
    
    instream.close()
    
    DataSets['UnImputated_Whole'] = UnImputated_Whole_Lib
    DataSets['UnImputated_Local'] = UnImputated_Local_Lib
    
    
    ## Data Imputation
    print "    Data Imputation ..."
    log.log(21, "    Data Imputation ...")
    
    Imputated_Whole_Lib = copy.deepcopy(UnImputated_Whole_Lib)
    Imputated_Local_Lib = copy.deepcopy(UnImputated_Local_Lib)
    
    
    ## Imputation parameters
    sub_Index = {}
    i_Index = {}
    j_Index = {}
    for res in wholeRes:
        sub_Index[res] = Imputation_region // res
        i_Index[res], j_Index[res] = GetNeighborhoodIndex(sub_Index[res])


    print "    Imputation for M -- ? and ? -- M Contacts ..."
    log.log(21, "    Imputation for M -- ? and ? -- M Contacts ...")
    pread = subprocess.Popen(M_M_cmd, stdout=subprocess.PIPE,bufsize=-1)
    instream = pread.stdout
    for line in instream:
        line = line.strip().split()
        mark = line[-1]
        if mark == 'Both':
            continue
        else:
            c1 = line[0].lstrip('chr')
            c2 = line[2].lstrip('chr')
            check1 = ((not chroms) or (c1.isdigit() and ('#' in chroms)) or (c1 in chroms))
            check2 = ((not chroms) or (c2.isdigit() and ('#' in chroms)) or (c2 in chroms))
            if not(check1 and check2):
                continue
            if c1 == c2:
                for res in wholeRes:
                    pos1 = int(line[1]) // res
                    pos2 = int(line[3]) // res
                    bin1 = pos1 + Hap_Bins_Pos[res]['M'+c1][0]
                    bin2 = pos2 + Hap_Bins_Pos[res]['M'+c2][0]
                    if mark == 'R1':
                        Imputated_Whole_Lib[res]['Matrix'][bin1][bin2] += 1
                    else:
                        Imputated_Whole_Lib[res]['Matrix'][bin2][bin1] += 1
                
                for res in localRes:
                    pos1 = int(line[1]) // res
                    pos2 = int(line[3]) // res
                    if mark == 'R1':
                        Imputated_Local_Lib[res]['M'+c1][pos1][pos2] += 1
                    else:
                        Imputated_Local_Lib[res]['M'+c2][pos2][pos1] += 1
            else:
                ## Imputation for Inter-chromosome contact.
                for res in wholeRes:
                    pos1 = int(line[1]) // res 
                    pos2 = int(line[3]) // res
                    if mark == 'R1':
                        bin1 = pos1 + Hap_Bins_Pos[res]['M'+c1][0]
                        M_bin2 = pos2 + Hap_Bins_Pos[res]['M'+c2][0]
                        P_bin2 = pos2 + Hap_Bins_Pos[res]['P'+c2][0]
                        
                        s_i = sub_Index[res]
                        Matrix = UnImputated_Whole_Lib[res]['Matrix']
                        
                        if bin1 < s_i or M_bin2 < s_i or P_bin2 < s_i:
                            continue
                        if bin1 + s_i + 1 > Matrix.shape[0] or \
                            M_bin2 + s_i + 1 > Matrix.shape[0] or \
                            P_bin2 + s_i + 1 > Matrix.shape[0]:
                                continue
                        
                        i_i = i_Index[res]
                        j_i = j_Index[res]
                        
                        M_M_sub = Matrix[bin1-s_i:bin1+s_i+1,
                                         M_bin2-s_i:M_bin2+s_i+1]
                        
                        M_P_sub = Matrix[bin1-s_i:bin1+s_i+1,
                                         P_bin2-s_i:P_bin2+s_i+1]
                        
                        M_M_subset = GetNeighborhoodContacts(M_M_sub, i_i, j_i)                        
                        M_P_subset = GetNeighborhoodContacts(M_P_sub, i_i, j_i)
                        
                        M_M_sum = M_M_subset.sum()
                        M_P_sum = M_P_subset.sum()
                        
                        if M_M_subset.sum() >= Imputation_min and (M_M_sum / (M_M_sum+M_P_sum)) > Imputation_ratio:
                            Imputated_Whole_Lib[res]['Matrix'][bin1][M_bin2] += 1
                        elif M_P_subset.sum() >= Imputation_min and (M_P_sum / (M_M_sum+M_P_sum)) > Imputation_ratio:
                            Imputated_Whole_Lib[res]['Matrix'][bin1][P_bin2] += 1
                        else:
                            continue
                    else:
                        bin2 = pos2 + Hap_Bins_Pos[res]['M'+c1][0]
                        M_bin1 = pos1 + Hap_Bins_Pos[res]['M'+c2][0]
                        P_bin1 = pos1 + Hap_Bins_Pos[res]['P'+c2][0]
                        
                        s_i = sub_Index[res]
                        Matrix = UnImputated_Whole_Lib[res]['Matrix']
                        
                        if bin2 < s_i or M_bin1 < s_i or P_bin1 < s_i:
                            continue
                        if bin2 + s_i + 1 > Matrix.shape[0] or \
                            M_bin1 + s_i + 1 > Matrix.shape[0] or \
                            P_bin1 + s_i + 1 > Matrix.shape[0]:
                                continue
                        
                        i_i = i_Index[res]
                        j_i = j_Index[res]
                        
                        M_M_sub = Matrix[M_bin1-s_i:M_bin1+s_i+1,
                                         bin2-s_i:bin2+s_i+1]
                        
                        M_P_sub = Matrix[P_bin1-s_i:P_bin1+s_i+1,
                                         bin2-s_i:bin2+s_i+1]
                        
                        M_M_subset = GetNeighborhoodContacts(M_M_sub, i_i, j_i)                        
                        M_P_subset = GetNeighborhoodContacts(M_P_sub, i_i, j_i)
                        
                        M_M_sum = M_M_subset.sum()
                        M_P_sum = M_P_subset.sum()
                        
                        if M_M_sum >= Imputation_min and (M_M_sum / (M_M_sum+M_P_sum)) > Imputation_ratio:
                            Imputated_Whole_Lib[res]['Matrix'][bin2][M_bin1] += 1
                        elif M_P_sum >= Imputation_min and (M_P_sum / (M_M_sum+M_P_sum)) > Imputation_ratio:
                            Imputated_Whole_Lib[res]['Matrix'][bin2][P_bin1] += 1
                        else:
                            continue
    
    instream.close()

    print "    Imputation for P -- ? and ? -- P Contacts ..."
    log.log(21, "    Imputation for P -- ? and ? -- P Contacts ...")
    pread = subprocess.Popen(P_P_cmd, stdout=subprocess.PIPE,bufsize=-1)
    instream = pread.stdout
    for line in instream:
        line = line.strip().split()
        mark = line[-1]
        if mark == 'Both':
            continue
        else:
            c1 = line[0].lstrip('chr')
            c2 = line[2].lstrip('chr')
            check1 = ((not chroms) or (c1.isdigit() and ('#' in chroms)) or (c1 in chroms))
            check2 = ((not chroms) or (c2.isdigit() and ('#' in chroms)) or (c2 in chroms))
            if not(check1 and check2):
                continue
            if c1 == c2:
                for res in wholeRes:
                    pos1 = int(line[1]) // res
                    pos2 = int(line[3]) // res
                    bin1 = pos1 + Hap_Bins_Pos[res]['P'+c1][0]
                    bin2 = pos2 + Hap_Bins_Pos[res]['P'+c2][0]
                    if mark == 'R1':
                        Imputated_Whole_Lib[res]['Matrix'][bin1][bin2] += 1
                    else:
                        Imputated_Whole_Lib[res]['Matrix'][bin2][bin1] += 1
                
                for res in localRes:
                    pos1 = int(line[1]) // res
                    pos2 = int(line[3]) // res
                    if mark == 'R1':
                        Imputated_Local_Lib[res]['P'+c1][pos1][pos2] += 1
                    else:
                        Imputated_Local_Lib[res]['P'+c2][pos2][pos1] += 1
            else:
                ## Imputation for Inter-chromosome contact.
                for res in wholeRes:
                    pos1 = int(line[1]) // res 
                    pos2 = int(line[3]) // res
                    if mark == 'R1':
                        bin1 = pos1 + Hap_Bins_Pos[res]['P'+c1][0]
                        M_bin2 = pos2 + Hap_Bins_Pos[res]['M'+c2][0]
                        P_bin2 = pos2 + Hap_Bins_Pos[res]['P'+c2][0]
                        
                        s_i = sub_Index[res]
                        Matrix = UnImputated_Whole_Lib[res]['Matrix']
                        
                        if bin1 < s_i or M_bin2 < s_i or P_bin2 < s_i:
                            continue
                        if bin1 + s_i + 1 > Matrix.shape[0] or \
                            M_bin2 + s_i + 1 > Matrix.shape[0] or \
                            P_bin2 + s_i + 1 > Matrix.shape[0]:
                                continue
                        
                        i_i = i_Index[res]
                        j_i = j_Index[res]
                        
                        P_P_sub = Matrix[bin1-s_i:bin1+s_i+1,
                                         P_bin2-s_i:P_bin2+s_i+1]
                        
                        M_P_sub = Matrix[bin1-s_i:bin1+s_i+1,
                                         P_bin2-s_i:P_bin2+s_i+1]
                        
                        P_P_subset = GetNeighborhoodContacts(M_M_sub, i_i, j_i)                        
                        M_P_subset = GetNeighborhoodContacts(M_P_sub, i_i, j_i)
                        
                        P_P_sum = P_P_subset.sum()
                        M_P_sum = M_P_subset.sum()
                        
                        if P_P_sum >= Imputation_min and (P_P_sum / (P_P_sum+M_P_sum)) > Imputation_ratio:
                            Imputated_Whole_Lib[res]['Matrix'][bin1][M_bin2] += 1
                        elif M_P_sum >= Imputation_min and (M_P_sum / (P_P_sum+M_P_sum)) > Imputation_ratio:
                            Imputated_Whole_Lib[res]['Matrix'][bin1][P_bin2] += 1
                        else:
                            continue
                    else:
                        bin2 = pos2 + Hap_Bins_Pos[res]['P'+c1][0]
                        M_bin1 = pos1 + Hap_Bins_Pos[res]['M'+c2][0]
                        P_bin1 = pos1 + Hap_Bins_Pos[res]['P'+c2][0]
                        
                        s_i = sub_Index[res]
                        Matrix = UnImputated_Whole_Lib[res]['Matrix']
                        
                        if bin2 < s_i or M_bin1 < s_i or P_bin1 < s_i:
                            continue
                        if bin2 + s_i + 1 > Matrix.shape[0] or \
                            M_bin1 + s_i + 1 > Matrix.shape[0] or \
                            P_bin1 + s_i + 1 > Matrix.shape[0]:
                                continue
                        
                        i_i = i_Index[res]
                        j_i = j_Index[res]
                        
                        P_P_sub = Matrix[P_bin1-s_i:P_bin1+s_i+1,
                                         bin2-s_i:bin2+s_i+1]
                        
                        M_P_sub = Matrix[M_bin1-s_i:M_bin1+s_i+1,
                                         bin2-s_i:bin2+s_i+1]
                        
                        P_P_subset = GetNeighborhoodContacts(P_P_sub, i_i, j_i)                        
                        M_P_subset = GetNeighborhoodContacts(M_P_sub, i_i, j_i)
                        
                        P_P_sum = P_P_subset.sum()
                        M_P_sum = M_P_subset.sum()
                        
                        if P_P_sum >= Imputation_min and (P_P_sum / (P_P_sum+M_P_sum)) > Imputation_ratio:
                            Imputated_Whole_Lib[res]['Matrix'][P_bin1][bin2] += 1
                        elif M_P_subset.sum() >= Imputation_min and (M_P_sum / (P_P_sum+M_P_sum)) > Imputation_ratio:
                            Imputated_Whole_Lib[res]['Matrix'][M_bin1][bin2] += 1
                        else:
                            continue
    
    instream.close()                   
        
    print "    Imputation Done"
    log.log(21, "    Imputation Done")
    
    DataSets['Imputated_Whole'] = Imputated_Whole_Lib
    DataSets['Imputated_Local'] = Imputated_Local_Lib
                
    print "    Matrix Correcting and save to cooler ..."
    log.log(21, "    Matrix Correcting and save to cooler ...")
    print "    Traditional Matrix starting ..."
    log.log(21, "    Traditional Matrix starting ...")
    Whole_Lib = {}
    Local_Lib = {}
    for res in wholeRes:
        Matrix = DataSets['Tradition_Whole'][res]['Matrix']
        Bins = DataSets['Tradition_Whole'][res]['Bins']
        Whole_Lib[res] = WholeMatrixToSparseDict(Bins, Matrix)
    
    for res in localRes:
        Local_Lib[res] = IntraMatrixToSparseDict(DataSets['Tradition_Local'][res])
    
    Tradition_cooler = os.path.join(OutPath,prefix+'Traditional_Multi.cool')
    
    NPZ2Cooler(datasets = Whole_Lib,
               outfil = Tradition_cooler,
               genomeSizes_file = genomeSize,
               chroms = chroms,
               onlyIntra = False,
               dtype = 'int')
    
    NPZ2Cooler(datasets = Local_Lib,
               outfil = Tradition_cooler,
               genomeSizes_file = genomeSize,
               chroms = chroms,
               onlyIntra = True,
               dtype = 'int')
    
    
    print "    ICE Balance for Traditional Matrix ..."
    log.log(21, "    ICE Balance for Traditional Matrix ...")
    
    for res in wholeRes:
        cmd = 'cooler balance --ignore-diags 1 --force {}::{}'.format(Tradition_cooler, res)
        a = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        a.communicate()
    
    for res in localRes:
        cmd = 'cooler balance --ignore-diags 1 --cis-only --force {}::{}'.format(Tradition_cooler, res)
        a = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        a.communicate()
        
    print "    UnImputated Matrix start ..."
    log.log(21, "    UnImputated Matrix start ...")
    Whole_Lib = {}
    Local_Lib = {}
    
    Hap_genomeSize = os.path.join(OutPath,'Hap_genomeSize')
    out = open(Hap_genomeSize,'w')
    Hap_chroms = []
    with open(genomeSize,'r') as f:
        for line in f:
            line = line.strip().split()
            line[0] = line[0].lstrip('chr')
            c = line[0]
            if ((not chroms) or (c.isdigit() and ('#' in chroms)) or (c in chroms)):
                Hap_chroms.append('M'+c)
                Hap_chroms.append('P'+c)
                out.writelines('M'+'\t'.join(line)+'\n')
                out.writelines('P'+'\t'.join(line)+'\n')
    out.close()
    
    for res in wholeRes:
        Matrix = DataSets['UnImputated_Whole'][res]['Matrix']
        Bins = DataSets['UnImputated_Whole'][res]['Bins']
        Whole_Lib[res] = WholeMatrixToSparseDict(Bins, Matrix)
    
    for res in localRes:
        Local_Lib[res] = IntraMatrixToSparseDict(DataSets['UnImputated_Local'][res])
    
    UnImputated_cooler = os.path.join(OutPath,prefix+'UnImputated_Haplotype_Multi.cool')
    
    NPZ2Cooler(datasets = Whole_Lib,
               outfil = UnImputated_cooler,
               genomeSizes_file = Hap_genomeSize,
               chroms = Hap_chroms,
               onlyIntra = False,
               dtype = 'int')
               
    NPZ2Cooler(datasets = Local_Lib,
               outfil = UnImputated_cooler,
               genomeSizes_file = Hap_genomeSize,
               chroms = Hap_chroms,
               onlyIntra = True,
               dtype = 'int')
    
    print "    Imputated Matrix start ..."
    log.log(21, "    Imputated Matrix start ...")
    print "    Two-Step balancing for Imputated Matrix ..."
    log.log(21, "    Two-Step balancing for Imputated Matrix ...")
    
    Balanced_Whole_Data = {}
    for res in wholeRes:
        Bins = DataSets['Tradition_Whole'][res]['Bins']
        Hap_Bins = DataSets['Imputated_Whole'][res]['Bins']
        T_M = DataSets['Tradition_Whole'][res]['Matrix']
        H_M = DataSets['Imputated_Whole'][res]['Matrix']
        Balanced_M = GenomeWideMatrixCorrection(Bins_Pos = Bins,
                                                Hap_Bins_Pos = Hap_Bins,
                                                T_M = T_M,
                                                H_M = H_M)        
        Balanced_Whole_Data[res] = WholeMatrixToSparseDict(Hap_Bins, Balanced_M)
    
    Balanced_Local_Data = {}
    Gap_Local = {}
    for res in localRes:
        Tra_Lib = DataSets['Tradition_Local'][res]
        Hap_Lib = DataSets['Imputated_Local'][res]
        Nor_Lib, Gap_Lib = IntraChromMatrixCorrection(Tra_Lib, Hap_Lib)
        Balanced_Local_Data[res] = IntraMatrixToSparseDict(Nor_Lib)
        Gap_Local[str(res)] = Gap_Lib
    
    Gap_fil = os.path.join(OutPath,prefix+'Imputated_Gap.npz')
    np.savez(Gap_fil, **Gap_Local)
    
    Imputated_cooler = os.path.join(OutPath,prefix+'Imputated_Haplotype_Multi.cool')
    
    NPZ2Cooler(datasets = Balanced_Whole_Data,
               outfil = Imputated_cooler,
               genomeSizes_file = Hap_genomeSize,
               chroms = Hap_chroms,
               onlyIntra = False,
               dtype = 'float')
    
    NPZ2Cooler(datasets = Balanced_Local_Data,
               outfil = Imputated_cooler,
               genomeSizes_file = Hap_genomeSize,
               chroms = Hap_chroms,
               onlyIntra = True,
               dtype = 'float')
    
    print "    Done !!"
    log.log(21, "    Done !!")
    
    return prefix, DataSets

    
def HaplotypeMatrixConstruction(OutPath, RepPath, genomeSize, wholeRes, localRes,
                                Imputation_region = 10000000, Imputation_min = 2,
                                Imputation_ratio = 0.9, chroms = ['#', 'X']):
    """
        Haplotype-Resolved Matrix Construction.
        Finally,
            We support the three type of Haplotype-resolved Contact Matrix.
            
            First, Traditonal Matrix

            Second, UnImputated Haplotype-resolved Matrix which is build by contact
            pairs that both mate can be indicate to Maternal or Paternal
            
            Third, Imputated Haplotype-resolved Matrix which is build by Contact pairs
            that One mate can be indicate to Maternal or Paternal. The Imputation Method
            is similar to 
            <Three-dimensional genome structures of single diploid human cells>
                --Science.
        
        All Matrix are saved to cooler format...
        
        Haplotype-resolved Matrix do not have balance weight. The Count is already
        corrected by HiCHap.
                
    """
    print "Building Replicate Matrix respectively !!!"
    log.log(21, "Building Replicate Matrix respectively !!!")
    
    CoolerPath = os.path.join(OutPath, 'Cooler')
    if not os.path.exists(CoolerPath):
        os.mkdir(CoolerPath)
    
    
    All_Data = {}
    if len(RepPath) == 1:
        prefix, DataSets = HaplotypeMatrixBuilding(OutPath = CoolerPath,
                                                   BedPath = RepPath[0],
                                                   wholeRes = wholeRes,
                                                   localRes = localRes,
                                                   Imputation_ratio = Imputation_ratio,
                                                   Imputation_region = Imputation_region,
                                                   Imputation_min = Imputation_min,
                                                   chroms = chroms)
        print "    No replicates to merge."
        log.log(21, "    No replicates to merge.")
        print "All Done !!!"
        log.log(21, "All Done !!!")
        return
        
    for rep_p in RepPath:
        prefix, DataSets = HaplotypeMatrixBuilding(OutPath = CoolerPath,
                                                   BedPath = rep_p,
                                                   genomeSize = genomeSize,
                                                   wholeRes = wholeRes,
                                                   localRes = localRes,
                                                   Imputation_ratio = Imputation_ratio,
                                                   Imputation_region = Imputation_region,
                                                   Imputation_min = Imputation_min,
                                                   chroms = chroms)
        if len(All_Data) == 0:
            All_Data = copy.deepcopy(DataSets)
        else:
            for res in wholeRes:
                All_Data['Tradition_Whole'][res]['Matrix'] += DataSets['Tradition_Whole'][res]['Matrix']
                All_Data['UnImputated_Whole'][res]['Matrix'] += DataSets['UnImputated_Whole'][res]['Matrix']
                All_Data['Imputated_Whole'][res]['Matrix'] += DataSets['Imputated_Whole'][res]['Matrix']
            for res in localRes:
                
                k1 = All_Data['Tradition_Local'][res].keys()
                for k in k1:
                    All_Data['Tradition_Local'][res][k] += DataSets['Tradition_Local'][res][k]
                    
                k2 = All_Data['UnImputated_Local'][res].keys()
                for k in k2:
                    All_Data['UnImputated_Local'][res][k] += DataSets['UnImputated_Local'][res][k]
                    
                k3 = All_Data['Imputated_Local'][res].keys()
                for k in k3:
                    All_Data['Imputated_Local'][res][k] += DataSets['Imputated_Local'][res][k]
                    
    del DataSets
    gc.collect()
    time.sleep(10)
    
    
    print "Merging the Replicates..."
    log.log(21, "Merging the Replicates ...")
    print "    Matrix Correcting and save to cooler ..."
    log.log(21, "    Matrix Correcting and save to cooler ...")
    prefix = 'Merged_'
    Whole_Lib = {}
    Local_Lib = {}
    for res in wholeRes:
        Matrix = All_Data['Tradition_Whole'][res]['Matrix']
        Bins = All_Data['Tradition_Whole'][res]['Bins']
        Whole_Lib[res] = WholeMatrixToSparseDict(Bins, Matrix)
    
    for res in localRes:
        Local_Lib[res] = IntraMatrixToSparseDict(All_Data['Tradition_Local'][res])
    
    Tradition_cooler = os.path.join(CoolerPath,prefix+'Traditional_Multi.cool')
    
    NPZ2Cooler(datasets = Whole_Lib,
               outfil = Tradition_cooler,
               genomeSizes_file = genomeSize,
               chroms = chroms,
               onlyIntra = False,
               dtype = 'int')
    
    NPZ2Cooler(datasets = Local_Lib,
               outfil = Tradition_cooler,
               genomeSizes_file = genomeSize,
               chroms = chroms,
               onlyIntra = True,
               dtype = 'int')
    
    print "    ICE Balance for Traditional Matrix ..."
    log.log(21, "    ICE Balance for Traditional Matrix ...")
    
    for res in wholeRes:
        cmd = 'cooler balance --ignore-diags 1 --force {}::{}'.format(Tradition_cooler, res)
        a = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        a.communicate()
        
    for res in localRes:
        cmd = 'cooler balance --ignore-diags 1 --cis-only --force {}::{}'.format(Tradition_cooler, res)
        b = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        b.communicate()
    
    print "    UnImputated Matrix start ..."
    log.log(21, "    UnImputated Matrix start ...")
    Whole_Lib = {}
    Local_Lib = {}
   
    Hap_genomeSize = os.path.join(CoolerPath,'Hap_genomeSize')
    out = open(Hap_genomeSize,'w')
    Hap_chroms = []
    with open(genomeSize,'r') as f:
        for line in f:
            line = line.strip().split()
            line[0] = line[0].lstrip('chr')
            c = line[0]
            if ((not chroms) or (c.isdigit() and ('#' in chroms)) or (c in chroms)):
                Hap_chroms.append('M'+c)
                Hap_chroms.append('P'+c)
                out.writelines('M'+'\t'.join(line)+'\n')
                out.writelines('P'+'\t'.join(line)+'\n')
    out.close()

    for res in wholeRes:
        Matrix = All_Data['UnImputated_Whole'][res]['Matrix']
        Bins = All_Data['UnImputated_Whole'][res]['Bins']
        Whole_Lib[res]  = WholeMatrixToSparseDict(Bins, Matrix)
    
    for res in localRes:
        Local_Lib[res] = IntraMatrixToSparseDict(All_Data['UnImputated_Local'][res])
    
    UnImputated_cooler = os.path.join(CoolerPath,prefix+'UnImputated_Haplotype_Multi.cool')

    NPZ2Cooler(datasets = Whole_Lib,
               outfil = UnImputated_cooler,
               genomeSizes_file = Hap_genomeSize,
               chroms = Hap_chroms,
               onlyIntra = False,
               dtype = 'int')
               
    NPZ2Cooler(datasets = Local_Lib,
               outfil = UnImputated_cooler,
               genomeSizes_file = Hap_genomeSize,
               chroms = Hap_chroms,
               onlyIntra = True,
               dtype = 'int')
    
    print "    Imputated Matrix start ..."
    log.log(21, "    Imputated Matrix start ...")
    print "    Two-Step balancing for Imputated Matrix ..."
    log.log(21, "    Two-Step balancing for Imputated Matrix ...")
    
    Balanced_Whole_Data = {}
    for res in wholeRes:
        Bins = All_Data['Tradition_Whole'][res]['Bins']
        Hap_Bins = All_Data['Imputated_Whole'][res]['Bins']
        T_M = All_Data['Tradition_Whole'][res]['Matrix']
        H_M = All_Data['Imputated_Whole'][res]['Matrix']
        Balanced_M = GenomeWideMatrixCorrection(Bins_Pos = Bins,
                                                Hap_Bins_Pos = Hap_Bins,
                                                T_M = T_M,
                                                H_M = H_M)        
        Balanced_Whole_Data[res] = WholeMatrixToSparseDict(Hap_Bins, Balanced_M)
    
    Balanced_Local_Data = {}
    Gap_Local = {}
    for res in localRes:
        Tra_Lib = All_Data['Tradition_Local'][res]
        Hap_Lib = All_Data['Imputated_Local'][res]
        Nor_Lib, Gap_Lib = IntraChromMatrixCorrection(Tra_Lib, Hap_Lib)
        Balanced_Local_Data[res] = IntraMatrixToSparseDict(Nor_Lib)
        Gap_Local[str(res)] = Gap_Lib
    
    Gap_fil = os.path.join(CoolerPath,prefix+'Imputated_Gap.npz')
    np.savez(Gap_fil, **Gap_Local)
    
    Imputated_cooler = os.path.join(CoolerPath,prefix+'Imputated_Haplotype_Multi.cool')
    
    NPZ2Cooler(datasets = Balanced_Whole_Data,
               outfil = Imputated_cooler,
               genomeSizes_file = Hap_genomeSize,
               chroms = Hap_chroms,
               onlyIntra = False,
               dtype = 'float')
    
    NPZ2Cooler(datasets = Balanced_Local_Data,
               outfil = Imputated_cooler,
               genomeSizes_file = Hap_genomeSize,
               chroms = Hap_chroms,
               onlyIntra = True,
               dtype = 'float')
    
    print "All Done !!!"
    log.log(21, "All Done !!!")
    
    
