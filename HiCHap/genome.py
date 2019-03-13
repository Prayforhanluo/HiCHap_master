# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 09:26:36 2018

@author: han-luo
"""
import numpy as np
import os, logging, re, subprocess, bisect, cPickle
import Bio
import Bio.Restriction
Bio.Restriction         #shut up the editor warning

from HapHiC.fastqPlus import Enzyme_Handle

log = logging.getLogger(__name__)



def buildIndex(genome, OutPath, threads):
    """
        Build Genome Index.
    """
    genome_fil = os.path.split(genome)[-1]
    bowtie_Index = os.path.join(OutPath,genome_fil.rstrip('.fa'))
    log.log(21,'Building Index %s...', bowtie_Index)
    
    buildcmd = ['bowtie2-build','--threads',str(threads), genome, bowtie_Index]
    log.log(21,'command : %s',' '.join(buildcmd))
    subprocess.call(buildcmd)
    log.log(21,'BowtieIndex %s, Done', bowtie_Index)
    

    
def enzymeFind(genome, enzymeName, OutPath):
    """
        Finds restriction sites.
        Create a enzyme fragment file.Format and rule like :
        eg : 
            seq : 1(chr)  CTAGATCAATTCGATCTTAC    enzyme : GATC
            the file will be :
            1 1 4    ->   [1,4)
            1 4 13   ->   [4,13) 
            1 13 20  ->   [13,20)
        
    Parameters
    ----------
    genome : str
        Original genome file
    
    enzymeName : str
        enzyme Site
        
    """
    
    # loading Genome        
    dict_genome = {}
    file_genome = open(genome,'r')
    genome_name = genome.split('/')[-1].rstrip('.fa')
    for line in file_genome:
        line = line.strip('\n')
        lists = list(line)
        if len(lists) > 1 and lists[0] == '>':
            chrs = (line.split('>')[1]).split()[0].lstrip('chr')
            dict_genome[chrs] = []
        else:
            dict_genome[chrs].extend(lists)
    # Finding enzyme site

    enzyme_site, cutsite = Enzyme_Handle(enzymeName)
    
    enzyme_file = os.path.join(OutPath,enzymeName+'_'+genome_name+'_fragments.txt')
    f = open(enzyme_file,'w')
    for k in sorted(dict_genome.keys()):
        seq_str = ''.join(dict_genome[k]).upper()
        pos_list = [int(m.start() + 1 + cutsite[0]) for m in re.finditer(enzyme_site, seq_str)]
        pos_list = [1] + pos_list + [len(seq_str)]
        for i in range(len(pos_list)-1):
            line = [k, str(pos_list[i]), str(pos_list[i+1])]
            f.writelines('\t'.join(line)+'\n')
    f.close()






def SNPs_integration(SNP_file,outPath):
    """
        Integrate The SNP information.
    """
    
    log.log(21,'Loading SNPs And  creat temp file under %s',outPath)
    Whole_SNPs = {}
    SNP_f = open(SNP_file,'r')
    for line in SNP_f:
        line = line.split()
        if line[0] not in Whole_SNPs.keys():
            Whole_SNPs[line[0]] = {}
            Whole_SNPs[line[0]]['pos'] = [int(line[1])]
            Whole_SNPs[line[0]]['ref'] = [line[2]]
            Whole_SNPs[line[0]]['m_alt'] = [line[3]]
            Whole_SNPs[line[0]]['p_alt'] = [line[4]]
        else:
            mid_pos = bisect.bisect_left(Whole_SNPs[line[0]]['pos'],int(line[1]))
            Whole_SNPs[line[0]]['pos'].insert(mid_pos,int(line[1]))
            Whole_SNPs[line[0]]['ref'].insert(mid_pos,line[2])
            Whole_SNPs[line[0]]['m_alt'].insert(mid_pos,line[3])
            Whole_SNPs[line[0]]['p_alt'].insert(mid_pos,line[4])
    
    for c in Whole_SNPs.keys():
        Whole_SNPs[c]['pos'] = np.array(Whole_SNPs[c]['pos'])
        Whole_SNPs[c]['ref'] = np.array(Whole_SNPs[c]['ref'])
        Whole_SNPs[c]['p_alt'] = np.array(Whole_SNPs[c]['p_alt'])
        Whole_SNPs[c]['m_alt'] = np.array(Whole_SNPs[c]['m_alt'])

    f = open(os.path.join(outPath,'Snps.pickle'),'wb')
    cPickle.dump(Whole_SNPs,f,2)
    f.close()
    #np.savez_compressed(os.path.join(outPath,'Snps.pickle'),**Whole_SNPs)
    log.log(21,'Done!,temp Snps file %s',os.path.join(outPath,'Snps.pickle'))




def OutputGenome(dict_genome,outfile):
    """
    """
    file_out = open(outfile,'w')
    for k in sorted(dict_genome.keys()):
        out_list = dict_genome[k]
        lens = len(out_list)
        n = int(lens) / 60 + 1
        file_out.write('>'+'chr'+k+' dna:chromosome chromosome:HapHiC:1:1:'\
                        +str(lens)+':1 REF'+'\n')
        for i in range(n):
            strs = ''.join(out_list[(i*60) : ((i+1)*60)])
            file_out.write(strs+'\n')


def buildRawGenome(Input, enzyme, OutPath, threads):
    """
        Build the raw genome Index, A simple way to non-Allelic HiC process.    
    
    """
    log.log(21,'Generate the genome size file')
    dict_genome = {}
    file_genome = open(Input,'r')
    for line in file_genome:
        line = line.strip('\n')
        lists = list(line)
        if len(lists) > 1 and lists[0] == '>':
            chrs = (line.split('>')[1]).split()[0].lstrip('chr')
            dict_genome[chrs] = []
        else:
            dict_genome[chrs].extend(lists)
    
    #Generate the genomeSzie file
    with open(os.path.join(OutPath,'genomeSize'),'w') as o:
        for key in sorted(dict_genome.keys()):
            o.writelines('\t'.join([key,str(len(dict_genome[key]))])+'\n')
    
    del dict_genome
    log.log(21,'split the genome %s by enzyme.',Input)
    enzymeFind(Input,enzyme,OutPath)
    
    buildIndex(Input, OutPath, threads)
    log.log(21,'genome %s done!',Input)



def rebuildGenome(genomePath, snpPath, enzyme, OutPath, threads):
    """
    Rebuild the genome by replacing the snp sites. And build BowtieIndex.
    
    Parameters
    ----------
    genomePath : str
        Original genome file.
        
    snpPath : str
        SNP file.
        
    enzyme : str
        enzyme site.
        
    OutPath : str
        Rebuild Genome Folder
    
    """
    #Load the snp
    log.log(21,'Loading SNP from temp Snp file...')
    
    Snps = cPickle.load(open(snpPath,'rb'))
        
    #Load the Genome
    log.log(21,'Loading Genome ...')
    dict_genome = {}
    file_genome = open(genomePath,'r')
    for line in file_genome:
        line = line.strip('\n')
        lists = list(line)
        if len(lists) > 1 and lists[0] == '>':
            chrs = (line.split('>')[1]).split()[0].lstrip('chr')
            dict_genome[chrs] = []
        else:
            dict_genome[chrs].extend(lists)
    
    #Generate the genomeSzie file
    log.log(21,'Generate the genome size file')
    with open(os.path.join(OutPath,'genomeSize'),'w') as o:
        for key in sorted(dict_genome.keys()):
            o.writelines('\t'.join([key,str(len(dict_genome[key]))])+'\n')
    
    #Replacing Maternal Genome
    log.log(21,'Repalcing Maternal Genome ...')
    
    for c in Snps.keys():
        pos = Snps[c]['pos']
        for index, i in enumerate(pos):
            dict_genome[c][i - 1] = Snps[c]['m_alt'][index]
        
    #Output Maternal Genome
    M_out_file = os.path.join(OutPath,'Maternal/Maternal.fa')
    log.log(21,'Output Maternal Genome %s ...',M_out_file)
    OutputGenome(dict_genome,M_out_file)
    
    #Replacing Paternal Genome
    log.log(21,'Replacing Paternal Genome ...')
    
    for c in Snps.keys():
        pos = Snps[c]['pos']
        for index, i in enumerate(pos):
            dict_genome[c][i - 1] = Snps[c]['p_alt'][index]
    
    #Output Paternal Genome
    P_out_file = os.path.join(OutPath,'Paternal/Paternal.fa')
    log.log(21,'Output Paternal Genome %s ...',P_out_file)
    OutputGenome(dict_genome,P_out_file)
    
    del dict_genome
    del Snps
    
    #Building Maternal Index and Fragment file
    M_Out_Path = os.path.split(M_out_file)[0]
    enzymeFind(M_out_file,enzyme,M_Out_Path)
    buildIndex(M_out_file,M_Out_Path, threads)
    
    #Building Paternal Index and Fragment file
    P_Out_Path = os.path.split(P_out_file)[0]
    enzymeFind(P_out_file,enzyme,P_Out_Path)
    buildIndex(P_out_file,P_Out_Path, threads)
    