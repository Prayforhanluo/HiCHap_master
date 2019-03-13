# -*- coding: utf-8 -*-
"""
Created on Thu Mar 08 16:42:14 2018

@author: han-luo
"""
from __future__ import division
import  os, logging, re, subprocess,pysam, sys
import multiprocessing
import time
import Bio
import Bio.Restriction
Bio.Restriction  # To shut up the editor warning

log = logging.getLogger(__name__)


def Enzyme_Handle(enzyme):
    """
    """
    legal_str = ['A','-','G','C','T']
    if hasattr(Bio.Restriction,enzyme):
        enzyme_site = eval('Bio.Restriction.%s.site' % enzyme)
        cutsite = eval('Bio.Restriction.%s.charac' % enzyme)[:2]
    else:
        log.log(21,'Enzyme %s could not find in Bio Database ...',enzyme)
        log.log(21,'Parsing the self-defined enzyme %s',enzyme)
        for i in enzyme:
            if i not in legal_str:
                log.error('Illegal string in enzyme ...')
                log.error('Exit')
                sys.exit(1)
        if '-' not in enzyme:
            log.error('No cutsite in enzyme, pls use - to index it.')
            log.error('Exit')
            sys.exit(1)
        else:
            enzyme_site = ''.join(enzyme.split('-'))
            cutsite = (enzyme.index('-'),-enzyme.index('-'))
    log.log(21,'Enzyme sequence is %s',enzyme_site)
    log.log(21,'Enzyme cut site is %s',cutsite)
    return enzyme_site, cutsite


def GetJuncSeqInfo(enzyme_site,cutsite):
    """
    """

    Dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    reverse = ''.join([Dict[i] for i in enzyme_site])
    
    if cutsite[-1]:
        jplus = enzyme_site[:cutsite[-1]] + enzyme_site[cutsite[0]:]
        jminus = reverse[:cutsite[-1]] + reverse[cutsite[0]:]
        jminus = jminus[::-1]
    else:
        jplus = enzyme_site[:None] + enzyme_site[cutsite[0]:]
        jminus = reverse[:None] + reverse[cutsite[0]:]
        jminus = jminus[::-1]    
    
    if jplus == jminus:
        return jplus,jminus,True
    else:
        return jplus,jminus,False


class Read(object):
    """
    """
    def __init__(self,mate,JuncSeqInfo):
        '''
        '''
        

        self.mate = mate         #pysam for a read
        self.JuncSeq = JuncSeqInfo     #JuncSeqInfo
        self.MIN_LEN = 10        #MIN_LENGTH for read
        self._updateReads()
    
    def _getName(self):
        """Get mate Name """
        return '@'+self.mate.query_name
        
    def _getSeq(self):
        """Get mate Seq"""
        return self.mate.seq
    
    def _getFlag(self):
        """Get mate flag """
        return '+'
    
    def _getQuality(self):
        """ Get mate quality"""
        return self.mate.qual
    
    def _updateReads(self):
        """
        """
        name = self._getName()          #first line for a read in fastq format(Name)
        seq = self._getSeq()              #second line(Seq)
        flag = self._getFlag()           #third line(flag)
        qua = self._getQuality()         #fourth line(quality)
        
        MIN_LEN = self.MIN_LEN
        JunLen = len(self.JuncSeq[0])
        
        if self.JuncSeq[-1]:
            site = [int(m.start())  for m in re.finditer(self.JuncSeq[0], seq)]
        else:
            tmp_site = [int(m.start())  for m in re.finditer(self.JuncSeq[0], seq)]
            if len(tmp_site) == 0:
                tmp_site = [int(m.start())  for m in re.finditer(self.JuncSeq[1], seq)]
            site  = tmp_site
            
        if len(site) == 0:
            self.type = 0                       #No cutting site,Never can be rescued
            self.new_mate = ''
        elif len(site) == 1:
            self.type = 1
            
            seq_part_1 = seq[0:site[0]]
            seq_part_2 = seq[site[0]+JunLen:]
            if len(seq_part_1) < MIN_LEN:
                new_mate = [name,seq_part_2,flag,qua[site[0]+JunLen:]]
                self.new_mate = '\n'.join(new_mate)+'\n'
                
            elif len(seq_part_2) < MIN_LEN:
                new_mate = [name,seq_part_1,flag,qua[0:site[0]]]
                self.new_mate = '\n'.join(new_mate)+'\n'
                
            else:
                new_mate = [name+'1',seq_part_1,flag,qua[0:site[0]],
                            name+'2',seq_part_2,flag,qua[site[0]+JunLen:]]
                
                self.new_mate = '\n'.join(new_mate)+'\n'
        else:
            self.type = 2
            self.new_mate = ''
            



def sub_cutting_Reads(bam,out_fil,JuncSeqInfo):
    """
        Extract the unmapped reads and rescue them for next mapping.
    """
    File = pysam.AlignmentFile(bam)    
    outStream = open(out_fil,'w')    
    for read in File:
        if read.is_unmapped:
            new_read = Read(read,JuncSeqInfo)
            outStream.writelines(new_read.new_mate)
        
    outStream.close()

def Cutting_Reads_To_ReMapping(bam_path,out_folder,enzyme,allel_Mark,threads):
    """
        Cut Fastq by haphic rescue mode to improve Data utilization.


 Mode:
        
    1.    R1, R2 both have no ligation sites.

            5'        3'  3'        5'
        R1  ---------->   <----------   R2    --->   No sites split pair
                                       
            
    2.    R1 split Mode   (cut 3'LS site as a candidate):
        
        R1  ----|----->   <----------   R2       | : Ligation site
              A    B           C
        
        A,C interaction pairs.   B as a candidate read for C.
        
        A <= MIN_LENGTH     ----->   R1 split Useless pair
        
        B <= MIN_LENGTH     ----->   R1 split Normal pair (Only keep A,C)
        
        A,B > MIN_LENGTH    ----->   R1 split Reusable pair(Keep A,C as a pair, B as a candidate)
        
    3.    R2 split Mode
        
        R1  ---------->   <----|------  R2
                A            B    C
        
        A,C interaction pairs.   B as a candidate read for A.
        
        C <= MIN_LENGTH     ----->   R2 split Useless pair
        
        B <= MIN_LENGTH     ----->   R2 split Normal pair (Only keep A,C)
        
        C,B > MIN_LENGTH    ----->   R2 split Reusable pair(Keep A,C as a pair, B as a candidate)
        
    4.    Both split Mode
        
        R1  ----|----->   <-----|-----  R2
              A    *         *    B
        
        A,B interaction pairs. No candidate read
        
        A <= MIN_LENGTH or B <= MIN_LENGTH    -----> Both split Useless pair
        
        A,B > MIN_LENGTH                      -----> Both split Normal pair
        
    5.    Confuse Mode
        
        R1  --|---|---> or <---|---|--  R2
        
            Confuse pairs
    """
    if allel_Mark == 'NonAllelic':
        chunks = [i for i in os.listdir(bam_path) if 'chunk' in i]
    else:
        chunks = [i for i in os.listdir(bam_path) if allel_Mark in i]
    enzyme_site, cutsite = Enzyme_Handle(enzyme)
    JuncSeqInfo = GetJuncSeqInfo(enzyme_site,cutsite)
    if JuncSeqInfo[-1]:
        log.log(21,'junction sequence is %s',JuncSeqInfo[0])
    else:
        log.log(21,'junction sequence plus is %s',JuncSeqInfo[0])
        log.log(21,'junction sequence minus is %s',JuncSeqInfo[1])
    pool = multiprocessing.Pool(threads)
    
    for fil in chunks:
        bam_fil = os.path.join(bam_path,fil)
        out_name = fil.replace('.bam','_unmapped.fq')
        out_fil = os.path.join(out_folder,out_name)
        
        pool.apply_async(sub_cutting_Reads,args=(bam_fil,out_fil,JuncSeqInfo))
    
    time.sleep(2)
    pool.close()
    pool.join()



    
    
def gzipWriter(filename):
    """
    Create a writing process with gzip or parallel gzip (pigz) attached to
    it.
    
    """
    filename = os.path.abspath(filename)
    
    with open(filename, 'wb') as outFile:
        if commandExists("pigz"):
            writer = ["pigz", "-c", "-4"]
        else:
            writer = ["gzip", "-c", "-1"]

        pwrite = subprocess.Popen(writer, stdin = subprocess.PIPE,
                                  stdout = outFile,shell = False,
                                  bufsize = -1)
    return pwrite
    
def commandExists(command):
    "Check if the bash command exists"
    command = command.split()[0]
    if subprocess.call(['which', command],stderr=subprocess.PIPE) != 0:
        return False
    return True

def add_mate_mark(line,mate):
    """
        Add mate mark on read name.
    """
    line = line.split()
    line[0] = line[0]+'_'+str(mate)
    line = ' '.join(line)+'\n'
    
    return line

def Normal_Reads_Split(fq, folder, splitBy,mate):
    """
        Split the reads into chunks.
    
    Parameters
    ----------
    fq : str
        mate of pair-end data
        
    folder : str
        Output Folder
    
    splitBy : int
        number of reads for a chunk
        default : 4 million
    
    mate : int
        number of read mate flag(1,2)
    """
    #log.log(21,'Split reads into chunks ...')
    log.log(21,'Split reads %s into chunks ...',fq)
    
    inFile = os.path.abspath(fq)
    
    parse = os.path.split(inFile)[1].split('.')[0].split('_')
    outFile = '_'.join(parse[:-1]) + '_chunk{0}_{1}.fastq.gz' 
    
    if inFile.endswith('.fastq.gz'):
        pread = subprocess.Popen(['gunzip', inFile, '-c'],
                                  stdout = subprocess.PIPE, bufsize = -1)
    else:
        pread = subprocess.Popen(['cat', inFile],
                                  stdout = subprocess.PIPE, bufsize = -1)

                             
    inStream = pread.stdout
    
    halted = False
    counters = []
    for counter in xrange(1000000):
        
        outProc = gzipWriter(os.path.join(folder,outFile).format(counter, parse[-1]))
        outStream = outProc.stdin
        for j in xrange(splitBy):
            
            line = inStream.readline()

            try:
                assert line[0] == '@'
                assert line[0] == '@'
            except AssertionError:
                raise IOError('{0} is not a fastq file'.format(fq))
            except IndexError:
                halted = True
                counters.append(j)
                break
            
            line = add_mate_mark(line,mate)
            
            fastq_entry = (line, inStream.readline(), inStream.readline(),
                           inStream.readline())
            
            outStream.writelines(fastq_entry)
        
        outProc.communicate()
        counters.append(splitBy)
        if halted:
            log.log(21,'Split %s Down',fq)
            return counters
        
        counters.append(splitBy)
    log.log(21,'Split %s Down',fq)
    return counters

