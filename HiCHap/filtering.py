# -*- coding: utf-8 -*-
"""
Created on Thu Mar 08 15:40:08 2018

@author: han-luo
"""
from __future__ import  division
from collections import  namedtuple
import os, logging, subprocess, time, gc
import heapq, multiprocessing
import numpy as np

log = logging.getLogger(__name__)


"""
    BamProcess makes bam to bed file. formats as :
    
    column    :    Means
         0    :    Pair Name
         1    :    R1 mate Reference
         2    :    R1 mate Strand
         3    :    R1 mate Position
         4    :    R1 mate Length
         5    :    R1 mate AS score
         6    :    R1 mate Fragment Middle point
         7    :    R1 mate SNP Matching num(Non-Allelic results in 0)
         
         8    :    R2 mate Reference
         9    :    R2 mate Strand
        10    :    R2 mate Position
        11    :    R2 mate Length
        12    :    R2 mate AS score
        13    :    R2 mate Fragment Middle point
        14    :    R2 mate SNP Matching num(Non-Allelic results in 0)
        
        (candidate mate if it is possible)
        15    :    Candidate mate Reference
        16    :    Candidate mate Strand
        17    :    Candidate mate Position
        18    :    Candidate mate Length
        19    :    Candidate mate AS score
        20    :    Candidate mate Fragment Middle point
        21    :    Candidate mate SNP Matching num(Non-Allelic results in 0)
        22    :    Candidate Index for which mate.

"""



def uniqueIndex(data):
    """Returns a binary index of unique elements in an array data.
    This method is very memory efficient, much more than np.unique!
    It grabs only 9 extra bytes per record :)
    
    This Function comes from hiclib.(Author: Anton Goloborodko)
    Nice Function. :)
    """

    args = np.argsort(data)
    index = np.zeros(len(data), bool)
    myr = range(0, len(data), len(data) // 50 + 1) + [len(data)]
    for i in xrange(len(myr) - 1):
        start = myr[i]
        end = myr[i + 1]
        dataslice = data.take(args[start:end], axis=0)
        ind = dataslice[:-1] != dataslice[1:]
        index[args[start:end - 1]] = ind
        if end != len(data):
            if data[args[end - 1]] != data[args[end]]:
                index[args[end - 1]] = True
        else:
            index[args[end - 1]] = True

    return index

def Chunk_sort(chunk,OutPath):
    """
    Sorting the chunk bed for duplicate removing.
        sorted by :
            chr1  : chromosome of mate 1     
            flag1 : strand     of mate 1     
            pos1  : position   of mate 1     
                
            chr1  : chromosome of mate 2     
            flag2 : strand     of mate 2     
            pos2  : position   of mate 2
            
        sort order is from top to bottom.
    """
    tmp = []
    f = open(chunk,'r')
    for line in f:
        line = line.split()
        tmp.append(line)
    f.close()
    sorted_tmp = sorted(tmp,key = lambda x : (x[1],x[2],int(x[3]),x[8],x[9],int(x[10])))
        
    chunk_name = os.path.split(chunk)[-1]
    chunk_name = chunk_name.replace('.bed','_Sorted.bed')
    new_chunk = os.path.join(OutPath,chunk_name)
    with open(new_chunk,'w') as o:
        for line in sorted_tmp:
            o.writelines('\t'.join(line)+'\n')
        
        #subprocess.call(['rm',chunk])
        
    return new_chunk

def Merge_Sorted_Files(key, *iterables):
    """
        Merge the sorted chunks.
    """
    Keyed = namedtuple("Keyed",["key","obj"])
    if key is None:
        keyed_iterables = iterables
    else:
        keyed_iterables = [(Keyed(key(obj.split()), obj) for obj in iterable)
                            for iterable in iterables]
    for element in heapq.merge(*keyed_iterables):
        yield element.obj
        
        
        
# A Customized HiC Filtering Class
class cFiltering(object):
    
    def __init__(self, bedPath, Out_Path, threads, Collection = True,Allelic = 'NonAllelic'):
        
        
        #---------------Initialization the parameters--------
        self.bedPath = bedPath
        self.Allelic = Allelic
        self.Out_Path = Out_Path            #Non-Allelic Out Path.
        self.step = 20000000                #dups remove pairs ID number in memory. 
        self.Collection = Collection        #Whether rm the chunk bed files.
        self.Threads = threads              #Threads Number        
        
    def _StrAsciID(self,string):
        """
            return the Ascii ID sum of a string.
            
        """
        return sum(ord(i) for i in string)
    
    def _GetpairID(self,line):
        """
            return a pair of interaticon ID for duplicates filtering.
            
        """
        line = line.strip().split('\t')
        c1 = self._StrAsciID(line[1].lstrip('chr'))
        f1 = line[2];p1 = line[3]
        
        c2 = self._StrAsciID(line[8].lstrip('chr'))
        f2 = line[9];p2 = line[10]
        
        return int(str(c1)+f1+p1+str(c2)+f2+p2)    
    
    def _gcsleep(self):
        """
            sleep for a second, run garbage collector, sleep again.
            I dont konw if it makes sense, but it definitely does not hurts.
        """
        for _ in range(3):
            time.sleep(0.1)
        gc.collect()
        for _ in range(3):
            time.sleep(0.1)
     
    
    
    def _Getprefix(self):
        """
            return the prefix of file name.
        """
        for bed in os.listdir(self.bedPath):
            if 'chunk' in bed:
                prefix = bed.split('chunk')[0]
                return prefix
        
        return 'tmp_'
    
    def _Mergechunks(self):
        """
            Get chunks list and return the mergeing Linux command. 
        """
        if self.Allelic != 'NonAllelic':
            bed_files = [os.path.join(self.bedPath,i) for i in os.listdir(self.bedPath) \
                        if self.Allelic in i and 'chunk' in i]
        else:
            bed_files = [os.path.join(self.bedPath,i) for i in os.listdir(self.bedPath) \
                        if 'chunk' in i]
        return ['cat'] + sorted(bed_files)
    
    
    def _Collection(self):
        """
            Remove the chunks bed. Only HiC_filtered bed retain.
            Return the Remove Linux command
        """
        if self.Allelic != 'NonAllelic':
            bed_files = [os.path.join(self.bedPath,i) for i in os.listdir(self.bedPath)\
                        if 'chunk' in i and self.Allelic in i]
        else:
            bed_files = [os.path.join(self.bedPath,i) for i in os.listdir(self.bedPath) \
                        if 'chunk' in i]
        return ['rm'] + bed_files
    
    def _MergeSortedFiles(self,key=None,*iterables):
        """
        """
        Keyed = namedtuple('Keyed', ['key', 'obj'])
        if key is None:
            keyed_iterables = iterables
        else:
            keyed_iterables = [(Keyed(key(obj.split()), obj) for obj in iterable)
                                for iterable in iterables]
        for element in heapq.merge(*keyed_iterables):
            yield element.obj
        
    
    def _HiC_Sorting(self,key = None):
        """
            Sorting file for HiC duplicates filtering.
        """
        threads = self.Threads
        Allelic = self.Allelic
        if Allelic == 'NonAllelic':
            Allelic = 'chunk'
        bedPath = self.bedPath
        OutPath = self.Out_Path
        
        sort_Pool = multiprocessing.Pool(threads)
        
        chunks = [os.path.join(bedPath,i) for i in os.listdir(bedPath) if Allelic in i]
        log.log(21,'Chunk Sorting ...')
        #print 'Chunk Sorting ...'
        for fil in chunks:
            sort_Pool.apply_async(Chunk_sort,args=(fil,OutPath))
        
        sort_Pool.close()
        sort_Pool.join()
        log.log(21,'Merging Sorted Chunks...')
        #print 'Merging Sorted Chunks ...'
        sorted_chunks = [os.path.join(OutPath,i) for i in os.listdir(OutPath) \
                            if Allelic in i and 'Sorted' in i]
        chunks = []
        for fil in sorted_chunks:
            out = open(fil,'r')
            chunks.append(out)
        
        prefix = self._Getprefix()
        
        if Allelic == 'chunk':
            Total_Out = os.path.join(OutPath,prefix+'NonAllelic.bed')
        else:            
            Total_Out = os.path.join(OutPath,prefix+Allelic+'.bed')
        with open(Total_Out,'w') as Out:
            Out.writelines(Merge_Sorted_Files(lambda x : (x[1],x[2],int(x[3]),x[8],x[9],int(x[10])),*chunks))
        self.Total_Out = Total_Out
        for IO in chunks:
            IO.close()
        collect = ['rm'] + sorted_chunks
        subprocess.call(collect)
        log.log(21,'Sorting Done')
        #print 'Sorting Done'
        
        
        
        
    
    def _Redundant_kind(self,line):
        """
            Pair Tagging if pair is Redundant pair.
            
            Self-Circle:  SC
                Pair deriving from self-circularized ligation product.
                The two reads are mapped to the same restriction fragment and face
                in opposite directions.
                
               R1                            R1
            <----                            ---->              
            ==============    or    ==============
                     ---->          <----
                     R2                R2
            
            DanglingEnds:  DE
                There can be may causes of such products,ranging from low ligation efficiency
                to poor strptavidin specificity.
                Both tow reads are mapped to the same restriction fragment and face
                towards each other.
                
            R1                                  R1
            ---->                            <----  
            ==============    or    ==============
                     <----          ---->
                        R2          R2
            
            UnkownMechanism:  UM
                Unkown sources.
                Both two reads are mapped to the same restriction fragment and same strand.
            
            ---->                   <----
            ==============    or    ==============
                     ---->                   <----
            
            ExtraDanglingEnds:  ED
                The two reads of these pairs are mapped to the different restriction
                fragments but face towards each other and are separated by less than
                the library size (500bp) interval.Such read pairs may contain true 
                contacts, but are largely contaminated, so we also remove these paris.
                
                R1                              R1
                ---->                        <---- 
            =========|======== or   ========|=========
                      <----            ---->
                          R2           R2
        """
        
        line = line.split()
        C1 = line[1];C2 = line[8]                           #chromosome
        strand1 = int(line[2]);strand2 = int(line[9])       #strand
        pos1 = int(line[3]);pos2 = int(line[10])            #position
        Frag1 = int(line[6]);Frag2 = int(line[13])          #Fragments Middle
        
        if C1 != C2:
            return False
        elif Frag1 == Frag2:
            
            if pos1 < pos2:
                if strand1 == 0 and strand2 == 16:
                    return 'DE'
                elif strand1 == 16 and strand2 == 0:
                    return 'SC'
                else:
                    return 'UM'
            else:
                if strand1 == 0 and strand2 == 16:
                    return 'SC'
                elif strand1 == 16 and strand2 == 0:
                    return 'DE'
                else:
                    return 'UM'
        else:
            if abs(pos1 - pos2) <= 500:
                if pos1 < pos2 and strand1 == 0 and strand2 == 16:
                    return 'ED'
                elif pos1 > pos2 and strand1 == 16 and strand2 == 0:
                    return 'ED'
                else:
                    return False
            else:
                return False
    
    
    
    def HiC_Filtering(self):
        """
            HiC Filtering for duplicates, same Fragments reads.
        """
        log.log(21,'HiC Filtering Start for %s Interactions' % self.Allelic)
        #print 'HiC Filtering Start for %s Interations' % self.Allelic
        log.log(21,'Sorting and HiC filtering...')
        
        self._HiC_Sorting()
        log.log(21,'Filtering start ...')
        #print 'Filtering start ...'
        self.Duplicates = 0
        self.SelfCircle = 0
        self.DanglingEnds = 0
        self.UnknownMechiansm = 0
        self.ExtraDanglingEnds = 0
        self.ValidPairs = 0
        prefix = self._Getprefix()
        if self.Allelic != 'NonAllelic':
            self.Outbed = os.path.join(self.Out_Path,prefix+self.Allelic+'_Valid.bed')
            Out_bed = open(os.path.join(self.Out_Path,prefix+self.Allelic+'_Valid.bed'),'w')
        else:
            self.Outbed = os.path.join(self.Out_Path,prefix+'Valid.bed')
            Out_bed = open(os.path.join(self.Out_Path,prefix+'Valid.bed'),'w')
        mark = 0                        #Initial value
        count = 0                       #Step num for Unique
        f = open(self.Total_Out,'r')
        for line in f:
            if (count+1) % 40000000 == 0:
                log.log(21,'%d Pairs have been filtered',count+1)
                #print '%d Pairs have been filtered' % (count + 1)
            count += 1
            
            ID = self._GetpairID(line)
            if ID != mark:
                mark = ID
                PM = self._Redundant_kind(line)
                
                if PM == False:
                    self.ValidPairs += 1
                    Out_bed.writelines(line)
                elif PM == 'SC':
                    self.SelfCircle += 1
                elif PM == 'DE':
                    self.DanglingEnds += 1
                elif PM == 'UM':
                    self.UnknownMechiansm += 1
                elif PM == 'ED':
                    self.ExtraDanglingEnds += 1
            else:
                self.Duplicates += 1
        
        log.log(21,'%d pairs have been filtered',count)
        #print '%d paris have been filtered' % count
        
        Static_results = ['Total Unique Pairs : %d' % count,
                          'Duplicate Pairs : %d' % self.Duplicates,
                          'Valid Pairs  is %d' % self.ValidPairs,
                          'Self Circle Pairs is %d' % self.SelfCircle,
                          'DanglingEnds Pairs is %d' % self.DanglingEnds,
                          'UnkownMechanism Pairs is %d' % self.UnknownMechiansm,
                          'ExtraDanglingEnds Pairs is %d' % self.ExtraDanglingEnds]
        
        Static_results = '\n'*2 + '\n'.join(Static_results) + '\n'
        
        log.log(21,'HiC Statics Results : %s',Static_results)
        #print 'HiC Statics Results : %s' % Static_results
        log.log(21,'HiC filtering for %s Interations Done! \n',self.Allelic)
        
        f.close()
        subprocess.call(['rm',self.Total_Out])
        Out_bed.close()
        
        if self.Collection == True:
            subprocess.call(self._Collection())
        


# A Customized Allelic Filtering Class
class aFiltering(object):
    
    def __init__(self,Maternal_bed,Paternal_bed,Out_Path,save_ID = False):
        
        
        #---------------Initialization the parameters--------
        self.Maternal_bed = Maternal_bed                #Maternal Filtered bed
        self.Paternal_bed = Paternal_bed                #Paternal Filtered bed
        self.Out_Path = Out_Path                        #Out
        self.ChunkStep = 10000000                       #chunk Step for Name sorting.
        self.MAX_DIFF_SCORE = 18                         #Max diff score 
        self.save_ID = save_ID                          #Save the pair ID.
    
    
    def _BigFileSorting(self,bed):
        """ 
            Sorting the Interation bed file(.bed) by reads name created by bamProcess.
        
        """
        
        from itertools import islice, count
        from contextlib2 import ExitStack
        from heapq import merge
        
        log.log(21,'Sorting file  %s',bed)
        #print 'Sorting file %s' % bed
        log.log(21,'chunk file in step line: %d and sorting',self.ChunkStep)
        #print 'chunk file in step line %d' % self.ChunkStep
        chunk_names = []
        
        OutPath = self.Out_Path
        if not os.path.exists(OutPath):
            os.mkdir(OutPath)
            
        prefix = os.path.split(bed)[1].replace('.bed','')
        with open(bed) as input_fil:
            for chunk_number in count(1):
                sorted_chunk = sorted(islice(input_fil,self.ChunkStep))
                if not sorted_chunk:
                    break
                
                chunk_name = prefix+'_chunks_{}.chk'.format(chunk_number)
                chunk_names.append(chunk_name)
                with open(os.path.join(OutPath,chunk_name),'w') as chunk_fil:
                    chunk_fil.writelines(sorted_chunk)
                #log.log(21,'chunk %d completed',chunk_number)
                #print 'chunk %d completed' % chunk_number
        
        log.log(21,'Merge the sorted chunks')
        #print 'Merging...'
        Out = os.path.join(OutPath,prefix + '_sorted.bed')
        with ExitStack() as stack, open(Out,'w') as output_fil:
            files = [stack.enter_context(open(os.path.join(OutPath,chunk))) for chunk in chunk_names]
            output_fil.writelines(merge(*files))
        
        log.log(21,'temp files collection')
        chunk_files = [os.path.join(OutPath,chunk) for chunk in chunk_names]
        subprocess.call(['rm']+chunk_files)
        subprocess.call(['rm']+[bed])
        log.log(21,'Done!')
        #print 'Done!'
        
        return Out
    
    def _intToString(self,lst):
        """
            Transform the int type to string type in the list.
        """
        return [str(i) for i in lst]
    
    def _is_Available_Candidate(self,Info):
        """
            Whether the Candidate mate can be used to Allelic selecting.
            1.
            
            -------|-------   ---------------
            R1            C          R2
            
            2.
            ---------------   -----|---------
                   R1         C            R2
                   
            For 1 :
                Candidate(C) must have the same fragment mapping results as R1.
            
            For 2 :
                Candidate(C) must have the same fragment mapping results as R2.
        """
        mask = False
        candidate = Info[-1]
        if candidate == 'R1':
            C1 = Info[1]
            C2 = Info[15]
            Frag1 = int(Info[6])
            Frag2 = int(Info[20])
            
            if C1 == C2 and Frag1 == Frag2:
                mask = True
                return mask
        else:
            C1 = Info[8]
            C2 = Info[15]
            Frag1 = int(Info[13])
            Frag2 = int(Info[20])
            
            if C1 == C2 and Frag1 == Frag2:
                mask = True
                return mask
        
        return mask
        
    
    
    
    
    def _sub_search(self,M_C, M_pos, M_score, M_SNPs,
                         P_C, P_pos, P_score, P_SNPs):
        """
            Using Parameters(mentioned in Allelic_Filtering) to finish Allelic selecting.
            
            Parameters:
            ----------
            
            M_C : Maternal reference.
            
            M_pos : Maternal mapping position.
            
            M_score : Maternal mapping punish score.
            
            M_SNPs : Maternal SNP sites matching counts.
            
            ...
            
            Paternal Parameters are the same.
            
        """
        
        #Mapping To the Same Position
        if M_C == P_C and abs(M_pos - P_pos) <= 5:
            if M_SNPs > P_SNPs:
                mark = 'M'
            elif M_SNPs < P_SNPs:
                mark = 'P'
            else:
                mark = 'N'
        
        #Mapping To the Diff Position
        else:
            if (M_score - P_score) >= self.MAX_DIFF_SCORE and M_SNPs >= 2 * P_SNPs:
                mark = 'M'
            elif (P_score - M_score) >= self.MAX_DIFF_SCORE and P_SNPs >= 2 * P_SNPs:
                mark = 'P'
            else:
                mark = 'N'
        
        return mark
        




    
    def _Both_Mapping_line_Process(self):
        """
            Using Allelic searching rules for a single mapping results.
            Pair mapped to  both two allelic genome succeed.
            
            Parameters
            ----------
                self
                
            Returns
            ----------
            mark1 + mark2 : str
                NN ---> Bi-Allelic
                NM ---> R2 Maternal
                NP ---> R2 Paternal
                MN ---> R1 Maternal
                PN ---> R1 Paternal
                MM ---> R1,R2 Maternal
                PP ---> R1,R2 Paternal
                MP,PM ---> Regroup
            
            lines : list
                Information without ['both','R1','R2'] for IO stream.
            
        """
        
        self.M_1_C = self.M_Info[1];self.M_2_C = self.M_Info[8]
        self.P_1_C = self.P_Info[1];self.P_2_C = self.P_Info[8]
        
        self.M_1_pos = int(self.M_Info[3]);self.M_2_pos = int(self.M_Info[10])
        self.P_1_pos = int(self.P_Info[3]);self.P_2_pos = int(self.P_Info[10])
        
        self.M_1_Frag = int(self.M_Info[6]);self.M_2_Frag = int(self.M_Info[13])
        self.P_1_Frag = int(self.P_Info[6]);self.P_2_Frag = int(self.P_Info[13])
        
        self.M_1_score = int(self.M_Info[5]);self.M_2_score = int(self.M_Info[12])
        self.P_1_score = int(self.P_Info[5]);self.P_2_score = int(self.P_Info[12])
        
        self.M_1_SNPs = int(self.M_Info[7]);self.M_2_SNPs = int(self.M_Info[14])
        self.P_1_SNPs = int(self.P_Info[7]);self.P_2_SNPs = int(self.P_Info[14])
        
        if len(self.M_Info) == len(self.P_Info) == 15:
            
            #For R1
            mark1 = self._sub_search(self.M_1_C,self.M_1_pos,self.M_1_score,self.M_1_SNPs,
                                     self.P_1_C,self.P_1_pos,self.P_1_score,self.P_1_SNPs)
            
            if mark1 == 'N' or mark1 == 'M':
                line1 = [self.M_1_C,self.M_1_Frag]
            else:
                line1 = [self.P_1_C,self.P_1_Frag]
                
            #For R2
            mark2 = self._sub_search(self.M_2_C,self.M_2_pos,self.M_2_score,self.M_2_SNPs,
                                     self.P_2_C,self.P_2_pos,self.P_2_score,self.P_2_SNPs)
            
            if mark2 == 'N' or mark2 == 'M':
                line2 = [self.M_2_C,self.M_2_Frag]
            else:
                line2 = [self.P_2_C,self.P_2_Frag]
            
            line = line1 + line2
            
        
        elif len(self.M_Info) > 15 and len(self.P_Info) == 15:            
            candidate_mark = self._is_Available_Candidate(self.M_Info)
            candidate = self.M_Info[-1]
            
            mark1 = self._sub_search(self.M_1_C,self.M_1_pos,self.M_1_score,self.M_1_SNPs,
                                     self.P_1_C,self.P_1_pos,self.P_1_score,self.P_1_SNPs)
            
            if mark1 == 'N' or mark1 == 'M':
                line1 = [self.M_1_C,self.M_1_Frag]
            else:
                line1 = [self.P_1_C,self.P_1_Frag]                         
            
            mark2 = self._sub_search(self.M_2_C,self.M_2_pos,self.M_2_score,self.M_2_SNPs,
                                     self.P_2_C,self.P_2_pos,self.P_2_score,self.P_2_SNPs)
            
            if mark2 == 'N' or mark2 == 'M':
                line2 = [self.M_2_C,self.M_2_Frag]
            else:
                line2 = [self.P_2_C,self.P_2_Frag]
                
                
            if candidate == 'R1' and candidate_mark == True:
                if mark1 == 'N':
                    self.M_1_C = self.M_Info[15]
                    self.M_1_pos = int(self.M_Info[17])
                    self.M_1_Frag = int(self.M_Info[20])
                    self.M_1_score = int(self.M_Info[19])
                    self.M_1_SNPs = int(self.M_Info[21])
                    mark1 = self._sub_search(self.M_1_C,self.M_1_pos,self.M_1_score,self.M_1_SNPs,
                                             self.P_1_C,self.P_1_pos,self.P_1_score,self.P_1_SNPs)
                    
                    if mark1 == 'M':
                        line1 = [self.M_1_C,self.M_1_Frag]
                    elif mark1 == 'P':
                        line1 = [self.P_1_C,self.P_1_Frag]
                    else:
                        pass
                        
                else:
                    pass
                
            elif candidate == 'R2' and candidate_mark == True:
                if mark2 == 'N':
                    self.M_2_C = self.M_Info[15]
                    self.M_2_pos = int(self.M_Info[17])
                    self.M_2_Frag = int(self.M_Info[20])
                    self.M_2_score = int(self.M_Info[19])
                    self.M_2_SNPs = int(self.M_Info[21])
                    mark2 = self._sub_search(self.M_2_C,self.M_2_pos,self.M_2_score,self.M_2_SNPs,
                                             self.P_2_C,self.P_2_pos,self.P_2_score,self.P_2_SNPs)
                    
                    if  mark2 == 'M':
                        line2 = [self.M_2_C,self.M_2_Frag]
                    elif mark2 == 'P':
                        line2 = [self.P_2_C,self.P_2_Frag]
                    else:
                        pass
                        
                else:
                    pass
            
            else:
                pass
            
            line = line1 + line2
            
        
        elif len(self.M_Info) == 15 and len(self.P_Info) >= 15:
            candidate_mark = self._is_Available_Candidate(self.P_Info)
            candidate = self.P_Info[-1]
            
            mark1 = self._sub_search(self.M_1_C,self.M_1_pos,self.M_1_score,self.M_1_SNPs,
                                     self.P_1_C,self.P_1_pos,self.P_1_score,self.P_1_SNPs)
            
            if mark1 == 'N' or mark1 == 'M':
                line1 = [self.M_1_C,self.M_1_Frag]
            else:
                line1 = [self.P_1_C,self.P_1_Frag]
                         
            mark2 = self._sub_search(self.M_2_C,self.M_2_pos,self.M_2_score,self.M_2_SNPs,
                                     self.P_2_C,self.P_2_pos,self.P_2_score,self.P_2_SNPs)
            
            if mark2 == 'N' or mark2 == 'M':
                line2 = [self.M_2_C,self.M_2_Frag]
            else:
                line2 = [self.P_2_C,self.P_2_Frag]
            
            
            if candidate == 'R1' and candidate_mark == True:
                if mark1 == 'N':
                    self.P_1_C = self.P_Info[15]
                    self.P_1_pos = int(self.P_Info[17])
                    self.P_1_Frag = int(self.P_Info[20])
                    self.P_1_score = int(self.P_Info[19])
                    self.P_1_SNPs = int(self.P_Info[21])
                    mark1 = self._sub_search(self.M_1_C,self.M_1_pos,self.M_1_score,self.M_1_SNPs,
                                             self.P_1_C,self.P_1_pos,self.P_1_score,self.P_1_SNPs)
                    
                    if mark1 == 'M':
                        line1 = [self.M_1_C,self.M_1_Frag]
                    elif mark1 == 'P':
                        line1 = [self.P_1_C,self.P_1_Frag]
                    else:
                        pass
                    
                else:
                    pass
            elif candidate == 'R2' and candidate_mark == True:
                if mark2 == 'N':
                    self.P_2_C = self.P_Info[15]
                    self.P_2_pos = int(self.P_Info[17])
                    self.P_2_Frag = int(self.P_Info[20])
                    self.P_2_score = int(self.P_Info[19])
                    self.P_2_SNPs = int(self.P_Info[21])
                    mark2 = self._sub_search(self.M_2_C,self.M_2_pos,self.M_2_score,self.M_2_SNPs,
                                             self.P_2_C,self.P_2_pos,self.P_2_score,self.P_2_SNPs)
                    
                    if mark2 == 'M':
                        line2 = [self.M_2_C,self.M_2_Frag]
                    elif mark2 == 'P':
                        line2 = [self.P_2_C,self.P_2_Frag]
                    else:
                        pass
                        
                else:
                    pass
            else:
                pass
            
            line = line1 + line2
        
        else:
            candidate_mark_M = self._is_Available_Candidate(self.M_Info)
            candidate_mark_P = self._is_Available_Candidate(self.P_Info)
            
            candidate = self.M_Info[-1]     #Both have candidate. And will for the same mate.
            
            mark1 = self._sub_search(self.M_1_C,self.M_1_pos,self.M_1_score,self.M_1_SNPs,
                                     self.P_1_C,self.P_1_pos,self.P_1_score,self.P_1_SNPs)
            
            
            if mark1 == 'N' or mark1 == 'M':
                line1 = [self.M_1_C,self.M_1_Frag]
            else:
                line1 = [self.P_1_C,self.P_1_Frag]

                         
            mark2 = self._sub_search(self.M_2_C,self.M_2_pos,self.M_2_score,self.M_2_SNPs,
                                     self.P_2_C,self.P_2_pos,self.P_2_score,self.P_2_SNPs)
            
            
            if mark2 == 'N' or mark2 == 'M':
                line2 = [self.M_2_C,self.M_2_Frag]
            else:
                line2 = [self.P_2_C,self.P_2_Frag]
                
            
            if candidate == 'R1' and mark1 == 'N':
                if candidate_mark_M == True:
                    self.M_1_C = self.M_Info[15]
                    self.M_1_pos = int(self.M_Info[17])
                    self.M_1_Frag = int(self.M_Info[20])
                    self.M_1_score = int(self.M_Info[19])
                    self.M_1_SNPs = int(self.M_Info[21])
                    
                if candidate_mark_P == True:
                    self.P_1_C = self.P_Info[15]
                    self.P_1_pos = int(self.P_Info[17])
                    self.P_1_Frag = int(self.P_Info[20])
                    self.P_1_score = int(self.P_Info[19])
                    self.P_1_SNPs = int(self.P_Info[21])
                                
                mark1 = self._sub_search(self.M_1_C,self.M_1_pos,self.M_1_score,self.M_1_SNPs,
                                         self.P_1_C,self.P_1_pos,self.P_1_score,self.P_1_SNPs)
                
                if mark1 == 'M':
                    line1 = [self.M_1_C,self.M_1_Frag]
                elif mark1 == 'P':
                    line1 = [self.P_1_C,self.P_1_Frag]
                else:
                    pass
            
            
            elif candidate == 'R2' and mark2 == 'N':
                if candidate_mark_M == True:
                    self.M_2_C = self.M_Info[15]
                    self.M_2_pos = int(self.M_Info[17])
                    self.M_2_Frag = int(self.M_Info[20])
                    self.M_2_score = int(self.M_Info[19])
                    self.M_2_SNPs = int(self.M_Info[21])
                    
                if candidate_mark_P == True:
                    self.P_2_C = self.P_Info[15]
                    self.P_2_pos = int(self.P_Info[17])
                    self.P_2_Frag = int(self.P_Info[20])
                    self.P_2_score = int(self.P_Info[19])
                    self.P_2_SNPs = int(self.P_Info[21])
                
                mark2 = self._sub_search(self.M_2_C,self.M_2_pos,self.M_2_score,self.M_2_SNPs,
                                         self.P_2_C,self.P_2_pos,self.P_2_score,self.P_2_SNPs)
    
                if mark2 == 'M':
                    line2 = [self.M_2_C,self.M_2_Frag]
                elif mark2 == 'P':
                    line2 = [self.P_2_C,self.P_2_Frag]
                else:
                    pass
            
            else:
                pass
            
            line = line1 + line2
        
        if self.save_ID == True:
            
            return mark1 + mark2, [self.M_Info[0]]+line
        
        else:
            return mark1 + mark2, line
    
    
    
    
    
    
    def _Specific_Mapping_line_Process(self,Info):
        """
            When pair Only mapped to one Allelic Genome.
            
            Complete the Allelic selected Process.
            
            Parameters
            ----------
            Info : list
                Mapping Info list.
                
            Return
            ---------
            mark : str
                both ---> both mates can be allelic selected.
                R1   ---> Only R1 can be allelic selected.
                R2   ---> Only R2 can be allelic selected.
            
            lines : list
                Information that will be written to the IO stream.
                
        """
        SNP_1 = int(Info[7])
        SNP_2 = int(Info[14])
        
        lines = [Info[1],Info[6],Info[8],Info[13]]
        
        if SNP_1 != 0 and SNP_2 != 0:
            mark = 'Both'
            lines += [mark]
        
        elif SNP_1 != 0 and SNP_2 == 0:
            try:
                candidate = Info[-1]
                SNP_3 = int(Info[21])
                candidate_mark = self._is_Available_Candidate(Info)
                
                if candidate == 'R2' and candidate_mark == True:
                    if SNP_3 != 0:
                        mark = 'Both'
                        lines = [Info[1],Info[6],Info[15],Info[20]]
                        lines += [mark]
                    else:
                        mark = 'R1'
                        lines += [mark]
                else:
                    mark = 'R1'
                    lines += [mark]
            except:
                mark = 'R1'
                lines += [mark]
        
        elif SNP_1 == 0 and SNP_2 != 0:
            try:
                candidate = Info[-1]
                SNP_3 = int(Info[21])
                candidate_mark = self._is_Available_Candidate(Info)

                if candidate== 'R1' and candidate_mark == True:
                    if SNP_3 != 0:
                        mark = 'Both'
                        lines = [Info[15],Info[20],Info[8],Info[13]]
                        lines += [mark]
                    else:
                        mark = 'R2'
                        lines += [mark]
                else:
                    mark = 'R2'
                    lines += [mark]
            except:
                mark = 'R2'
                lines += [mark]
        else:
            #When normal sides cant figure out. Use the Candidate side.
            try:
                candidate = Info[-1]
                candidate_mark = self._is_Available_Candidate(Info)
                SNP_3 = int(Info[21])
                if candidate == 'R1' and SNP_3 != 0 and candidate_mark == True:
                    lines = [Info[15],Info[20],Info[8],Info[13]]
                    mark = 'R1'
                    lines += [mark]
                elif candidate == 'R2' and SNP_3 != 0 and candidate_mark == True:
                    lines = [Info[1],Info[6],Info[15],Info[20]]
                    mark = 'R2'
                    lines += [mark]
                else:
                    mark = 'N'
            except:
                mark = 'N'
        
        if self.save_ID == True:
            return mark, [Info[0]] + lines
        
        else:            
            return mark, lines
                
    
    
    
    
    def Allelic_Filtering(self):
        """
        
        Allelic Tagging.
            
            Pair is Specifical mapped to Paternal Genome:
                if Pair has Paternal SNP:
                    R1 has SNP and R2 has SNP --->   Paternal, Both
                    
                    R1 has SNP and R2 no SNP  --->   Paternal, R1
                    
                    R1 no SNP and R2 has SNP --->    Paternal, R2
                else:
                    Bi-Allelic Pair
            
            Pair is Specifical mapped to Maternal Genome:
                if Pair has Maternal SNP:
                    R1 has SNP and R2 has SNP --->   Maternal, Both
                    
                    R1 has SNP and R2 no SNP  --->   Maternal, R1
                    
                    R1 no SNP and R2 has SNP --->    Maternal, R2
                else:
                    Bi-Allelic Pair
            
            Pair mapped to Both Maternal Genome and Paternal Genome.
                
                =======      =======
                R1                R2
            
            Allelic Tagging is same for two mates.
            
            First,We will check whether mate mapped to the same position on Maternal
            genome and Paternal genome.
                
            Parameters:
            ----------
            pos_M    :    Maternal position
            pos_P    :    Paternal position
            MAX_DIFF_SCORE : the max diff score allowed bewteen two mapping results.
            score_M  :   Mapping punish score on Maternal Genome (AS)
            score_P  :   Mapping punish score on Paternal Genome (AS)
            SNP_M    :   SNP matching sites on Maternal Genome for this mate.
            SNP_P    :   SNP matching sites on Paternal Genome for this mate.
            
            Both Mapping Rules:
            -------------------
            
            if pos_M == pos_P,We first check the SNP matching sites counts. 
                SNP_M > SNP_P  ---->  Maternal mate.
                SNP_M < SNP_P  ---->  Paternal mate.
                when SNP_M == SNP_P == 0,We check the difference score.
                    check the difference of score.
                    score_M - score_P > MAX_DIFF_SCORE ----> Maternal mate
                    score_P - score_M > MAX_DIFF_SCORE ----> Paternal mate
            
            if pos_M != pos_P, We first check the difference score. if difference 
            score satisfiy the MAX_DIFF_SCORE.We select the best score position is 
            true position.Then We check the SNP Macthing sites. if nonzero, selected
            the corresponding allelic. 
            
            
            When the rules do not separate the Allelic Reads(Any mates). And the Candidates(C) is
            available. We can Use candidate mate to excavate more potential Allelic
            Information.
            
                ========      ========
                C(M)               C(P)
                
            For Candidates:
                We will replace the traditional mate with the corresponding candidate
                mate. And Repeat the the above process.
                
                
        """
        
        log.log(21,'Sorting for HiC results ...')
        #print 'Sorting for HiC results'
        
        M_bed = self._BigFileSorting(self.Maternal_bed)
        P_bed = self._BigFileSorting(self.Paternal_bed)
        
        log.log(21,'Allelic Filtering ...')
        #print 'Allelic Filtering'
        
        #M_bed = self.Maternal_bed
        #P_bed = self.Paternal_bed
        M_fil = open(M_bed,'r')
        P_fil = open(P_bed,'r')
        #prefix = self.Maternal_bed.split('Maternal')[0]
        #prefix = os.path.split(self.Maternal_bed)[-1].replace('.bed','')
        prefix = os.path.split(self.Maternal_bed)[-1].split('Maternal')[0]+'Valid'
        
        #Out Put
        Bi_Allelic_Out = open(os.path.join(self.Out_Path,prefix+'_Bi_Allelic.bed'),'w')
        M_M_Out = open(os.path.join(self.Out_Path,prefix+'_M_M.bed'),'w')
        P_P_Out = open(os.path.join(self.Out_Path,prefix+'_P_P.bed'),'w')
        M_P_Out = open(os.path.join(self.Out_Path,prefix+'_M_P.bed'),'w')
        P_M_Out = open(os.path.join(self.Out_Path,prefix+'_P_M.bed'),'w') 

        #Statics dic
        
        self.Static_dict = {}
        
        Bi_Allelic = 0
        Both_M = 0
        Both_P = 0
        Single_M = 0
        Single_P = 0
        Regroup = 0   
        Speci_M = 0
        Speci_P = 0
        Speci_M_single = 0
        Speci_M_both = 0
        Speci_P_single = 0
        Speci_P_both = 0
        
        M_line = M_fil.readline()
        P_line = P_fil.readline()
        count = 0 
        
        while True:
            count += 1
            if count % 40000000 == 0:
                log.log(21,'%d pairs completed',count)
                #print '%d pairs completed' % count
            
            self.M_Info = M_line.strip().split()
            self.P_Info = P_line.strip().split()
            
            #All file come to the END. STOP
            if len(self.M_Info) == 0 and len(self.P_Info) == 0:
                break
            
            #Only mapped to Paternal and appear in the END line.
            elif len(self.M_Info) == 0 and len(self.P_Info) != 0:
                Speci_P += 1
                mark, lines = self._Specific_Mapping_line_Process(self.P_Info)
                if mark == 'Both':
                    Both_P += 1
                    Speci_P_both += 1
                    P_P_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                elif mark == 'R1' or mark == 'R2':
                    Single_P += 1
                    Speci_P_single += 1
                    P_P_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                else:
                    Bi_Allelic += 1
                    Bi_Allelic_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                
                P_line = P_fil.readline()

            #Only mapped to Maternal and appear in the END line.
            elif len(self.M_Info) != 0 and len(self.P_Info) == 0:
                Speci_M += 1
                mark, lines = self._Specific_Mapping_line_Process(self.M_Info)
                if mark == 'Both':
                    Both_M += 1
                    Speci_M_both += 1
                    M_M_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                elif mark == 'R1' or mark == 'R2':
                    Single_M += 1
                    Speci_M_single += 1
                    M_M_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                else:
                    Bi_Allelic += 1
                    Bi_Allelic_Out.writelines('\t'.join(self._intToString(lines))+'\n')

                M_line = M_fil.readline()
            
            else:
                M_Name = self.M_Info[0];P_Name = self.P_Info[0]
                
                #Only mapped to Maternal and appears in the middle of file.
                if M_Name < P_Name:
                    Speci_M += 1
                    mark, lines = self._Specific_Mapping_line_Process(self.M_Info)
                    
                    if mark == 'Both':
                        Both_M += 1
                        Speci_M_both += 1
                        M_M_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                    elif mark == 'R1' or mark == 'R2':
                        Single_M += 1
                        Speci_M_single += 1
                        M_M_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                    else:
                        Bi_Allelic += 1
                        Bi_Allelic_Out.writelines('\t'.join(self._intToString(lines))+'\n')

                    M_line = M_fil.readline()
                
                #Only mapped to Paternal and appears in the middle of file.
                elif M_Name > P_Name:
                    Speci_P += 1
                    mark, lines = self._Specific_Mapping_line_Process(self.P_Info)
                    
                    if mark == 'Both':
                        Both_P += 1
                        Speci_P_both += 1
                        P_P_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                    elif mark == 'R1' or mark == 'R2':
                        Single_P += 1
                        Speci_P_single += 1
                        P_P_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                    else:
                        Bi_Allelic += 1
                        Bi_Allelic_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                        
                    P_line = P_fil.readline()
                
                #Pair mapped to the two genome.
                else:
                    mark, lines = self._Both_Mapping_line_Process()
                    if mark == 'NN':
                        Bi_Allelic += 1
                        Bi_Allelic_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                        
                    elif mark == 'NM':
                        Single_M += 1
                        M_M_Out.writelines('\t'.join(self._intToString(lines)+['R2'])+'\n')
                    
                    elif mark == 'MN':
                        Single_M += 1
                        M_M_Out.writelines('\t'.join(self._intToString(lines)+['R1'])+'\n')
                        
                    elif mark == 'MM':
                        Both_M += 1
                        M_M_Out.writelines('\t'.join(self._intToString(lines)+['Both'])+'\n')
                        
                    elif mark == 'NP':
                        Single_P += 1
                        P_P_Out.writelines('\t'.join(self._intToString(lines)+['R2'])+'\n')
                        
                    elif mark == 'PN':
                        Single_P += 1
                        P_P_Out.writelines('\t'.join(self._intToString(lines)+['R1'])+'\n')
                        
                    elif mark == 'PP':
                        Both_P += 1
                        P_P_Out.writelines('\t'.join(self._intToString(lines)+['Both'])+'\n')
                        
                    elif mark == 'MP':
                        Regroup += 1
                        M_P_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                        
                    elif mark == 'PM':
                        Regroup += 1
                        P_M_Out.writelines('\t'.join(self._intToString(lines))+'\n')
                    
                    M_line = M_fil.readline()
                    P_line = P_fil.readline()
        
        log.log(21,"%d pair completed", count-1)
        #print '%d pairs completed' % (count - 1)
        
        #Close the file
        Bi_Allelic_Out.close()
        M_M_Out.close()
        P_P_Out.close()
        M_P_Out.close()
        P_M_Out.close()
        
        #Statics dict
        self.Static_dict['Total_valid_pairs'] = count - 1
        self.Static_dict['Bi_Allelic_pairs'] = Bi_Allelic
        self.Static_dict['Maternal_Allelic_pairs'] = Both_M + Single_M
        self.Static_dict['Paternal_Allelic_pairs'] = Both_P + Single_P
        self.Static_dict['Maternal_both_sides_pairs'] = Both_M
        self.Static_dict['Paternal_both_sides_pairs'] = Both_P
        self.Static_dict['Maternal_single_side_pairs'] = Single_M
        self.Static_dict['Paternal_single_side_pairs'] = Single_P
        self.Static_dict['Speci_Maternal_Mapping_pairs'] = Speci_M
        self.Static_dict['Speci_Paternal_Mapping_pairs'] = Speci_P
        self.Static_dict['Speci_Maternal_both_sides_pairs'] = Speci_M_both
        self.Static_dict['Speci_Paternal_both_sides_pairs'] = Speci_P_both
        self.Static_dict['Speci_Maternal_single_sides_pairs'] = Speci_M_single
        self.Static_dict['Speci_Paternal_single_sides_pairs'] = Speci_P_single
        self.Static_dict['Recombination_pairs'] = Regroup
        self.Static_dict['Allelic_Ratio'] = float(Both_M+Both_P+Single_M+Single_P) / (count-1)
        
        Statics_out = ['Total pairs number is %d' % (count - 1),
                       'Bi-Allelic pairs number is %d' % Bi_Allelic,
                       'Maternal Allelic pairs number is %d' % (Both_M+Single_M),
                       'Paternal Allelic pairs number is %d' % (Both_P+Single_P),
                       'Both sides for Maternal number is %d' % Both_M,
                       'Both sides for Paternal number is %d' % Both_P,
                       'Single sides for Maternal number is %d' % Single_M,
                       'Single sides for Paternal number is %d' % Single_P,
                       'Recombination number is %d' % Regroup,
                       '\n',
                       '==========Specifical Mapping Statics==============',
                       
                       'Specifical Mapping pairs for Maternal number is %d' % Speci_M,
                       '        Both sides for Speci Maternal number is %d' % Speci_M_both,
                       '        Single sides for Speci Maternal number is %d' % Speci_M_single,
                       'Specifical Mapping pairs for Paternal number is %d' % Speci_P,
                       '        Both sides for Speci Paternal number is %d' % Speci_P_both,
                       '        Single sides for Speci Paternal number is %d' % Speci_P_single,
                       '\n',
                       'Allelic Ratio : %.4f' % self.Static_dict['Allelic_Ratio']]
        
        Statics_out = '\n'*2 + '\n'.join(Statics_out) + '\n'
        log.log(21,'Allelic Filtering Results : %s',Statics_out)
        #print 'Allelic Filtering Results : %s' % Statics_out
        
                        
                
                        