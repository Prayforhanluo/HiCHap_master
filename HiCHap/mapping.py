# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 15:10:17 2018

@author: han-luo
"""
from __future__ import division
import os, subprocess, multiprocessing, logging, re, sys, time
import xml.etree.ElementTree as ET

log = logging.getLogger(__name__)


def Getchunks(fastq):
    """
        Get chunks files and id from chunk Folder
    """
    chunks =[]
    num = 0
    reg = re.compile(r"(?<=_chunk)\d+")
    for fil in os.listdir(fastq):
        try:
            tmp_num = int(reg.search(fil).group(0))
        except:
            continue
            #log.warning('Unexpected file %s ', fil)
        if tmp_num < num :
            pass
        else:
            num = tmp_num
        
        chunks.append(fil)
    
    Cell_prefix = fil.split('_chunk')[0]
    
    return chunks, num + 1, Cell_prefix
    

#-------------------------Global Mapping Functions----------------------------------
    
def WS_sub_mapping(bowtieIndex, bowtiePath, threads, fq, OutPath, task_count):
    
    """
        sub-mapping  for WS station
        
    Parameters
    ----------
    bowtieIndex : str
        bowtieIndex if you already build the genome Index.
    
    bowtiePath : str
        bowtie2 path.
    
    threads : int
        threads
    
    fq : str
        Sequencing data for any chunk of one mate.
    
    OutPath : str
        OutFolder for mapping results.
    
    task_count : int
        id of sub-mapping task
    """
    
    fq_predix = os.path.join(OutPath,os.path.split(fq)[-1].split('.')[0])
    genome_prefix = os.path.split(bowtieIndex)[-1]
    out_sam_name = fq_predix + '_' + genome_prefix + '.sam'
    out_sam_file = os.path.join(OutPath,out_sam_name)
    out_bam_file = out_sam_file.rstrip('.sam') + '.bam'
    prefix = out_sam_file.rstrip('.sam')
    
    mapcmd = [bowtiePath, '-x', bowtieIndex, '-p', str(threads), '-U', fq, '|', 
              'samtools', 'view', '-b','-S', '|', 
              'samtools', 'sort', '-n','-T', prefix, '-o', out_bam_file]
    #log.log(21,'%s',' '.join(mapcmd))
    map_res = subprocess.Popen(' '.join(mapcmd),shell=True,stderr=subprocess.PIPE,bufsize=-1)
    
    map_res.communicate()
    
#    out_bam_file = out_sam_file.rstrip('.sam') + '.bam'
#    prefix = out_sam_file.rstrip('.sam')
#    sam_precess = ['samtools', 'view', '-b','-S',out_sam_file, '|', 'samtools', 'sort', '-n','-T', prefix, '-o', out_bam_file]
#    #log.log(21,'%s',' '.join(sam_precess))
#    samtools_res = subprocess.Popen(' '.join(sam_precess),shell = True,stdout=subprocess.PIPE,
#                                    stderr=subprocess.PIPE,bufsize=-1)
#    samtools_res.communicate()
    
    log.log(21,'%s mapping to %s done! count_id : %d ',fq_predix,genome_prefix,task_count)
    


def WS_mapping(fastq, threads, bowtiePath, OutPath, bowtieIndex = None):
    """
        WS(work station) Mapping mode.
        submit the mapping task in automatic way of system.
    
    Parameters
    ----------
    
    fastq : str
        chunk fq.gz Folder created by update.
        
    threads : str
        Total number of threads for four sub-mapping tasks.

    bowtiePath : str
        bowtie2 Path
        
    bowtieIndex : list or tuple
        List or tuple for Index of Maternal & Paternal Genome
    
    OutPath : str
        sam OutPath.
               
    """
    if not (isinstance(bowtieIndex,list) or isinstance(bowtieIndex,tuple)):
        log.error('bowtieIndex arrangement error ...')
        log.error('Exit ...')
        sys.exit(1)
    elif len(bowtieIndex) == 2:
        log.log(21,'two Index were found, arrange for M/Paternal genome ...')
        Maternal_Index = bowtieIndex[0]
        Paternal_Index = bowtieIndex[1]
        log.log(21,'Maternal Genome Index is %s ',bowtieIndex[0])
        log.log(21,'Paternal Genome Index is %s ',bowtieIndex[1])
    
        sub_threads = threads // 4
        chunks, num, Cell = Getchunks(fastq)
        log.log(21,'%d chunk paris were found under %s',num,fastq)
        log.log(21,'total number of sub-mapping tasks is %d ', num*4)    
    
        mapping_Pool = multiprocessing.Pool(4)
        count = 0
        for fil in chunks:
            fq = os.path.join(fastq,fil)
            count += 1
            mapping_Pool.apply_async(WS_sub_mapping,args=(Maternal_Index,
                                                          bowtiePath,
                                                          sub_threads,
                                                          fq,
                                                          OutPath,
                                                          count))
            count += 1
            mapping_Pool.apply_async(WS_sub_mapping,args=(Paternal_Index,
                                                          bowtiePath,
                                                          sub_threads,
                                                          fq,
                                                          OutPath,
                                                          count))
            
        time.sleep(2)
        mapping_Pool.close()
        mapping_Pool.join()        
    elif len(bowtieIndex) == 1:
        log.log(21,'Only one Index was found, Non-Allelic HiC mapping ...')
        genome_Index = bowtieIndex[0]
        log.log(21,'Genome Index is %s ', genome_Index)
        
        sub_threads = threads // 4
        chunks, num, Cell = Getchunks(fastq)
        log.log(21,'%d chunk paris were found under %s',num,fastq)
        log.log(21,'total number of sub-mapping tasks is %d ', num*2)
        
        mapping_Pool = multiprocessing.Pool(4)
        count = 0
        for fil in chunks:
            fq = os.path.join(fastq,fil)
            count += 1
            mapping_Pool.apply_async(WS_sub_mapping, args = (genome_Index,
                                                             bowtiePath,
                                                             sub_threads,
                                                             fq,
                                                             OutPath,
                                                             count))
        time.sleep(2)
        mapping_Pool.close()
        mapping_Pool.join()
        
    else:
        log.error('bowtieIndex Number not correct, only one Index or two Index supported. ...')
        log.error('Exit')
        sys.exit(1)
    
    log.log(21,'Multiprocess mapping Done!'
                'Try to continue with filtering task,haphic filtering -h for help information')
    
    
 
def PBS_sub_mapping(fq, bowtieIndex, bowtiePath, threads, OutPath, logPath, mem,task_count):
    """
        PBS scripts for sub mapping.
    
    Parameters
    ----------
    
    fq : str
        Sequencing data for any chunk of one mate.
    
    bowtieIndex : str
        Genome Index Path if you already build the genome Index.
    
    bowtiePath : str
        bowtie2 Path.
        
    threads : int
        threads.
    
    OutPath : str
        OutFolder for mapping results.
        
    logPath : str
        OutFolder for mapping log.
    
    mem : int
        memory  resource for each task        
                    
    task_count : int
        id of sub-mapping task
        
    """
    cell = os.path.split(fq)[-1].split('_chunk')[0]
    fq_prefix = os.path.split(fq)[-1].split('.')[0]
    genome_prefix = os.path.split(bowtieIndex)[-1]
    out_bam_name = fq_prefix + '_' + genome_prefix + '.bam'
    prefix = out_bam_name.rstrip('.bam')
    out_bam_file = os.path.join(OutPath,out_bam_name)
    
    PBS_scripts = ('echo "bowtie2 -x %s -U %s -p %s  | samtools view -b -S | '
                    ' samtools sort -n -T %s -@ %s -o %s" '
                    '| qsub -N %s '
                    '-l nodes=1:ppn=%s -l mem=%sgb -d ./ '
                    '-e %s -o %s'
                    %(bowtieIndex,fq,threads,prefix,threads,out_bam_file,
                      cell ,threads,mem,logPath,logPath))
    
    subprocess.Popen(PBS_scripts,shell=True,stdout=subprocess.PIPE)
    
    log.log(21,'%s To %s has been submitted. count_id : %d ',fq_prefix, genome_prefix, task_count)
    
#def PBS_Pooler(keyword):
#    """
#        Check the tasks state in Pool.
#    
#    Parameters
#    ----------
#    keyword : str
#        key word of task name. 
#    
#    return the number of task which has key word.      
#    """
#    
#    caller = subprocess.check_output(['qstat'])
#    
#    state = caller.split('\n')[2:-1]
#    
#    num = 0
#    for i in state:
#        i = i.split()
#        if keyword in i[1]:
#            num += 1
#        else:
#            pass
#    
#    return num

def PBS_Pooler(keyword):
    """
        Check the tasks state in Pool.
    
    Parameters
    ----------
    keyword : str
        key word of task name.
    
    return the number of task which has keyword.
    """
    
    jobs_XML = os.popen('qstat -xl') # check the jobs full information as XML format
    job_num = 0
    try:
        tree = ET.parse(jobs_XML)
    except:
        return 0
    root = tree.getroot()
    
    for sub_job in root.getchildren():
        tmp_name = sub_job.find('Job_Name').text
        if keyword in tmp_name:
            job_num += 1
        else:
            continue
    return job_num



def Pool_controller(keyword, num_task):
    """
        controller of Whether submit the new task.
    """
    while True:
        if PBS_Pooler(keyword) < num_task:
            return True
        else:
            time.sleep(5)

def check_Non_Allelic_out(path,chunk_num):
    """
        check out whether the mapping output are completed.
        we will check whether there are 2 mates mapping sam/bam for a chunk id.
        and whether the size of mapping results are zero.
        
    Parameters
    ----------
    
    path : str
        mapping out path
    
    chunk_num : int
        The number of chunks.    
        
    """
    files = os.listdir(path)
    error_list = []
    for i in range(chunk_num):
        tmp = []
        for fil in files:
            if 'chunk'+str(i)+'_' in fil:
                tmp.append(fil)
        if len(tmp) != 2:
            R1 = 'chunk'+str(i)+'_1'
            R2 = 'chunk'+str(i)+'_2'
            R1_tmp = False
            R2_tmp = False
            for fil in tmp:
                if R1 in fil:
                    R1_tmp = True
                if R2 in fil:
                    R2_tmp = True
            if not R1_tmp:
                error_list.append((i,'1'))
                log.warning('chunk%s_1 mapping results is lost...',i)
            if not R2_tmp:
                error_list.append((i,'2'))
                log.warning('chunk%s_2 mapping results is lost...',i)
    
    for fil in files:
        if os.path.getsize(os.path.join(path,fil)) < 100:
            tmp = fil.split('chunk')[1].replace('.bam','').split('_')[:-1]
            error_list.append(tuple(tmp))
            log.warning('mapping results %s is empty',fil)

    return error_list
    
    
def check_Allelic_out(path,chunk_num):
    """
        check out whether the mapping output are completed.
        We will check whether there are 4 Allelic mapping sam/bam for a chunk id.
        and whether the size of mapping results are zero.
    
    Parameters
    ----------
    
    path : str
        mapping out path
    
    chunk_num : int
        The number of chunks.
        
    """
    files = os.listdir(path)
    error_list = []
    for i in range(chunk_num):
        tmp = []
        for fil in files:
            if 'chunk'+str(i)+'_' in fil:
                tmp.append(fil)
        if len(tmp) != 4:
            M_1 = 'chunk'+str(i)+'_1_Maternal'
            M_2 = 'chunk'+str(i)+'_2_Maternal'
            P_1 = 'chunk'+str(i)+'_1_Paternal'
            P_2 = 'chunk'+str(i)+'_2_Paternal'
            M_1_tmp = False
            M_2_tmp = False
            P_1_tmp = False
            P_2_tmp = False
            for fil in tmp:
                if M_1 in fil:
                    M_1_tmp = True
                if M_2 in fil:
                    M_2_tmp = True
                if P_1 in fil:
                    P_1_tmp = True
                if P_2 in fil:
                    P_2_tmp = True
            if not M_1_tmp:
                error_list.append((i,'1','Maternal'))
                log.warning('chunk%s_1 mapping To Maternal is lost...',i)
            if not M_2_tmp:
                error_list.append((i,'2','Maternal'))
                log.warning('chunk%s_2 mapping To Maternal is lost...',i)
            if not P_1_tmp:
                error_list.append((i,'1','Paternal'))
                log.warning('chunk%s_1 mapping To Paternal is lost...',i)
            if not P_2_tmp:
                error_list.append((i,'2','Paternal'))
                log.warning('chunk%s_2 mapping To Paternal is lost...',i)
            
    
    for fil in files:
        if os.path.getsize(os.path.join(path,fil)) < 100:
            tmp = fil.split('chunk')[1].rstrip('.bam').split('_')
            error_list.append(tuple(tmp))
            log.error('mapping results %s is empty',fil)
      
    return error_list
    

def PBS_controller(fastq, num_task, threads, OutPath, logPath, bowtiePath, mem, bowtieIndex = None):   
    """
        PBS tasks controller.
    Parameters
    ----------
    fastq : str
        chunk fq.gz Folder created by update.
        
    num_task : int
        Total sub-tasks number in task Pools.
    
    threads : int
        number of threads for a sub mapping task.
    
    OutPath : str
        OutFolder for mapping results.
    
    logPath : str
        OutFolder for mapping log.
    
    bowtiePath : str
        bowtie2 Path

    mem : int
        memory for each sub mapping tasks.        
        
    bowtieIndex : list or tuple
        List or tuple for Index of Maternal & Paternal Genome.
    
    """
    if not (isinstance(bowtieIndex,list) or isinstance(bowtieIndex,tuple)):
        log.error('Allelic bowtieIndex arrangement error ...')
        log.error('Exit ...')
        sys.exit(1)
    elif len(bowtieIndex) == 2:
        log.log(21,'two Index were found, arrange for M/Paternal genome ...')
        Maternal_Index = bowtieIndex[0]
        Paternal_Index = bowtieIndex[1]
        log.log(21,'Maternal Genome Index is %s ',bowtieIndex[0])
        log.log(21,'Paternal Genome Index is %s ',bowtieIndex[1])    
    
        chunks, chunk_num, Cell = Getchunks(fastq)
        log.log(21,'%d chunk paris were found under %s',chunk_num,fastq)
        log.log(21,'total number of sub-mapping tasks is %d ', chunk_num * 4)    
    
        count = 0
        for fil in chunks:
            fq = os.path.join(fastq,fil)
            count += 1
            if Pool_controller(Cell, num_task):
                PBS_sub_mapping(fq=fq,
                                bowtieIndex = Maternal_Index,
                                bowtiePath= bowtiePath,
                                threads = threads,
                                OutPath = OutPath,
                                logPath=logPath,
                                mem = mem,
                                task_count= count)
            time.sleep(2)
            count += 1        
            if Pool_controller(Cell, num_task):
                PBS_sub_mapping(fq=fq,
                                bowtieIndex = Paternal_Index,
                                bowtiePath= bowtiePath,
                                threads = threads,
                                OutPath = OutPath,
                                logPath=logPath,
                                mem = mem,
                                task_count= count)
            time.sleep(2)
    
        while True:
            if PBS_Pooler(Cell) == 0:
                logging.log(21,'all chunks have been mapped To the Genome.')
                break
            else:
                time.sleep(2)
    
        # re-building mapping task for error mapping results.
        while True:
            error_list = check_Allelic_out(OutPath,chunk_num)
            if len(error_list) == 0:
                log.log(21,'all chunks mapping results are completed.')
                break
            else:
                log.log(21,'rebuilding the mapping for escaped mapping tasks')
                count = 0
                for i in error_list:
                    mask = 'chunk'+str(i[0])+'_'+str(i[1])
                    fq = [fil for fil in chunks if mask in fil][0]
                    fq = os.path.join(fastq,fq)
                    count += 1
                    if i[2] == 'Maternal':
                        Index = Maternal_Index
                    else:
                        Index = Paternal_Index
            
                    if Pool_controller(Cell, num_task):
                        PBS_sub_mapping(fq=fq,
                                        bowtieIndex = Index,
                                        bowtiePath= bowtiePath,
                                        threads = threads,
                                        OutPath = OutPath,
                                        logPath=logPath,
                                        mem = mem,
                                        task_count= count)
                    time.sleep(2)
        
            log.log(21,'all rebuilding tasks have submitted.')
    
            while True:
                if PBS_Pooler(Cell) == 0:
                    logging.log(21,'all rebuilding tasks have been mapped To the Genome.')
                    break
                else:
                    time.sleep(5)
    elif len(bowtieIndex) == 1:
        log.log(21,'Only one Index was found, Non-Allelic HiC mapping ...')
        genome_Index = bowtieIndex[0]
        log.log(21,'Genome Index is %s ', genome_Index)
        chunks, chunk_num, Cell = Getchunks(fastq)
        log.log(21,'%d chunk paris were found under %s',chunk_num,fastq)
        log.log(21,'total number of sub-mapping tasks is %d ', chunk_num*2)
        
        count = 0
        for fil in chunks:
            fq = os.path.join(fastq,fil)
            count += 1
            if Pool_controller(Cell, num_task):
                PBS_sub_mapping(fq = fq,
                                bowtieIndex=genome_Index,
                                bowtiePath=bowtiePath,
                                threads=threads,
                                OutPath=OutPath,
                                logPath=logPath,
                                mem=mem,
                                task_count=count)
            time.sleep(2)
        while True:
            if PBS_Pooler(Cell) == 0:
                logging.log(21,'all chunks have been mapped To the Genome.')
                break
            else:
                time.sleep(5)
        
        # re-building mapping task for error mapping results
        while True:
            error_list = check_Non_Allelic_out(OutPath,chunk_num)
            if len(error_list) == 0:
                log.log(21,'all chunks mapping results are completed.')
                break
            else:
                log.log(21,'rebuilding the mapping for escaped mapping tasks')
                count = 0
                for i in error_list:
                    mask = 'chunk'+str(i[0])+'_'+str(i[1])
                    fq = [fil for fil in chunks if mask in fil][0]
                    fq = os.path.join(fastq,fq)
                    count += 1
                    if Pool_controller(Cell, num_task):
                        PBS_sub_mapping(fq = fq,
                                        bowtieIndex = genome_Index,
                                        bowtiePath = bowtiePath,
                                        threads = threads,
                                        OutPath = OutPath,
                                        logPath = logPath,
                                        mem = mem,
                                        task_count = count)
                    time.sleep(2)
            
            log.log(21,'all rebuilding tasks have submitted.')
            
            while True:
                if PBS_Pooler(Cell) == 0:
                    logging.log(21,'all rebuilding tasks have been mapped To the Genome.')
                    break
                else:
                    time.sleep(5)
                                        
    else:
        log.error('bowtieIndex Number not correct, only one Index or two Index supported. ...')
        log.error('Exit')
        sys.exit(1)

    
#----------------------------------Rescue Mapping ----------------------------------

def Rescue_sub_WS_mapping(bowtieIndex, bowtiePath, threads, fq, OutPath):
    """
        sub-rescue-mapping for WS station.
    
    Parameters
    ----------
    bowtieIndex : str
        bowtieIndex if you already build the genome Index.
    
    bowtiePath : str
        bowtie2 Path.
    
    threads : int
        threads
    
    fq : str
        Sequencing data for any chunk of one mate.
    
    OutPath : str
        OutFolder for mapping results
        
    """
    fq_prefix = os.path.join(OutPath,os.path.split(fq)[-1].split('.')[0])
    out_bam_file = os.path.join(OutPath,fq_prefix+'.bam')
    
    mapcmd = [bowtiePath, '-x', bowtieIndex, '-p', str(threads), '-U', fq, '|', 
              'samtools', 'view', '-b','-S', '|', 
              'samtools', 'sort', '-n','-T', fq_prefix, '-o', out_bam_file] 
              
    map_res = subprocess.Popen(' '.join(mapcmd), shell=True,stderr=subprocess.PIPE,bufsize=-1)
    
    #print ' '.join(mapcmd)
    map_res.communicate()
    

 
def Rescue_WS_mapping(fastq,threads,bowtiePath,OutPath,bowtieIndex = None):
    """
        Resuce mapping for WS station.
        
    """
    if not (isinstance(bowtieIndex,list) or isinstance(bowtieIndex,tuple)):
        log.error('BowtieIndex arrangement error ...')
        log.error('Exit ...')
        sys.exit(1)
    elif len(bowtieIndex) == 2:
        log.log(21,'two Index were found, arrange for M/Paternal genome ...')
        log.log(21,'Allelic Rescue Mapping start ...')
        Maternal_Index = bowtieIndex[0]
        Paternal_Index = bowtieIndex[1]
        log.log(21,'Maternal Genome Index is %s', bowtieIndex[0])
        log.log(21,'Paternal Genome Index is %s', bowtieIndex[1])
        
        sub_threads = threads // 4
        Maternal_chunks = [i for i in os.listdir(fastq) if 'Maternal' in i]
        Paternal_chunks = [i for i in os.listdir(fastq) if 'Paternal' in i]
        
        log.log(21,'Rescue Mapping ...')
        mapping_Pool = multiprocessing.Pool(4)
        for fil in Maternal_chunks:
            fq = os.path.join(fastq,fil)
            
            mapping_Pool.apply_async(Rescue_sub_WS_mapping,args=(Maternal_Index,
                                                                 bowtiePath,
                                                                 sub_threads,
                                                                 fq,
                                                                 OutPath))
        time.sleep(1)
        for fil in Paternal_chunks:
            fq = os.path.join(fastq,fil)
            mapping_Pool.apply_async(Rescue_sub_WS_mapping,args=(Paternal_Index,
                                                                 bowtiePath,
                                                                 sub_threads,
                                                                 fq,
                                                                 OutPath))
        time.sleep(1)
        mapping_Pool.close()
        mapping_Pool.join()
        log.log(21,'Rescue Mapping Done!')
        
    elif len(bowtieIndex) == 1:
        log.log(21,'Only one Index was found, Non-Allelic Rescue Mapping start...')
        genome_Index = bowtieIndex[0]
        log.log(21,'Genome Index is %s ', genome_Index)
        
        sub_threads = threads // 4
        fastq_chunks = [i for i in os.listdir(fastq) if 'chunk' in i]

        log.log(21,'Rescue Mapping')
        mapping_Pool = multiprocessing.Pool(4)
        for fil in fastq_chunks:
            fq = os.path.join(fastq,fil)
            mapping_Pool.apply_async(Rescue_sub_WS_mapping,args=(genome_Index,
                                                                 bowtiePath,
                                                                 sub_threads,
                                                                 fq,
                                                                 OutPath))
        time.sleep(1)
        mapping_Pool.close()
        mapping_Pool.join()
        log.log(21,'Rescue Mapping Done!')
    else:
        log.error('bowtieIndex Number not correct, only one Index or two Index supported. ...')
        log.error('Exit ...')
        sys.exit(1)


def Rescue_check_out(fq,OutPath):
    """
        Check Out Whether the ReMapping results.  
    """
    errorlist = []
    Info_fq = os.listdir(fq)
    Info_bam = os.listdir(OutPath)
    for i in Info_fq:
        tmp = i.split('.')[0]+'.bam'
        if tmp in Info_bam:
            if os.path.getsize(os.path.join(OutPath,tmp)) < 100:
                error = i.split('chunk')[-1].split('_')[:3]
                errorlist.append(error)
                log.warning('chunk%s_%s to %s ReMapping result is empty ...',error[0],error[1],error[2])
                #print 'chunk%s_%s to %s ReMapping result is empty ...' % (error[0],error[1],error[2])
            else:
                pass
        else:
            error = i.split('chunk')[-1].split('_')[:3]
            errorlist.append(error)
            log.warning('chunk%s_%s to %s ReMapping result is lost ...',error[0],error[1],error[2])
            #print 'chunk%s_%s to %s ReMapping result is lost ...'% (error[0],error[1],error[2])
    return errorlist





def Rescue_sub_PBS_mapping(fq,bowtieIndex,bowtiePath,threads,OutPath,logPath,mem,task_count):
    """
        PBS scripts for sub Rescue mapping.
        
    Parameters
    ----------
    fq : str
        Sequencing data for any chunk of one mate.
    
    bowtieIndex : str
        Genome Index Path if you already build the genome Index.
    
    bowtiePath : str
        bowtie2 Path.
        
    threads : int
        threads.
    
    OutPath : str
        OutFolder for mapping results.
        
    logPath : str
        OutFolder for mapping log.
    
    mem : int
        memory  resource for each task        
                    
    task_count : int
        id of sub-mapping task
        
    """
    cell = os.path.split(fq)[-1].split('_chunk')[0]
    fq_prefix = os.path.split(fq)[-1].split('.')[0]
    out_bam_file = os.path.join(OutPath,fq_prefix+'.bam')
    
    PBS_scripts = ('echo "bowtie2 -x %s -U %s -p %s  | samtools view -b -S | '
                    ' samtools sort -n -T %s -@ %s -o %s" '
                    '| qsub -N %s '
                    '-l nodes=1:ppn=%s -l mem=%sgb -d ./ '
                    '-e %s -o %s'
                    %(bowtieIndex,fq,threads,fq_prefix,threads,out_bam_file,
                      cell ,threads,mem,logPath,logPath))
    
    subprocess.Popen(PBS_scripts,shell=True,stdout=subprocess.PIPE)
    log.log(21,'%s  is remapping . count_id : %d ',fq_prefix, task_count)
    

def Rescue_PBS_controller(fastq,num_task, threads, OutPath, logPath, bowtiePath, mem, bowtieIndex = None):
    """
       Rescue PBS tasks controller.
    Parameters
    ----------
    fastq : str
        chunk fq.gz Folder created by update.
        
    num_task : int
        Total sub-tasks number in task Pools.
    
    threads : int
        number of threads for a sub mapping task.
    
    OutPath : str
        OutFolder for mapping results.
    
    logPath : str
        OutFolder for mapping log.
    
    bowtiePath : str
        bowtie2 Path

    mem : int
        memory for each sub mapping tasks.        
        
    bowtieIndex : list or tuple
        List or tuple for Index of  Genome.    
        
    """
    
    if not (isinstance(bowtieIndex,list) or isinstance(bowtieIndex,tuple)):
        log.error('bowtieIndex arrangement error')
        log.error('Exit ...')
        sys.exit(1)
    elif len(bowtieIndex) == 2:
        log.log(21,'two Index were found, arrange for M/Paternal genome ...')
        Maternal_Index = bowtieIndex[0]
        Paternal_Index = bowtieIndex[1]
        log.log(21,'Maternal Genome Index is %s', bowtieIndex[0])
        log.log(21,'Paternal Genome Index is %s', bowtieIndex[1])
        
        count = 0
        Maternal_chunks = [i for i in os.listdir(fastq) if 'Maternal' in i]
        Paternal_chunks = [i for i in os.listdir(fastq) if 'Paternal' in i]
        Cell = Maternal_chunks[0].split('_chunk')[0]
        for fil in Maternal_chunks:
            fq = os.path.join(fastq,fil)
            count += 1
            if Pool_controller(Cell, num_task):
                Rescue_sub_PBS_mapping(fq=fq,
                                       bowtieIndex=Maternal_Index,
                                       bowtiePath=bowtiePath,
                                       threads=threads,
                                       OutPath=OutPath,
                                       logPath=logPath,
                                       mem=mem,
                                       task_count=count)
            time.sleep(1)
        
        for fil in Paternal_chunks:
            fq = os.path.join(fastq,fil)
            count += 1
            if Pool_controller(Cell,num_task):
                Rescue_sub_PBS_mapping(fq=fq,
                                       bowtieIndex=Paternal_Index,
                                       bowtiePath=bowtiePath,
                                       threads=threads,
                                       OutPath=OutPath,
                                       logPath=logPath,
                                       mem=mem,
                                       task_count=count)
            time.sleep(1)
        
        while True:
            if PBS_Pooler(Cell) == 0:
                log.log(21,'Rescue Mapping have Done,Check the reults ...')
                break
            else:
                time.sleep(5)
        
        #re-building ReMapping task for error ReMapping results
        while True:
            error_list = Rescue_check_out(fastq,OutPath)
            if len(error_list) == 0:
                log.log(21,'all chunks mapping results are completed.')
                break
            else:
                log.log(21,'rebuilding the mapping for escaped mapping tasks')
                count = 0
                for i in error_list:
                    mask = 'chunk'+str(i[0])+'_'+str(i[1])
                    if i[2] == 'Maternal':
                        fq = [fil for fil in Maternal_chunks if mask in fil][0]
                        fq = os.path.join(fastq,fq)                        
                        bowtieIndex = Maternal_Index
                    else:
                        fq = [fil for fil in Paternal_chunks if mask in fil][0]
                        fq = os.path.join(fastq,fq)
                        bowtieIndex = Paternal_Index
                    count += 1
                    if Pool_controller(Cell, num_task):
                        Rescue_sub_PBS_mapping(fq = fq,
                                               bowtieIndex = bowtieIndex,
                                               bowtiePath = bowtiePath,
                                               threads = threads,
                                               OutPath = OutPath,
                                               logPath = logPath,
                                               mem = mem,
                                               task_count = count)
                    time.sleep(2)
            log.log(21,'all rebuilding tasks have submitted')
            
            while True:
                if PBS_Pooler(Cell) == 0:
                    logging.log(21,'all rebuilding tasks have been mapped To the Genome.')
                    break
                else:
                    time.sleep(5)
    
    elif len(bowtieIndex) == 1:
        log.log(21,'Only one Index was found, Non-Allelic Rescue Mapping start...')
        genome_Index = bowtieIndex[0]
        log.log(21,'Genome Index is %s', genome_Index)
        fastq_chunks = [i for i in os.listdir(fastq) if 'chunk' in i]
        Cell = fastq_chunks[0].split('_chunk')[0]
        
        count = 0
        for fil in fastq_chunks:
            fq = os.path.join(fastq,fil)
            count += 1
            if Pool_controller(Cell,num_task):
                Rescue_sub_PBS_mapping(fq=fq,
                                       bowtieIndex=genome_Index,
                                       bowtiePath=bowtiePath,
                                       threads=threads,
                                       OutPath=OutPath,
                                       logPath=logPath,
                                       mem=mem,
                                       task_count=count)
            time.sleep(1)
        
        while True:
            if PBS_Pooler(Cell) == 0:
                log.log(21,'Rescue Mapping have Done,Check the results')
                break
            else:
                time.sleep(5)
        #rebuilding the ReMapping task for error task
        while True:
            error_list = Rescue_check_out(fastq,OutPath)
            if len(error_list) == 0:
                log.log(21,'all chunks mapping results are completed.')
                break
            else:
                log.log(21,'rebuilding the mapping for escaped mapping tasks')
                count = 0
                for i in error_list:
                    mask = 'chunk'+str(i[0])+'_'+str(i[1])
                    fq = [fil for fil in fastq_chunks if mask in fil][0]
                    fq = os.path.join(fastq,fq)
                    bowtieIndex = genome_Index
                    count += 1
                    if Pool_controller(Cell, num_task):
                        Rescue_sub_PBS_mapping(fq = fq,
                                               bowtieIndex = bowtieIndex,
                                               bowtiePath = bowtiePath,
                                               threads = threads,
                                               OutPath = OutPath,
                                               logPath = logPath,
                                               mem = mem,
                                               task_count = count)
                    time.sleep(2)
            log.log(21,'all rebuilding tasks have submitted')
            
            while True:
                if PBS_Pooler(Cell) == 0:
                    logging.log(21,'all rebuilding tasks have been mapped To the Genome.')
                    break
                else:
                    time.sleep(5)

