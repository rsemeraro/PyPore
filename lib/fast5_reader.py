#!/bin/env python
import os
import sys
import re
import time
from StringIO import StringIO
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio import SeqIO, bgzf
import h5py_cache
from dateutil import parser as par

global FastQFlag

FastQFlag = False

def Pathcheck(f5, fn):
    Visits = []
    DataPaths = []
    f5.visit(Visits.append)
    if 'fail' not in fn:
        Paths = ['Raw/Reads', 'UniqueGlobalKey/channel_id', 'UniqueGlobalKey/tracking_id',
                 'Analyses/Basecall_1D_000/Summary/basecall_1d_template',
                 'Analyses/Basecall_1D_000/BaseCalled_template/Fastq']
    else:
        Paths = ['Raw/Reads', 'UniqueGlobalKey/channel_id', 'UniqueGlobalKey/tracking_id']
    for i in Paths:
        try:
            Visits.index(i)
            DataPaths.append(i)
        except Exception:
            log.warn('%s is not in the %s' % (i, f5))
            continue
    return DataPaths

def get_content(f, fqf):
    QScore = 6.99
    FqInfos = None
    ReadLength = 0
    GCcont = 0
    #log.debug("[Fast5Parser] - Reading file: %s \r" % (f))
    ReadNumber = int(re.sub(r'.*read', 'read', f).split('_')[1])
    with h5py_cache.File(f, 'r', chunk_cache_mem_size=1024**2) as f5:
        DP = Pathcheck(f5, f)
        Key0 = DP[0] + '/Read_' + str(ReadNumber)
        Key1 = DP[1]
        Key2 = DP[2]
        Mux = f5[Key0].attrs['start_mux']
        IdRead = f5[Key0].attrs['read_id']
        Channel = f5[Key1].attrs['channel_number']
        StartTime = f5['Analyses/Segmentation_000'].attrs['time_stamp'] 
        if 'fail' not in f:
            Key3 = DP[3]
            Key4 = DP[4]
            ReadLength = f5[Key3].attrs['sequence_length']
            QScore = f5[Key3].attrs['mean_qscore']
            Fqstring = SeqIO.read(StringIO(f5[Key4][()]), "fastq")
            GCcont = GC(str(Fqstring.seq))
            if fqf == True:
                Fqstring = None
                FqInfos = f5[Key4][()]
            else:
                Fqstring = None
                FqInfos = None       
        f5.flush()       
    dt = par.parse(StartTime)
    StartTConv = time.mktime(dt.timetuple())
    ReturnList = []
    ReturnLine = str(Channel), str(Mux), str(StartTConv), str(ReadLength), str(QScore), str(GCcont)
    ReturnList.append(ReturnLine)
    ReturnList.append(FqInfos)
    return ReturnList


def yielder(newgen):
    for ls in newgen:
        sls = SeqIO.read(StringIO(ls), "fastq")
        yield sls

def parsing_func(datas, Norder, t_dir, flag):
    FastQFlag = flag
    datas = datas.split(' ')
    read_lists = []
    MinTime=time.time()
    ChrDMux = {str(el): {str(il): [] for il in range(1, 5)} for el in range(1, 513)}
    for ds in datas:
        try:
            res = get_content(ds, FastQFlag)
            ch = res[0][0]
            mu = res[0][1]
            ChrDMux[ch][mu] = ChrDMux[ch].get(mu, []) + [res[0][2:]]        
            NewTime = float(res[0][2])
            if NewTime < MinTime:
                MinTime = NewTime
            if FastQFlag == True:    
                if float(res[0][4]) >= 7.0:               
                    read_lists.append(res[1])
        except:
            pass    
    return ChrDMux, MinTime, read_lists