#!/bin/env python
import os
import sys
import re
import time
from io import StringIO
from Bio.SeqUtils import GC
from Bio import SeqIO, bgzf
import h5py
from dateutil import parser as par


def Pathcheck(fn):
    Visits = []
    DataPaths = []
    with h5py.File(fn, 'r') as f5:
        f5.visit(Visits.append)
    if 'fail' not in fn:
        Paths = ['Raw/Reads', 'UniqueGlobalKey/channel_id', 'UniqueGlobalKey/tracking_id',
                 'Analyses/Basecall_1D_000/Summary/basecall_1d_template',
                 'Analyses/Basecall_1D_000/BaseCalled_template/Fastq']              
    else:
        Paths = ['Raw/Reads', 'UniqueGlobalKey/channel_id', 'UniqueGlobalKey/tracking_id']
    for ii,i in enumerate(Paths):
        try:
            Visits.index(i)
            DataPaths.append(i)
            if ii == 0:
                DataPaths[0] = Visits[Visits.index(i)]              
        except Exception:
            log.warn('%s is not in the %s' % (i, f5))
            continue
    return DataPaths


def get_content(f, fqf, DP):
    QScore = 6.99
    FqInfos = None
    ReadLength = 0
    GCcont = 0
    ReturnList = []                 
    ReadNumber = int(re.sub(r'.*read', 'read', f).split('_')[1])  
    with h5py.File(f, 'r') as f5:
        Key0b = DP[0]
        Key0 = '%(Key0b)s/Read_%(ReadNumber)s' % locals()
        Key1 = DP[1]
        Key2 = DP[2]
        Mux = f5[Key0].attrs['start_mux']
        IdRead = f5[Key0].attrs['read_id'].decode('UTF-8')
        Channel = f5[Key1].attrs['channel_number'].decode('UTF-8')
        StartTime = f5['Analyses/Segmentation_000'].attrs['time_stamp'].decode('UTF-8') 
        if 'fail' not in f:
            Key3 = DP[3]
            Key4 = DP[4]
            ReadLength = f5[Key3].attrs['sequence_length']
            QScore = f5[Key3].attrs['mean_qscore']
            Fqstring = f5[Key4][()].decode('UTF-8').split('\n')[1]
            GCcont = GC(Fqstring)
            if fqf == True:
                Fqstring = None
                FqInfos = f5[Key4][()].decode('UTF-8')
            else:
                Fqstring = None
                FqInfos = None       
        f5.flush()       
    dt = par.parse(StartTime)
    StartTConv = time.mktime(dt.timetuple())
    ReturnLine = str(Channel), str(Mux), str(StartTConv), str(ReadLength), str(QScore), str(GCcont)
    ReturnList.append(ReturnLine)
    ReturnList.append(FqInfos)
    return ReturnList


def yielder(newgen):
    for ls in newgen:
        sls = SeqIO.read(StringIO(ls), "fastq")
        yield sls


class parsing_func(object):
    def __init__(self, datas, Norder, t_dir, flag, RefStart):
        self.datas = datas
        self.Norder = Norder
        self.t_dir = t_dir
        self.flag = flag
        self.RefStart = RefStart

    def __call__(self):
        FastQFlag = self.flag
        datas = self.datas
        read_lists = []
        MinTime = time.time()
        FailedReads = 0
        ChrList = []
        try:
            file_path_check = next(i for i in datas if 'pass' in i)
        except:
            try:
                file_path_check = next(i for i in datas if 'fail' in i)
            except:
                file_path_check = datas[0]   
        p_check = Pathcheck(file_path_check)
        for ds in datas:
            try:
                res = get_content(ds, FastQFlag, p_check)
                ch = res[0][0]
                mu = res[0][1]
                ChrList.append(res[0][1:])      
                NewTime = float(res[0][2])
                if NewTime < MinTime:
                    MinTime = NewTime
                if FastQFlag == True:    
                    if float(res[0][4]) >= 7.0:               
                        read_lists.append(res[1])     
            except: 
                del ds
        ChrList[:] = [ChrList[z] for z in (y[0] for y in sorted(enumerate(zip(*ChrList)[1]), key=lambda z: z[1]))]
        TimeVec = map(lambda x: int(time.strftime('%H', time.localtime(float(x) - float(self.RefStart)))), zip(*ChrList)[1])
        for e in range(len(TimeVec)):
            hour = TimeVec[e]
            if e > 1:
                if hour < TimeVec[e - 1]:
                    TimeVec[e] = hour + 24
        ChrList = zip(TimeVec, list(zip(*ChrList)[2]), list(zip(*ChrList)[3]), list(zip(*ChrList)[0]), list(zip(*ChrList)[4]))
        ReadPerChannel = len(ChrList)
        BasesPerChannel = sum(map(int, zip(*ChrList)[1]))
        MuxProductivity = {str(il): [] for il in range(1, 5)}
        for k in MuxProductivity.keys():
            MuxProd = {}
            MuxBase = {}
            MucList = map(lambda y:ChrList[y], [i for i, x in enumerate(zip(*ChrList)[3]) if x == k])
            for hr, val, qual, muc, gcs in MucList:
                if float(qual) < 7.0:
                    FailedReads += 1
                if str(hr) in MuxProd:
                    MuxProd[str(hr)] += int(val)
                    MuxBase[str(hr)] += 1
                    
                else:
                    MuxProd[str(hr)] = int(val)
                    MuxBase[str(hr)] = 1        
            MuxProductivity[k].append(dict(zip(MuxBase.keys(), zip(MuxBase.values(), MuxProd.values()))))
        if FastQFlag == True:
            file_out = os.path.join(self.t_dir, 'tmp.' + str(self.Norder) + '.fastq.gz')
            Gzout = yielder(read_lists)
            with bgzf.BgzfWriter(file_out, "wb") as outgz:
                SeqIO.write(sequences=Gzout, handle=outgz, format="fastq")                                   
        ObjectOut = []
        ObjectOut.append(str(self.Norder))
        ObjectOut.append(str(ReadPerChannel))
        ObjectOut.append(str(BasesPerChannel))
        ObjectOut.append(str(FailedReads))
        ObjectOut.append(ChrList)
        ObjectOut.append(MuxProductivity) 
        return ObjectOut
