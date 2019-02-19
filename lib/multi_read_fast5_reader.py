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
    DataPaths = []
    if F_Flag == False:
        Paths = map(lambda x:fn + '/' + x,['Raw', 'channel_id',
                 'Analyses/Basecall_1D_000/Summary/basecall_1d_template',
                 'Analyses/Basecall_1D_000/BaseCalled_template/Fastq'])         
    else:
        Paths = map(lambda x:fn + '/' + x,['Raw', 'channel_id'])
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


def yielder(newgen):
    for ls in newgen:
        sls = SeqIO.read(StringIO(ls), "fastq")
        yield sls


class channel_parser(object):
    def __init__(self, datas, Norder, RefStart):
        self.datas = datas
        self.Norder = Norder
        self.RefStart = RefStart

    def __call__(self):
        FailedReads = 0
        l = self.datas
        l[:] = [l[z] for z in (y[0] for y in sorted(enumerate(zip(*l)[1]), key=lambda z: z[1]))]
        TimeVec = map(lambda x: int(time.strftime('%H', time.localtime(float(x) - float(self.RefStart)))), zip(*l)[1])
        for e in range(len(TimeVec)):
            hour = TimeVec[e]
            if e > 1:
                if hour < TimeVec[e - 1]:
                    TimeVec[e] = hour + 24
        l = zip(TimeVec, list(zip(*l)[2]),
                               list(zip(*l)[3]), 
                               list(zip(*l)[0]),
                               list(zip(*l)[4]))
        ReadPerChannel = len(l)
        BasesPerChannel = sum(map(int, zip(*l)[1]))
        MuxProductivity = {str(il): [] for il in range(1, 5)}
        for k in MuxProductivity.keys():
            MuxProd = {}
            MuxBase = {}
            if len(l) > 1:
                MucList = map(lambda y:l[y], [i for i, x in enumerate(zip(*l)[3]) if x == k])
            else:
                MucList = l       
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
        Outdata = [self.Norder, ReadPerChannel, BasesPerChannel, FailedReads, l, MuxProductivity]            
        return Outdata

class mf5_reader(object):
    def __init__(self, multi_fast5, Norder, t_dir, fqflag, Fa_Flag = False):
        self.multi_fast5 = multi_fast5
        self.Norder = Norder
        self.t_dir = t_dir
        self.fqflag = fqflag
        self.Fa_Flag = Fa_Flag
    
    def get_content(self, f, fqf):
        QScore = 6.99
        FqInfos = None
        ReadLength = 0
        GCcont = 0
        ReturnList = []
        DP = Pathcheck(str(f))
        Key0 = DP[0]
        Key1 = DP[1]
        Key2 = DP[2]
        Mux = f5[Key0].attrs['start_mux']
        Channel = f5[Key1].attrs['channel_number'].decode('UTF-8')
        StartTime = f5[str(f) + '/Analyses/Segmentation_000'].attrs['time_stamp'].decode('UTF-8') 
        if F_Flag == False:
            Key3 = DP[3]
            ReadLength = f5[Key2].attrs['sequence_length']
            QScore = f5[Key2].attrs['mean_qscore']
            Fqstring = f5[Key3][()].decode('UTF-8').split('\n')[1]
            GCcont = GC(Fqstring)
            if fqf == True:
                Fqstring = None
                FqInfos = f5[Key3][()].decode('UTF-8')
            else:
                Fqstring = None
                FqInfos = None          
        dt = par.parse(StartTime)
        StartTConv = time.mktime(dt.timetuple())
        ReturnLine = str(Channel), str(Mux), str(StartTConv), str(ReadLength), str(QScore), str(GCcont)
        ReturnList.append(ReturnLine)
        ReturnList.append(FqInfos)      
        return ReturnList

    def __call__(self):
        global F_Flag
        F_Flag = self.Fa_Flag
        FastQFlag = self.fqflag
        read_lists = []
        MinTime = time.time()
        FailedReads = 0
        ChrList = [[] for _ in range(512)]
        global Visits
        Visits = []
        ChannelDict = {str(il): () for il in range(1, 513)}
        Outdata = []
        global f5
        f5 = h5py.File(self.multi_fast5, 'r')
        f5.visit(Visits.append)
        reads_list_to_read=f5.keys()
        RefStart = time.mktime(par.parse(f5[str(reads_list_to_read[0]) + '/tracking_id'].attrs['exp_start_time'].decode('UTF-8')).timetuple())
        for r in reads_list_to_read:
            res = self.get_content(r, FastQFlag)
            Outdata.append(res)       
        if FastQFlag == True:
            file_out = os.path.join(self.t_dir, 'tmp.' + str(self.Norder) + '.fastq.gz')
            Gzout = yielder(read_lists)
            with bgzf.BgzfWriter(file_out, "wb") as outgz:
                SeqIO.write(sequences=Gzout, handle=outgz, format="fastq")
        Outdata.append(RefStart)        
        return Outdata
