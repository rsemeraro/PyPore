#!/bin/env python
import os, random
from io import StringIO
from Bio import SeqIO, bgzf
import h5py
from logging_module import log
def get_content(ff):
    Key4 = 'Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
    for f in ff:
        with h5py.File(f, 'r') as f5:
            try:
                Fqstring = SeqIO.read(StringIO(f5[Key4][()].decode('UTF-8')), "fastq")
                yield Fqstring
            except:
                log.warning("File %s seems doesn't contain any sequence!" % f)
                pass    
                

class parsing_func(object):
    def __init__(self, datas, Norder, t_dir):
        self.datas = datas
        self.Norder = Norder
        self.t_dir = t_dir

    def __call__(self):            
        file_out = os.path.join(self.t_dir, 'tmp.' + str(self.Norder) + '.fastq.gz')
        Gzout = get_content(self.datas)
        with bgzf.BgzfWriter(file_out, "wb") as outgz:
            SeqIO.write(sequences=Gzout, handle=outgz, format="fastq")
        return file_out
