#!/bin/env python
import os, random
from StringIO import StringIO
from Bio import SeqIO
from Bio import bgzf
import h5py_cache

def get_content(ff):
    Key4 = 'Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
    for f in ff.split(' '):
        with h5py_cache.File(f, 'r', chunk_cache_mem_size=1024**2) as f5:
            try:
                Fqstring = SeqIO.read(StringIO(f5[Key4][()]), "fastq")
                yield Fqstring
            except:
                pass    

def parsing_func(datas, Norder, t_dir):             
    file_out = os.path.join(t_dir, 'tmp.' + str(random.randint(1,21000000)) + '.fastq.gz')
    Gzout = get_content(datas)
    with bgzf.BgzfWriter(file_out, "wb") as outgz:
        SeqIO.write(sequences=Gzout, handle=outgz, format="fastq")
    return file_out