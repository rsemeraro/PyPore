from mpi4py import MPI
from subprocess import Popen, PIPE
import sys, os, shutil, ntpath, warnings
from numpy import sort
from logging_module import log

if __name__ == '__main__':

    global size
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if not len(sys.argv) > 4:
        sys.exit(1)   
    work_dir = sys.argv[1]
    file_folder = sys.argv[2]
    prefix = sys.argv[3]
    out_dir = sys.argv[4]
    ver = sys.argv[5]
    tmp_dir = os.path.join(work_dir, '.sts_temp')

    if rank == 0:
        #############################
        ### Define main functions ###
        #############################

        def verbosity(ver):
            if ver == 0:
                lib.logging_module.log.setLevel(lib.logging_module.logging.WARN)
            elif ver == 1:
                lib.logging_module.log.setLevel(lib.logging_module.logging.INFO)
            elif ver == 2:
                lib.logging_module.log.setLevel(lib.logging_module.logging.DEBUG)
            elif ver == 3:
                sys.stdout = open(os.devnull, 'w')
                lib.logging_module.log.setLevel(lib.logging_module.logging.DEBUG)
                lib.logging_module.ch.setLevel(lib.logging_module.logging.DEBUG)


        def dir_surfer():
            for subdir, dirs, files in os.walk(file_folder):
                for fi in files:
                    if fi.endswith('fast5') and 'fail' not in subdir:
                        f = os.path.join(subdir, fi)
                        yield f
                        
        def split_seq(seq, size):
            nbase = 100000
            if len(seq) < 100000:
                nbase = int(round(len(seq), -4))/15
            if size > 1:
                n = nbase/size
                list_base = [' '.join(seq[i * n:(i + 1) * n]) for i in range((len(seq) + n - 1) // n )]
            else:
                n = nbase
            list_base = [' '.join(seq[i * n:(i + 1) * n]) for i in range((len(seq) + n - 1) // n )]   
            newseq = []
            ls = len(list_base)
            splitsize = 1.0/size*ls
            for i in range(size):
                newseq.append(list_base[int(round(i*splitsize)):int(round((i+1)*splitsize))])  
            return newseq       
        
        def fastq_writer(tmp_dir):
            finalfile = os.path.join(work_dir, (prefix + '.fastq.gz'))
            file_list = os.listdir(tmp_dir)
            file_list.sort(key=lambda x: int(ntpath.basename(x)[4:-9]))
            with open(finalfile, 'w') as outfile:
    	        for infile in map(lambda x: os.path.join(tmp_dir, x), file_list):
        	        shutil.copyfileobj(open(infile), outfile)
            for b in map(lambda x: os.path.join(tmp_dir, x), file_list):
                try:
                    os.remove(b)
                except Exception:
                    log.warning("File %s doesn't exist!" % b)
            if not os.path.exists(os.path.join(work_dir, out_dir)):
                os.makedirs(os.path.join(work_dir, out_dir))            
            try:
                shutil.move(finalfile, os.path.join(work_dir, out_dir))
            except shutil.Error:
                log.error("Unable to move %s" % finalfile)
            shutil.rmtree(tmp_dir, ignore_errors=True)
            log.info('Finish')    
                               
    if rank == 0:
        log.info('Start Main Parser')
        verbosity(ver)
        file_l = list(dir_surfer())
        if not os.path.exists(os.path.join(tmp_dir)):
            os.makedirs(os.path.join(tmp_dir))   
        chunks = split_seq(file_l, size)
    else:
        file_l = None
        chunks  = None
         
    # ===================================
    # This should be node-related work
    # ===================================
    
    def mpi_caller(m_f, Norder):
        import fastqparser as fqr
        fqr.parsing_func(m_f, Norder, tmp_dir)
           
    # ===================================
    # End of node-related work
    # ===================================  
    results = [] 
    my_files = comm.scatter(chunks, root=0) 
    for job in enumerate(my_files):
          results.append(mpi_caller(job[1], job[0]))
    results = comm.gather(results, root=0)
    
    if rank == 0:
        fastq_writer(tmp_dir)
