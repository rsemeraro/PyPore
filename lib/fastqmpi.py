import multiprocessing
import sys, os, shutil, ntpath, warnings
from numpy import sort
from logging_module import log
warnings.filterwarnings('ignore', category=RuntimeWarning)

####################
## Define classes ##
####################


class Consumer(multiprocessing.Process):
    def __init__(self, task_queue, result_queue, main='.'):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.main = main

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                self.task_queue.task_done()
                break
            next_task.main = self.main
            next_task.sub = proc_name
            answer = next_task()
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return


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


def dir_surfer(skip_flag):
    for subdir, dirs, files in os.walk(file_folder):
        for fi in files:
            if skip_flag == False:
                if fi.endswith('fast5') and 'fail' not in subdir:
                    f = os.path.join(subdir, fi)
                    yield f
            else:
                if fi.endswith('fast5'):
                    f = os.path.join(subdir, fi)
                    yield f       
                

def fastq_writer(tmp_dir):
    finalfile = os.path.join(work_dir, (prefix + '.fastq.gz'))
    file_list = os.listdir(tmp_dir)
    file_list.sort(key=lambda x: int(ntpath.basename(x)[4:-9]))
    with open(finalfile, 'wb') as outfile:
        for infile in map(lambda x: os.path.join(tmp_dir, x), file_list):
            shutil.copyfileobj(open(infile, 'rb'), outfile, 1024*1024*10)
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

 

def run(arguments):   
    if not len(arguments) > 4:
        sys.exit(1)
    global work_dir, file_folder, prefix, out_dir  
    work_dir = arguments[0]
    file_folder = arguments[1]
    prefix = arguments[2]
    out_dir = arguments[3]
    ver = arguments[4]
    tmp_dir = os.path.join(work_dir, '.sts_temp')
    size = int(arguments[5])
    f_flag = eval(arguments[6])[0]

    fail_flag = False
    if f_flag in ['y', 'yes']:
        fail_flag = True

    import fastqparser as fqr
    
    log.info('Start Main Parser')
    verbosity(ver)

    lis=[[] for _ in range(512)]

    for e,x in enumerate(dir_surfer(fail_flag)):
        i = x.split('_')[-2]
        idx = int(i)-1
        lis[idx].append(x)

    if not os.path.exists(os.path.join(tmp_dir)):
        os.makedirs(os.path.join(tmp_dir))

    Q = multiprocessing.JoinableQueue()
    R = multiprocessing.Queue()

    process_consumers_number = [Consumer(Q, R)
                                for ix in range(size)]

    for w in process_consumers_number:
        w.start()


    for e,chunk in enumerate(lis):
        if chunk != []:
            Q.put(fqr.parsing_func(chunk, e, tmp_dir))
        else: log.warn('Channel %s inactive!' % (e))


    for i in range(size):
        Q.put(None)
    Q.join()

    while not Q.empty():
        Q.get()

    if not os.path.exists(os.path.join(work_dir, out_dir)):
        os.makedirs(os.path.join(work_dir, out_dir))
                           
    fastq_writer(tmp_dir)
