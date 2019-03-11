#! /usr/bin/env python
import argparse
import sys
import tempfile
import atexit
import shutil
import os
from lib.logging_module import log

class readable_dir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        prospective_dir = values
        if not os.path.isdir(prospective_dir):
            log.error("readable_dir:{0} is not a directory".format(prospective_dir))
            sys.exit(1)
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace, self.dest, prospective_dir)
        else:
            log.error("readable_dir:{0} is not a directory".format(prospective_dir))
            sys.exit(1)


def ArgsReader(args):
    global work_dir
    global FileFolder
    global prefix
    global th
    global out_dir
    global ver
    work_dir = os.getcwd()
    FileFolder = args.input_directory
    prefix = args.label
    th = args.threads_number
    ver = args.verbose
    multi_read_flag = False
    fastq_flag = False
    summary_flag = False
    if args.fastq[0] in ['y', 'yes']:
        fastq_flag = True
    if args.output_dir is not None:
        out_dir = str(args.output_dir[0])
    else:
        out_dir = 'PyPore_results'
    albacore_file = args.albacore_summary
    if albacore_file != None:
        summary_flag = True
    if args.multi_read_fast5[0] in ['y', 'yes']:
        multi_read_flag = True    

    ags=[work_dir, FileFolder, prefix, out_dir, str(ver), str(fastq_flag), str(th), str(summary_flag), str(albacore_file), str(multi_read_flag)]
    import lib.seq_routines as sq
    sq.run(ags)

def run(argsin):
    ldir = tempfile.mkdtemp()
    atexit.register(lambda dir=ldir: shutil.rmtree(ldir))
    parser = argparse.ArgumentParser(add_help=False,
                                     prog='seqstats',
                                     usage='''\r      \n
                                  ________________
                                 |                |
                                 |    #####       |
                                 |    ##  ##      |
                                 |    #####       |
                                 |    ## #####    |
                                 |    ## ##  ##   |
                                 |       #####    |
                                 |       ##       |
                                 |       ##       |
                                 |________________|

                                       PyPore    
                                      SeqStats

usage: %(prog)s [-i <input_directory>] [-l <my_label>] [options]''',
                                 epilog='''

PyPore. Written by Roberto Semeraro, Department of Clinical and Sperimental Medicine, 
University of Florence. For bug report or suggestion write to robe.semeraro@gmail.com''',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     )
    g = parser.add_argument_group(title='mandatory arguments',
                                  description='''-i,  --input_directory                                      path to the file folder
-l,  --label                                                  label for output file
    ''')
    g.add_argument('-i', '--input_directory', action=readable_dir, default=ldir,
                   help=argparse.SUPPRESS)
    g.add_argument('-l', '--label', action='store', metavar='',
                   help=argparse.SUPPRESS)
    f = parser.add_argument_group(title='options',
                                  description='''-fq, --fastq                                                 generate fastq.gz [no]
-a,  --albacore_summary            if provided, it is used to gather the sequencing
                                   informations (fastq, not fast5, folder required)  
-m,  --multi_read_fast5                  activate the multi read fast5 reading [no]                                                              
-o,  --output_dir                  if omitted, generates a results directory in the 
                                                                   current position                                                                                               
-n,  --threads_number                              number of processing threads [1]
-v,  --verbose                     increase verbosity: 0 = warnings only, 1 = info, 
                                    2 = debug, 3 = debug on terminal. Default is no
                                                    verbosity. No number means info
-h,  --help                                         show this help message and exit
    ''')
    f.add_argument('-fq', '--fastq', action="store", nargs=1, metavar='',
                   choices=['y', 'n', 'yes', 'no'], default = ['no'], help=argparse.SUPPRESS)
    f.add_argument('-a', '--albacore_summary', action="store", nargs=1, metavar='',
                   default = None, help=argparse.SUPPRESS)
    f.add_argument('-m', '--multi_read_fast5', nargs=1, metavar='',
                   choices=['y', 'n', 'yes', 'no'], default = ['no'], help=argparse.SUPPRESS)                   
    f.add_argument('-o', '--output_dir', action="store", nargs=1, metavar='',
                   help=argparse.SUPPRESS)
    f.add_argument('-n', '--threads_number', action="store", type=int, default=1, metavar='',
                   help=argparse.SUPPRESS)
    f.add_argument('-h', '--help', action="help",
                   help=argparse.SUPPRESS)
    f.add_argument('-v', '--verbose', const=1, default=1, type=int, nargs="?",
                   help=argparse.SUPPRESS, choices=range(0, 4))

    try:
        args = parser.parse_args(argsin)
        if not args.input_directory or not args.label:
            parser.print_help()
            log.error('Missing mandatory arguments!')
            sys.exit(1)
        else:
            ArgsReader(args)
    except IOError as msg:
        parser.error(str(msg))
