#! /usr/bin/env python
import argparse
import os, sys
from subprocess import Popen
from logging_module import log

def MainReader(args):
    global work_dir
    global fast_Q_file
    global ref
    global stats_trigger
    global prefix
    global out_dir
    global als
    global th
    work_dir = os.getcwd()
    fast_Q_file = ' '.join(map(str, args.inputs))
    ref = args.reference
    stats_trigger = str(args.alignment_stats[0])
    if args.output_dir is not None:
        out_dir = str(args.output_dir[0])
    else:
        out_dir = 'PyPore_results'
    prefix = str(args.prefix)
    als = args.aligner
    th = args.threads_number
    ags=[work_dir, fast_Q_file, ref, stats_trigger, prefix, out_dir, als, str(th)]
    if os.name == 'nt':
        import lib.alg_routines_win as sq
    else:   
        import lib.alg_routines_unix as sq
    sq.run(ags)
        

######################
## Argument Parsing ##
######################

def run(argsin):
    parser = argparse.ArgumentParser(add_help=False,
                                     prog='alignment',
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
                                      Alignment
                                                
usage: %(prog)s [-i <input.fq>] [-r <reference.fa>] [-l <my_label>] [options]''',
                                     epilog='''
                                     
PyPore. Written by Roberto Semeraro, Department of Clinical and Sperimental Medicine,
University of Florence. For bug report or suggestion write to robe.semeraro@gmail.com''',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     )
    g = parser.add_argument_group(title='mandatory arguments',
                                  description='''-i, --inputs                                       path to the read file (query.fq)
-r, --reference                                     indexed reference sequence file
-l, --prefix                                              prefix for alignment file
    ''')
    g.add_argument('-i', '--inputs', action='store', metavar='',
                   nargs='*', help=argparse.SUPPRESS)
    g.add_argument('-r', '--reference', action='store', metavar='',
                   help=argparse.SUPPRESS)
    g.add_argument('-l', '--prefix', action='store', metavar='',
                   help=argparse.SUPPRESS)
    f = parser.add_argument_group(title='options',
                                  description='''-a, --aligner           aligment module (b = bwa, m = minimap2, n = ngmlr)  [m b n]
-s, --alignment_stats                                 generate alignment stats [no]
-o, --output_dir                          if omitted, generates a results directory
                                                            in the current position                         
-n, --threads_number                               number of processing threads [1]
-v, --verbose           increase verbosity: 0 = only warnings, 1 = info, 2 = debug,
                                        3 = debug on terminal. No number means info
                                                            Default is no verbosity
-h, --help                                          show this help message and exit
    ''')
    f.add_argument('-a', '--aligner', action="store", nargs='*', default=['m'], metavar='',
                   help=argparse.SUPPRESS, choices=['m', 'b', 'n'])
    f.add_argument('-s', '--alignment_stats', action="store", nargs='*', default=['n'], metavar='',
                   help=argparse.SUPPRESS, choices=['y', 'n', 'yes', 'no'])
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
        if not args.inputs or not args.reference or not args.prefix:
            parser.print_help()
            log.error('Missing mandatory arguments!')
            sys.exit(1)
        else:
            MainReader(args)
    except IOError as msg:
        parser.error(str(msg))
