#! /usr/bin/env python
import argparse
import sys

def run_subtool(mod,ags, ps):
    if mod == 'seqstats':
        import seqstats as sq 
    elif mod == 'fastqgen':
        import fastqgen as sq 
    elif mod == 'alignment':
        import alignment as sq 
    else:
        ps.print_help()
        print("pypore: error: invalid command: '%s' (choose from 'seqstats', 'fastqgen', 'alignment')" % (mod))
        sys.exit(1)
    sq.run(ags)

def main():
    parser = argparse.ArgumentParser(add_help=False,
                                     prog='PyPore',
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


usage: pypore [command] [options]''',
                                     epilog='''

PyPore. Written by Roberto Semeraro, Department of Clinical and Sperimental Medicine,
University of Florence. For bug report or suggestion write to robe.semeraro@gmail.com''',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     )
    g = parser.add_argument_group(title='Commands',
                                  description='''seqstats             analyze the entries in the directory given by user, generating
                     an interactive (HTML) experimental summary and, optionally, a FASTQ file.
fastqgen             convert FAST5 in FASTQ (faster than seqstats)
alignment            align fastq files and generate an interactive (HTML) alignment
                     report       
    ''')


    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit(1)
    else:
        mod=sys.argv[1]
        ags=[]
        if len(sys.argv) > 2:
        	ags=map(str, sys.argv[2:])
        run_subtool(mod, ags, parser)
if __name__ == "__main__":
    main()
