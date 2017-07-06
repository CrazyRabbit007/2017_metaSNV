#!/usr/bin/env python
import os
import sys
import time
import argparse
import glob
import subprocess
from shutil import copyfile
from multiprocessing import Pool
from functools import partial

try:
    import numpy as np
except ImportError:
    sys.stderr.write("Numpy is necessary to run this script.\n")
    sys.exit(1)
try:
    import pandas as pd
except ImportError:
    sys.stderr.write("Pandas is necessary to run this script.\n")
    sys.exit(1)


basedir = os.path.dirname(os.path.abspath(__file__))

############################################################
### Parse Commandline Arguments

def get_arguments():
    '''
    Get commandline arguments and return namespace
    '''
    ## Initialize Parser
    parser = argparse.ArgumentParser(prog='metaSNV_universal.py', description='metaSNV extraction of universal genes', epilog='''Note:''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Not Showns:
    parser.add_argument('--version', action='version', version='%(prog)s 2.0', help=argparse.SUPPRESS)
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)
    
    # REQUIRED  arguments:
    parser.add_argument('--filt', metavar=': Filtered frequency files', help="Folder containing /pop/*.filtered.freq", required = True)
    parser.add_argument('--universal', metavar=': Universal genes ID', help="Universal genes ID for the database used", required = True)
    
    # OPTIONAL  arguments:
    parser.add_argument('--n_threads',metavar=': Number of Processes',default=1,type=int, help="Number of jobs to run simmultaneously.")
    
    return parser.parse_args()

############################################################
### Basic functions
    
def file_check():
    ''' Check if required files exist (True / False)'''
    args.projdir = '/'.join(args.filt.rstrip('/').split('/')[:-1])
    args.pars = args.filt.rstrip('/').split('/')[-1].strip('filtered')
    
    print("Checking for necessary input files...")
    if os.path.isfile(args.universal):
        print("found: '{}'".format(args.universal))
    else:
        sys.exit("\nERROR: No such file '{}'".format(args.universal))


def print_arguments():
    print("")
    print("Checking required arguments:")
    if args.filt:
        print("Filtered folder : {}".format(args.filt))
    if args.universal:
        print("Universal ids : {}".format(args.universal))
    if args.n_threads:
        print("Number of parallel processes : {}".format(args.n_threads) )
    print("")


############################################################
### Compute Statistics


def filter_universal(filt_file, universal, outdir):
    
    species = int(filt_file.split('/')[-1].split('.')[0])
    
    data = pd.read_table(filt_file, index_col=0, na_values=['-1'])
    data = data.sort_index()
    
    def check_universal(x):
        y = x.split(':')[1]
        return y not in universal
    
    index_drop = [index for index in data.index if check_universal(index)]
    data = data.drop(index_drop)
    
    data.to_csv(outdir + '/pop/' + '%s.filtered.freq' % species, sep='\t')
    

############################################################
### Script

if __name__ == "__main__":
    
    #####################
    
    args = get_arguments()
    print_arguments()
    file_check()
    
    if args.debug:
        print_arguments()

    #####################

    outdir = args.projdir + '/' + args.filt.rstrip('/').split('/')[-1] + '_universal'

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        if not os.path.exists(outdir + '/pop/'):
            os.makedirs(outdir + '/pop/')

    allFiles = glob.glob(args.filt + '/pop/*.freq')
        
    print "Universal only..."
    
    # Load external info : universal genes ID
    with open(args.universal) as f:
        universal = f.read().splitlines()
    
    p = Pool(processes = args.n_threads)
    partial_filt = partial(filter_universal,  universal = universal, outdir = outdir)
    p.map(partial_filt, allFiles)
    p.close()
    p.join()












