#!/usr/bin/env python
import os
import sys
import time
import argparse
import glob
from shutil import copyfile
from multiprocessing import Pool
from functools import partial

try:
    import numpy as np
except ImportError:
    sys.stderr.write("Numpy is necessary to run postfiltering.\n")
    sys.exit(1)
try:
    import pandas as pd
except ImportError:
    sys.stderr.write("Pandas is necessary to run postfiltering.\n")
    sys.exit(1)


basedir = os.path.dirname(os.path.abspath(__file__))


############################################################
### Parse Commandline Arguments

def get_arguments():
    '''
    Get commandline arguments and return namespace
    '''
    ## Initialize Parser
    parser = argparse.ArgumentParser(prog='metaSNV_filtomics.py', description='metaSNV joind metaT and metaG', epilog='''Note:''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Not Showns:
    parser.add_argument('--version', action='version', version='%(prog)s 2.0', help=argparse.SUPPRESS)
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

    # REQUIRED  arguments:
    parser.add_argument('--filtered', metavar='filtered freq', help="Folder containing the filtered /pop/*.filtered.freq to filter", required = True)
    
    # OPTIONAL  arguments:
    parser.add_argument('-b', metavar='FLOAT', type=float, default=40.0, help="Coverage breadth: minimal horizontal genome coverage percentage per sample per species")
    parser.add_argument('-p',metavar='FLOAT', type=float, help="FILTERING STEP II: required proportion of informative samples (coverage non-zero) per position", default=0.50)
    parser.add_argument('--n_threads',metavar=': Number of Processes',default=1,type=int, help="Number of jobs to run simmultaneously.")

    return parser.parse_args()

############################################################
### Basic functions

def file_check():
    ''' Check if required files exist (True / False)'''
    args.projdir = '/'.join(args.filtered.rstrip('/').split('/')[:-1])
    
    args.coverage_file = args.projdir+'/'+args.projdir.split('/')[-1]+'.all_cov.tab'
    args.percentage_file = args.projdir+'/'+args.projdir.split('/')[-1]+'.all_perc.tab'

    print("Checking for necessary input files...")
    if os.path.isfile(args.coverage_file) and os.path.isfile(args.percentage_file):
        print("found: '{}' \nfound:'{}'".format(args.coverage_file, args.percentage_file))
    else:
        sys.exit("\nERROR: No such file '{}',\nERROR: No such file '{}'".format(args.coverage_file, args.percentage_file))


def print_arguments():
    ## Print Defined Arguments:
    print("")
    print("Checking parameters:")
    if args.filtered:
        print("Filtered folder {}".format(args.filtered) )
    if args.b:
        print("MetaG threshold: percentage covered (breadth) {}".format(args.b) )
    if args.p:
        print("MetaG threshold: Min. proportion of covered samples in samples_of_interest {}".format(args.p) )
    if args.n_threads:
        print("Number of parallel processes : {}".format(args.n_threads) )
    print("")


def filter_metaG(filt_file, horizontal_coverage, outdir):
    
    species = int(filt_file.split('/')[-1].split('.')[0])
    
    data = pd.read_table(filt_file, index_col=0, na_values=['-1'])
    data = data.sort_index()
    
    samples_G = [sample for sample in data.columns if '_G.filtered' in sample]
    samples_drop = list(horizontal_coverage.loc[species,samples_G][horizontal_coverage.loc[species,samples_G] < args.b].index)
    data = data.drop(samples_drop, axis=1)
    
    def filt_proportion(x):
        n = np.count_nonzero(np.isnan(x))
        return n>(len(x)*(1-args.p))
    
    samples_G = [sample for sample in data.columns if '_G.filtered' in sample]
    index_drop = [index for index in data.index if filt_proportion(data.loc[index, samples_G])]
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

    outdir = args.projdir + '/' + args.filtered.rstrip('/').split('/')[-1] + '_filtomics'

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(outdir + '/pop/'):
        os.makedirs(outdir + '/pop/')

    allFiles = glob.glob(args.filtered + '/pop/*.freq')

    print "Filtering..."
    
    # Load external info : Coverage, genomes size, genes size
    horizontal_coverage = pd.read_table(args.percentage_file, skiprows=[1], index_col=0)

    p = Pool(processes = args.n_threads)
    partial_filt = partial(filter_metaG,  horizontal_coverage = horizontal_coverage, outdir = outdir)
    p.map(partial_filt, allFiles)
    p.close()
    p.join()












