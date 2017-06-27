#!/usr/bin/env python
import os
import sys
import time
import argparse
import glob
from shutil import copyfile

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
    parser = argparse.ArgumentParser(prog='metaSNV_post.py', description='metaSNV joind metaT and metaG', epilog='''Note:''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Not Showns:
    parser.add_argument('--version', action='version', version='%(prog)s 2.0', help=argparse.SUPPRESS)
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

    # REQUIRED  arguments:
    parser.add_argument('--filtG', metavar='metaG filtered freq', help="Folder containing the metaG /pop/*.filtered.freq to merge", required = True)
    parser.add_argument('--filtT', metavar='metaT filtered freq', help="Folder containing the metaT /pop/*.filtered.freq to merge", required = True)
    parser.add_argument('--out',metavar='Output directory', help="Output directory, where 'tara_join_metaSNV' will be", required = True)
    
    # OPTIONAL  arguments:
    parser.add_argument('--match',action='store_true', help="Subselect expressed genes")

    return parser.parse_args()

############################################################
### Basic functions

def file_check():
    ''' Check if required files exist (True / False)'''
    args.projdir_G = '/'.join(args.filtG.rstrip('/').split('/')[:-1])
    args.projdir_T = '/'.join(args.filtT.rstrip('/').split('/')[:-1])
    
    args.pars_G = args.filtG.rstrip('/').split('/')[-1].strip('filtered')
    args.pars_T = args.filtT.rstrip('/').split('/')[-1].strip('filtered')
    
    args.coverage_file_G = args.projdir_G+'/'+args.projdir_G.split('/')[-1]+'.all_cov.tab'
    args.percentage_file_G = args.projdir_G+'/'+args.projdir_G.split('/')[-1]+'.all_perc.tab'
    args.all_samples_G = args.projdir_G+'/'+'all_samples'
    
    args.coverage_file_T = args.projdir_T+'/'+args.projdir_T.split('/')[-1]+'.all_cov.tab'
    args.percentage_file_T = args.projdir_T+'/'+args.projdir_T.split('/')[-1]+'.all_perc.tab'
    args.all_samples_T = args.projdir_T+'/'+'all_samples'

    print("Checking for necessary input files...")
    print("\nMeta G...")
    if os.path.isfile(args.coverage_file_G) and os.path.isfile(args.percentage_file_G) and os.path.isfile(args.all_samples_G):
        print("found: '{}' \nfound:'{}' \nfound:'{}'".format(args.coverage_file_G, args.percentage_file_G, args.all_samples_G))
    else:
        sys.exit("\nERROR: No such file '{}',\nERROR: No such file '{}',\nERROR: No such file '{}'".format(args.coverage_file_G, args.percentage_file_G, args.all_samples_G))
    print("\nMeta T...")
    if os.path.isfile(args.coverage_file_T) and os.path.isfile(args.percentage_file_T) and os.path.isfile(args.all_samples_T):
        print("found: '{}' \nfound:'{}' \nfound:'{}'".format(args.coverage_file_T, args.percentage_file_T, args.all_samples_T))
    else:
        sys.exit("\nERROR: No such file '{}',\nERROR: No such file '{}',\nERROR: No such file '{}'".format(args.coverage_file_T, args.percentage_file_T, args.all_samples_T))

    print("\nChecking for output directory...")
    if os.path.exists(args.out):
        print("found: '{}'".format(args.out))
    else:
        sys.exit("\nERROR: No such directory '{}'".format(args.out))


def print_arguments():
    ## Print Defined Arguments:
    print("")
    print("Checking parameters:")
    if args.filtG:
        print("Filtered folder, metaG : {}".format(args.filtG) )
    if args.filtT:
        print("Filtered folder, metaT : {}".format(args.filtT) )
    if args.out:
        print("Output directory : {}".format(args.out) )
    if args.match:
        print("Matching expressed genes : {}".format(args.match) )
    print("")

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

    outdir = args.out + '/tara_joined_metaSNV'

    if args.match:
        outdir_filt = outdir + '/filtered.G' + args.pars_G + '.T' + args.pars_T + '.matched_genes'
    else:
        outdir_filt = outdir + '/filtered.G' + args.pars_G + '.T' + args.pars_T

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(outdir_filt):
        os.makedirs(outdir_filt)
    if not os.path.exists(outdir_filt + '/pop/'):
        os.makedirs(outdir_filt + '/pop/')

    allFreq_G = glob.glob(args.filtG + '/pop/*.freq')
    species_G = [f.split('/')[-1].split('.')[0] for f in allFreq_G]

    allFreq_T = glob.glob(args.filtT + '/pop/*.freq')
    species_T = [f.split('/')[-1].split('.')[0] for f in allFreq_T]

    species = list(set(species_G).intersection(species_T))

    print('\nWe find : {} species. '.format(len(species)))
    print(species)

    coverage_joined = pd.concat([pd.read_table(args.coverage_file_G, index_col=0),pd.read_table(args.coverage_file_T, index_col=0)], axis = 1)
    coverage_joined.to_csv(outdir + '/' + outdir.split('/')[-1] + '.all_cov.tab', sep='\t')
    percentage_joined = pd.concat([pd.read_table(args.percentage_file_G, index_col=0),pd.read_table(args.percentage_file_T, index_col=0)], axis = 1)
    percentage_joined.to_csv(outdir + '/' + outdir.split('/')[-1] + '.all_perc.tab', sep='\t')
    copyfile(args.projdir_G+'/bed_header',outdir+'/bed_header')

    for sp in species:
        
        dataG = pd.read_table(args.filtG + '/pop/' + sp + '.filtered.freq', index_col=0, na_values=['-1'])
        index_G = [i.split(':') for i in list(dataG.index)]
        dataG = dataG.set_index(pd.Index([item[0] + ':' + item[1] + ':' + item[2] + ':' + item[3] for item in index_G]))
        dataG = dataG.sort_index()

        dataT = pd.read_table(args.filtT + '/pop/' + sp + '.filtered.freq', index_col=0, na_values=['-1'])
        index_T = [i.split(':') for i in list(dataT.index)]
        dataT = dataT.set_index(pd.Index([item[0] + ':' + item[1] + ':' + item[2] + ':' + item[3] for item in index_T]))
        dataT = dataT.sort_index()
        
        if args.match:
            genes_T = set([item[1] for item in index_T])
            genes_T = [x for x in genes_T if x != '-']
            subset_G = [item.split(':')[1] in genes_T for item in list(dataG.index)]
            dataG_sub = dataG[subset_G]
            subset_T = [item.split(':')[1] in genes_T for item in list(dataT.index)]
            dataT_sub = dataT[subset_T]
        else:
            dataG_sub = dataG
            dataT_sub = dataT

        
        data_joined = pd.concat([dataG_sub, dataT_sub], axis=1)
        
        data_joined.to_csv(outdir_filt + '/pop/' + '%s.filtered.freq' % sp, sep='\t')

        print(sp)
        print('Subsetted for annotated genes detected in metaT :')
        positionsG = set([item[0] + ':' + item[1] + ':' + item[2] for item in [i.split(':') for i in list(dataG_sub.index)]])
        print('Selected {} metaG position'.format(len(positionsG)))
        positionsT = set([item[0] + ':' + item[1] + ':' + item[2] for item in [i.split(':') for i in list(dataT_sub.index)]])
        print('Selected {} metaT position'.format(len(positionsT)))
        positionsBoth = len(positionsG.intersection(positionsT))
        print('Among which {} are shared position'.format(positionsBoth))








