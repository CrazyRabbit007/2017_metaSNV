#!/usr/bin/env python

# This code is part of the metagenomic SNP calling pipeline (metaSNP)
# Copyright (c) 2016 Robin Munch
# Licenced under the GNU General Public License (see LICENSE) 


import os
import sys
import argparse
import numpy as np
import pandas as pd

#############################
# Parse Commandline Arguments
#############################
def get_arguments():
    '''
    Get commandline arguments and return namespace
    '''
    ## Initialize Parser
    parser = argparse.ArgumentParser(
			prog='metaSNP_distances.py', 
			description='Compute pairwise distances between samples from species frequency tables.',
			epilog='''Note:''', 
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# REQUIRED  arguments:
    parser.add_argument("in_freq",  type=str, help="Input species frequency table (metaSNP_filtering.py)") #, metavar='FILE')
    parser.add_argument("outdir",  type=str, help="Write output to outdir") #, metavar='FILE')

    return parser.parse_args()

    ##############################

def l1nonans(d1,d2):
    return np.abs(d1 - d2).mean()

def alleledist(d1,d2, threshold=.6):
    return (np.abs(d1 - d2) > threshold).mean()


if __name__ == "__main__":
    args = get_arguments()

f = args.in_freq
sampleName=os.path.basename(args.in_freq)

data = pd.read_table(f, index_col=0, na_values=['-1']).T
dist = [[l1nonans(data.iloc[i], data.iloc[j]) for i in range(len(data))]
            for j in range(len(data))]
dist = pd.DataFrame(dist, index=data.index, columns=data.index)

dist.to_csv(args.outdir+'/'+'%s.mann.dist' % sampleName, sep='\t')

dist = [[alleledist(data.iloc[i], data.iloc[j]) for i in range(len(data))]
            for j in range(len(data))]
dist = pd.DataFrame(dist, index=data.index, columns=data.index)

dist.to_csv(args.outdir+'/'+'%s.allele.dist' % sampleName, sep='\t')

