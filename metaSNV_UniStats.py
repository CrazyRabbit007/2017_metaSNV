#!/usr/bin/env python
import os
import sys
import time
import argparse
import glob

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
    parser = argparse.ArgumentParser(prog='metaSNV_stats.py', description='metaSNV descriptive statistics computation', epilog='''Note:''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Not Showns:
    parser.add_argument('--version', action='version', version='%(prog)s 2.0', help=argparse.SUPPRESS)
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)
    
    # REQUIRED  arguments:
    parser.add_argument('--filt', metavar=': Filtered frequency files', help="Folder containing /pop/*.filtered.freq", required = True)
    parser.add_argument('--bedfile', metavar=': Annotation file', help="Bedfile containing genomes annotations", required = True)
    
    # OPTIONAL  arguments:
    parser.add_argument('-b', metavar='FLOAT', type=float, default=40.0, help="Coverage breadth treshold used in filtering")
    parser.add_argument('-d', metavar='FLOAT', type=float, default=5.0, help="Coverage depth treshold used in filtering")
    
    return parser.parse_args()

############################################################
### Basic functions
    
def file_check():
    ''' Check if required files exist (True / False)'''
    args.projdir = '/'.join(args.filt.rstrip('/').split('/')[:-1])
    args.pars = args.filt.rstrip('/').split('/')[-1].strip('filtered')
    args.coverage_file = args.projdir+'/'+args.projdir.split('/')[-1]+'.all_cov.tab'
    args.percentage_file = args.projdir+'/'+args.projdir.split('/')[-1]+'.all_perc.tab'
    
    print("Checking for necessary input files...")
    if os.path.isfile(args.coverage_file) and os.path.isfile(args.percentage_file) and os.path.isfile(args.bedfile):
        print("found: '{}' \nfound:'{}' \nfound:'{}'".format(args.coverage_file, args.percentage_file, args.bedfile))
    else:
        sys.exit("\nERROR: No such file '{}',\nERROR: No such file '{}',\nERROR: No such file '{}'".format(args.coverage_file, args.percentage_file, args.bedfile))


def print_arguments():
    print("")
    print("Checking required arguments:")
    if args.filt:
        print("Filtered folder : {}".format(args.filt) )
    if args.bedfile:
        print("Filtered folder : {}".format(args.bedfile) )
    if args.pars:
        print("Suffixe : {}".format(args.pars) )
    print("Options:")
    if args.b:
        print("threshold: percentage covered (breadth) {}".format(args.b) )
    if args.d:
        print("threshold: average coverage (depth) {}".format(args.d) )
    print("")


############################################################
### Compute Statistics


def computeStats(args):

    def get_multiallelic(positions):
        seen = set()
        seen2 = set()
        seen3 = set()
        seen_add = seen.add
        seen2_add = seen2.add
        seen3_add = seen3.add
        for item in positions:
            if item in seen2:
                seen3_add(item)
            elif item in seen:
                seen2_add(item)
            else:
                seen_add(item)
        return list(seen), list(seen2), list(seen3)

    def get_stats(data, stats_df, species, type = ''):

        positions = list(data.index.get_level_values('position'))
        stats_df.loc[species,'Number of ' + type + 'SNVs'] = len(positions)

        multiallelics = get_multiallelic(positions)
        stats_df.loc[species,type+'Quadriallelic'] = len(multiallelics[2])
        stats_df.loc[species,type+'Triallelic'] = len(multiallelics[1]) - len(multiallelics[2])
        stats_df.loc[species,type+'Biallelic'] = len(multiallelics[0]) - len(multiallelics[1])
    
        return stats_df

    print "Computing Statistics"

    allFreq = glob.glob(args.filt + '/pop/*.freq')
    index = [int(f.split('/')[-1].split('.')[0]) for f in allFreq]
    columns = ['Number of variable positions','Number of SNVs','Number of S SNVs', 'Number of N SNVs',
    'Biallelic','Triallelic','Quadriallelic',
    'S Biallelic','S Triallelic','S Quadriallelic',
    'N Biallelic','N Triallelic','N Quadriallelic',
    'Length genome','Avg Hcov','Avg Vcov','SNVs per Kb']
    stats_df = pd.DataFrame(index=index, columns=columns)
    
    # Load external info : Coverage, genomes size, genes size
    horizontal_coverage = pd.read_table(args.percentage_file, skiprows=[1], index_col=0)
    vertical_coverage = pd.read_table(args.coverage_file, skiprows=[1], index_col=0)
    bedfile_tab = data = pd.read_table(args.bedfile, index_col=0, header =None,usecols=[1,5])
    
    for f in allFreq:
        species = int(f.split('/')[-1].split('.')[0])
        
        print species
        
        # Read the freq table as data.frame
        data = pd.read_table(f, index_col=0, na_values=['-1'])
        pre_index = [i.split(':') for i in list(data.index)]
        pos_index = [item[0] + ':' + item[1] + ':' + item[2] for item in pre_index]
        # Setting index for Non-synonymous vs Synonymous
        sig_index = [item[4] for item in pre_index]
        sig_index = [sig.split('[')[0] for sig in sig_index]
        data = data.set_index(pd.MultiIndex.from_arrays([pos_index, sig_index], names=['position', 'significance']))
        data = data.sort_index()
        
        list_genes = set([item[1] for item in pre_index])
        print list_genes
        genome_length = sum([int(x) for x in bedfile_tab.loc[list_genes,5].tolist()])

        stats_df.loc[species,'Length genome'] = genome_length
        stats_df.loc[species,'Avg Hcov'] = horizontal_coverage.loc[species, :][horizontal_coverage.loc[species, :]>=args.b].mean()
        stats_df.loc[species,'Avg Vcov'] = vertical_coverage.loc[species, :][vertical_coverage.loc[species, :]>=args.d].mean()
        
        positions = list(data.index.get_level_values('position'))
        stats_df.loc[species,'Number of variable positions'] = len(set(positions))
        
        stats_df = get_stats(data, stats_df, species)
        
        data_N = data.xs('N', level='significance')
        data_S = data.loc[(slice(None), ['.','S']),:]
        stats_df = get_stats(data_N, stats_df, species, 'N ')
        stats_df = get_stats(data_S, stats_df, species, 'S ')
        
        sp_dict = {Key : None for Key in data.columns}
        data_richness = data
        data_richness[data_richness == 1] = 0
        data_richness[data_richness > 0] = 1
        
        for col in data_richness.columns:
            positions_filtered = [i for (i, v) in zip(list(data_richness.index.get_level_values('position')), -np.isnan(data_richness.loc[:,col].values)) if v]
            sp_dict[col] = (np.nansum(data_richness.loc[:,col].values) + len(set(positions_filtered))) / len(set(positions_filtered))
        
        sp_richness = pd.DataFrame.from_dict(sp_dict, orient = 'index')
        sp_richness.to_csv(args.outdir + '/' + str(species) + args.pars + '.stats.tab', sep='\t')
    
    stats_df.loc[:,'SNVs per Kb'] = stats_df.loc[:,'Number of variable positions']/(stats_df.loc[:,'Length genome']/1000)
    stats_df.to_csv(args.projdir + '/' +args.projdir.split('/')[-1] + args.pars + '.stats.tab', sep='\t')
    

############################################################
### Script

if __name__ == "__main__":

    #####################

    args = get_arguments()
    file_check()
    print_arguments()
    
    if args.debug:
        print_arguments()
    
    #####################

    args.outdir = args.projdir + '/' + 'Stats' + args.pars
    
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    computeStats(args)











