#!/usr/bin/env python
import os
import sys
import time
import argparse
import glob
import subprocess
from multiprocessing import Pool

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
try:
    from Bio import SeqIO
except ImportError:
    sys.stderr.write("BioPython SeqIO is necessary to run postfiltering.\n")
    sys.exit(1)


basedir = os.path.dirname(os.path.abspath(__file__))

############################################################
### Parse Commandline Arguments

def get_arguments():
    '''
        Get commandline arguments and return namespace
        '''
    ## Initialize Parser
    parser = argparse.ArgumentParser(prog='metaSNV_pnps.py', description='metaSNV pnps ration computation', epilog='''Note:''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Not Showns:
    parser.add_argument('--version', action='version', version='%(prog)s 2.0', help=argparse.SUPPRESS)
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)
    
    # REQUIRED  arguments:
    parser.add_argument('--filt', metavar=': Filtered frequency files', help="Folder containing /pop/*.filtered.freq", required = True)
    parser.add_argument('--fasta', metavar=': fasta database', help="fasta database used for mapping", required = True)
    parser.add_argument('--annotations', metavar=': annotation file', help="annotation file used at the first steps", required = True)
    parser.add_argument('--nucl_table', metavar=': nucleotide table', help="nucleotide_table.csv, found in ~/DEV_METASNV/", required = True)
    
    return parser.parse_args()

############################################################
### Basic functions

def file_check():
    ''' Check if required files exist (True / False)'''
    args.projdir = '/'.join(args.filt.rstrip('/').split('/')[:-1])
    args.pnps = args.projdir + '/pnps/'
    args.temp = args.projdir + '/pnps/temp/'
    
    print("Checking for necessary input files...")
    if os.path.isfile(args.fasta) and os.path.isfile(args.annotations) and os.path.isfile(args.nucl_table)):
        print("found: '{}' \nfound:'{}' \nfound:'{}'".format(args.fasta, args.annotations, args.nucl_table))
    else:
        sys.exit("\nERROR: No such file '{}',\nERROR: No such file '{}',\nERROR: No such file '{}'".format(args.fasta, args.annotations, args.nucl_table))


def print_arguments():
    print("")
    print("Checking required arguments:")
    if args.filt:
        print("Filtered folder : {}".format(args.filt) )
    if args.fasta:
        print("Fasta : {}".format(args.fasta) )
    if args.annotations:
        print("Annotations : {}".format(args.annotations) )
    if args.nucl_table:
        print("Nucleotide table : {}".format(args.nucl_table) )
    print("")


############################################################
### Compute Pn/Ps

def get_seq(file, gene):
    """Return the sequence corresponding to the gene identifier as a string"""
    for seq_record in SeqIO.parse(file, "fasta"):
        seq_gene = str(seq_record.id).split(':')[1]
        if seq_gene == gene:
            return str(seq_record.seq)


def get_codon_usage(string,codons):
    """Return the number of each codon from a string"""
    if len(string) % 3 == 0:
        split = [string[i:i+3] for i in range(0, len(string), 3)]
        return {x:split.count(x) for x in codons}
    else:
        return {x:np.NaN for x in codons}


def get_exp_ratio(file,codons):
    """Getting the expected ratio"""
    exp_S=0
    exp_N=0
    for seq_record in SeqIO.parse(file, "fasta"):
        seq = str(seq_record.seq)
        if len(seq) % 3 == 0:
            split = [seq[i:i+3] for i in range(0, len(seq), 3)]
            dict = {x:split.count(x) for x in codons}
        else:
            dict = {x:np.NaN for x in codons}
        temp_df = pd.DataFrame.from_dict(dict,orient="index")
        # Edit values
        exp_S += np.nansum(temp_df[0]*sig_table["Synonymous"])
        exp_N += np.nansum(temp_df[0]*sig_table["Non_Synonymous"])
    return exp_N/exp_S


def get_obs_ratio(data,i):
    obs_N = np.nansum(data.xs('N', level='significance')[i])
    obs_S = np.nansum(data.xs('S', level='significance')[i])
    return obs_N/float(obs_S)


def pnps(f):
    species = int(f.split('/')[-1].split('.')[0])
    print species
    
    samples = [line.split('/')[-1].replace('\n','') for line in f]
    species_df = pd.DataFrame(index=species, columns=[samples])
    
    cmd1 = ['cut -f1', f ,' | cut -f2 -d\':\' | grep -v - | grep -v ^$ | sort -u > ', args.temp + species + '.genes']
    subprocess.call(cmd1)
    cmd2 = ['grep -f ', args.temp + species + '.genes', args.annotations, ' | awk \'BEGIN { FS = OFS = "\t" } { print $3,$7-1,$8,$3 ":" $2,0, $9}\' > ', args.temp + species + '.bed']
    subprocess.call(cmd2)
    cmd3 = ['bedtools getfasta -fi ', args.fasta, ' -bed ', args.temp + species + '.bed', ' -s -name > ', args.temp + species + '.fa']
    subprocess.call(cmd3)
    
    fastafile = args.temp + species + '.fa'
    
    data = pd.read_table(f, index_col=0, na_values=['-1'])
    pre_index = [i.split(':') for i in list(data.index)]
    # Setting index for each genes
    index1 = [item[1] for item in pre_index]
    # Setting index for Non-synonymous vs Synonymous
    index2 = [item[4].split('[')[0] for item in pre_index]
    index = pd.MultiIndex.from_arrays([index1, index2], names=['genes', 'significance'])
    data = data.set_index(index)
    data = data.sort_index()
    genes = list(data.index.levels[0])
    genes = [x for x in genes if x != '-']
    data = data.loc[genes]
    
    exp = get_exp_ratio(fastafile,codons)
    for i in data.columns:
        species_df.loc[sp][i] = get_obs_ratio(data,i)/exp
    
    species_df.to_csv(projdir + '/pnps/' + 'summary.pnps', sep='\t')

    return all_df


def Allpnps(args):
    
    sig_table = pd.read_table(args.nucl_table, index_col=0)
    codons = ["GCA","GCC","GCG","GCT","TGC","TGT","GAC","GAT","GAA","GAG","TTC","TTT","GGA","GGC","GGG","GGT","CAC","CAT","ATA","ATC","ATT","AAA","AAG","CTA","CTC","CTG","CTT","TTA","TTG","ATG","AAC","AAT","CCA","CCC","CCG","CCT","CAA","CAG","AGA","AGG","CGA","CGC","CGG","CGT","AGC","AGT","TCA","TCC","TCG","TCT","TAA","TAG","TGA","ACA","ACC","ACG","ACT","GTA","GTC","GTG","GTT","TGG","TAC","TAT"]
    
    allFreq = glob.glob(args.filt + '/pop/*.freq')
    allSpeceis =

    for f in allFreq:
        pnps(f)






    
    








