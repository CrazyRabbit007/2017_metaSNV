#!/usr/bin/env python
import os
import sys
import time
import argparse
import glob
import subprocess
import shutils
from collections import Counter
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

def get_exp_ratio(file,genes):
    """Getting the expected ratios from a species fasta file"""
    
    genes_exp = pd.DataFrame(index = genes, columns = ['Exp_ratio'])
    all_dict=Counter({})
    
    for seq_record in SeqIO.parse(file, "fasta"):
        seq = str(seq_record.seq)
        seq_gene = str(seq_record.id).split(':')[1]
        if len(seq) % 3 == 0:
            split = [seq[i:i+3] for i in range(0, len(seq), 3)]
            gene_dict = Counter({x:split.count(x) for x in codons})
        else:
            gene_dict = Counter({x:np.NaN for x in codons})
        
        gene_df = pd.DataFrame.from_dict(gene_dict,orient="index")
        if (np.nansum(gene_df[0]*sig_table["Synonymous"]) == 0) or (np.nansum(gene_df[0]*sig_table["Non_Synonymous"]) == 0):
            gene_ratio = np.NaN
        else:
            gene_ratio = np.nansum(gene_df[0]*sig_table["Non_Synonymous"])/np.nansum(gene_df[0]*sig_table["Synonymous"])
        genes_exp.iloc[seq_gene][0] = gene_ratio
        all_dict += gene_dict

    all_df = pd.DataFrame.from_dict(all_dict,orient="index")
    if (np.nansum(all_df[0]*sig_table["Synonymous"]) == 0) or (np.nansum(all_df[0]*sig_table["Non_Synonymous"]) == 0):
        all_ratio = np.NaN
    else:
        all_ratio = np.nansum(all_df[0]*sig_table["Non_Synonymous"])/np.nansum(all_df[0]*sig_table["Synonymous"])

    return all_ratio, genes_exp

def compute_gene_obs_ratio(data,i,genes):
    genes_obs = []
    for g in genes :
        if 'S' not in data.loc[g][i].index:
            genes_obs.append(np.NaN)
        elif np.nansum(data.loc[g,"S"][i]) == 0:
            genes_obs.append(np.NaN)
        else:
            if 'N' not in data.loc[g][i].index:
                genes_obs.append(np.NaN)
            else:
                genes_obs.append(np.nansum(data.loc[g,"N"][i])/np.nansum(data.loc[g,"S"][i]))
    return genes_obs

def compute_all_obs_ratio(data,i):
    obs_N = np.nansum(data.xs('N', level='significance')[i])
    obs_S = np.nansum(data.xs('S', level='significance')[i])
    return obs_N/float(obs_S)


def pnps(f, all_df):
    species = int(f.split('/')[-1].split('.')[0])
    print species
    
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
    data = data.set_index(pd.MultiIndex.from_arrays([index1, index2], names=['genes', 'significance']))
    data = data.sort_index()
    genes = list(data.index.levels[0])
    genes = [x for x in genes if x != '-']
    data = data.loc[genes]
    
    species_df = pd.DataFrame(index=genes, columns=data.columns)
    
    all_exp_ratio, genes_exp_df = get_exp_ratio(fastafile,genes)
    
    for i in data.columns:
        all_df.loc[sp][i] = compute_all_obs_ratio(data,i)/all_exp_ratio
        species_df.loc[:,i] = compute_genes_obs_ratio(data,i)/genes_exp_df.loc[:,0]
    
    species_df.to_csv(args.pnps + '/' + species + '.pnps', sep='\t')

    return all_df


def Allpnps(args):
    
    sig_table = pd.read_table(args.nucl_table, index_col=0)
    codons = ["GCA","GCC","GCG","GCT","TGC","TGT","GAC","GAT","GAA","GAG","TTC","TTT","GGA","GGC","GGG","GGT","CAC","CAT","ATA","ATC","ATT","AAA","AAG","CTA","CTC","CTG","CTT","TTA","TTG","ATG","AAC","AAT","CCA","CCC","CCG","CCT","CAA","CAG","AGA","AGG","CGA","CGC","CGG","CGT","AGC","AGT","TCA","TCC","TCG","TCT","TAA","TAG","TGA","ACA","ACC","ACG","ACT","GTA","GTC","GTG","GTT","TGG","TAC","TAT"]
    
    allFreq = glob.glob(args.filt + '/pop/*.freq')
    allSpecies = [sp.split('/')[-1].strip('.filtered.freq') for sp in allFreq]
    with open(args.projdir + 'all_samples') as samples:
        allSamples = samples.read().splitlines()
    allSamples = [samples.split('/')[-1] for samples in allSamples]

    all_df = pd.DataFrame(index = allSpecies, columns = allSamples)

    for f in allFreq:
        all_df = pnps(f, all_df)

    all_df.to_csv(args.pnps + '/Summary.pnps', sep='\t')


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
    
    if not os.path.exists(args.pnps):
        os.makedirs(args.pnps)
    if not os.path.exists(args.temp):
        os.makedirs(args.temp)
    
    Allpnps(args)

    shutil.rmtree(args.temp)




    
    








