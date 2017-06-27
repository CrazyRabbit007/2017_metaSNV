#!/usr/bin/env python

##########
# IMPORT #
##########

import os       # interacting with operating system
import sys		# interact with system commandline 
import getopt   # manage command line options and arguments

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


#############
# FUNCTIONS #
#############

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
    obs_N = np.nansum([base > 0 for base in data.xs('N', level='significance')[i]])
    obs_S = np.nansum([base > 0 for base in data.xs('S', level='significance')[i]])
    return obs_N/float(obs_S)


###############
# ARGS & VARS #
###############

inputdir = sys.argv[1]

print inputdir

species = [sp.replace('.filtered.freq','') for sp in os.listdir(inputdir+'/')]

projdir = inputdir.replace('filtered/pop', '')

codons = ["GCA","GCC","GCG","GCT","TGC","TGT","GAC","GAT","GAA","GAG","TTC","TTT","GGA","GGC","GGG","GGT","CAC","CAT","ATA","ATC","ATT","AAA","AAG","CTA","CTC","CTG","CTT","TTA","TTG","ATG","AAC","AAT","CCA","CCC","CCG","CCT","CAA","CAG","AGA","AGG","CGA","CGC","CGG","CGT","AGC","AGT","TCA","TCC","TCG","TCT","TAA","TAG","TGA","ACA","ACC","ACG","ACT","GTA","GTC","GTG","GTT","TGG","TAC","TAT"]

sig_table = pd.read_table(projdir + 'nucleotide_table.csv', index_col=0)

##########
# SCRIPT #
##########
sample_file=str(sys.argv[2])
file = open(projdir+sample_file, 'r')
samples = [line.split('/')[-1].replace('\n','') for line in file]

species_df = pd.DataFrame(index=species, columns=[samples])

for sp in species:
    
    print('*******')
    print(sp)
    print('*******')
    inputfile = inputdir+sp+'.filtered.freq'
    fastafile = projdir + 'pnps/fasta/' + sp + '.fa'
    data = pd.read_table(inputfile, index_col=0, na_values=['-1'])
    pre_index = [i.split(':') for i in list(data.index)]
    # Setting index for each position
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


species_df.to_csv(projdir + '/pnps/' + 'schloiss_summary.pnps', sep='\t')




