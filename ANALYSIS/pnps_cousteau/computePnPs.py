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

def get_codon_usage(string,codons):
    """Return the number of each codon from a string"""
    if len(string) % 3 == 0:
        split = [string[i:i+3] for i in range(0, len(string), 3)]
        return {x:split.count(x) for x in codons}
    else:
        return {x:np.NaN for x in codons}


def get_seq(file, gene):
    """Return the sequence corresponding to the gene identifier as a string"""
    for seq_record in SeqIO.parse(file, "fasta"):
        seq_gene = str(seq_record.id).split(':')[1]
        if seq_gene == gene:
            return str(seq_record.seq)


def get_exp(file,codons,g):
    # Retrieve sequence
    string = get_seq(file,g)
    dict = get_codon_usage(string, codons)
    temp_df = pd.DataFrame.from_dict(dict,orient="index")
    if (sum(temp_df[0]*sig_table["Synonymous"]) == 0) or (sum(temp_df[0]*sig_table["Non_Synonymous"]) == 0):
        return  np.NaN
    else:
        return sum(temp_df[0]*sig_table["Non_Synonymous"])/sum(temp_df[0]*sig_table["Synonymous"])


def compute_obs_ratio(data,g,i):
    if 'S' not in data.loc[g][i].index:
        return np.NaN
    elif np.nansum(data.loc[g,"S"][i]) == 0:
        return np.NaN
    else:
        if 'N' not in data.loc[g][i].index:
            return np.NaN
        else:
            return np.nansum(data.loc[g,"N"][i])/np.nansum(data.loc[g,"S"][i])


###############
# ARGS & VARS #
###############

inputfile = sys.argv[1]

print inputfile

species = inputfile.split('/')[-1].replace('.filtered.freq', '')

projdir = inputfile.replace('filtered/pop/'+species+'.filtered.freq', '')

codons = ["GCA","GCC","GCG","GCT","TGC","TGT","GAC","GAT","GAA","GAG","TTC","TTT","GGA","GGC","GGG","GGT","CAC","CAT","ATA","ATC","ATT","AAA","AAG","CTA","CTC","CTG","CTT","TTA","TTG","ATG","AAC","AAT","CCA","CCC","CCG","CCT","CAA","CAG","AGA","AGG","CGA","CGC","CGG","CGT","AGC","AGT","TCA","TCC","TCG","TCT","TAA","TAG","TGA","ACA","ACC","ACG","ACT","GTA","GTC","GTG","GTT","TGG","TAC","TAT"]

sig_table = pd.read_table(projdir + 'nucleotide_table.csv', index_col=0)
fastafile = projdir + 'pnps/fasta/' + species + '.fa'

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


##########
# SCRIPT #
##########

species_df = pd.DataFrame(index=genes, columns=[data.columns])

for g in genes:
    print('*******')
    print(g)
    print('*******')
    obs = [compute_obs_ratio(data,g,i) for i in data.columns]
    exp = get_exp(fastafile,codons,g)
    species_df.loc[g] = obs/exp

species_df.to_csv(projdir + '/pnps/' + '%s.pnps' % species, sep='\t')




