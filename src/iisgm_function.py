#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Required libraries

import pandas as pd
import os
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, ClusterWarning, set_link_color_palette
from warnings import simplefilter

simplefilter("ignore", ClusterWarning)
set_link_color_palette(['C0'])


# In[ ]:


# ITERATION 1

def read_vcf(file):
    
    """ Reads .vcf file (its path should be specified as input).
        Returns variant dataframe. """
    
    with open(file, 'r') as f:
        
        for line in f:                                               # Get header line
            if line.startswith('#CHROM'): hd = line; break
                
        hd = hd[1:].strip().split('\t')                              # Remove # from header and convert to list
        df = pd.read_csv(f, delimiter='\t', names=hd, comment='#')   # Save vcf as dataframe ignoring comments
        
    return df


# In[ ]:


# ITERATION 2

def filter_vcf(df):
    
    """ Removes low quality variants and heterozygous genotypes. 
        Returns filtered dataframe with extra columns indicating the haploid genotype, 
        the allele in the position and the nature of the variant (SNP or INDEL). """
    
    f_df = df.copy()                                                 # Copy dataframe to be filtered
    f_df.drop(f_df[f_df.FILTER != 'PASS'].index, inplace=True)       # Drop low quality sites
    
    for inf in f_df.itertuples():                                    # Drop "heterozygous"
        l_line = inf[-1]
        if l_line[0] != l_line[2]: f_df.drop(inf[0], inplace=True)   # If the alleles are not equal, drop
    
    f_df['GT'] = f_df.apply(lambda x: x[-1][0], axis=1)              # Get haploid genotype
    f_df['VAR'] = f_df.apply(lambda x: x.ALT.split(',')[int(x.GT)-1], axis=1)
    f_df['TYPE'] = ['SNP' if (len(inf.REF) == len(inf.VAR)) and (inf.VAR != '*')
                    else 'INDEL' for inf in f_df.itertuples()]       # Classify variants as SNP or INDEL
    
    return f_df


# In[ ]:


# ITERATION 3

def presence_matrix(direct, indels=False):
    
    """ Reads all .vcf files from specified path.
        Turns mutation dict to presence/absence matrix with mutations (POS_REF_VAR)
        as columns and organism IDs as rows. """
    
    dir_l = [direct + f for f in os.listdir(direct) if f.endswith('.vcf')] # Get vcf files in directory
    
    d = defaultdict(dict)                                                  # Create empty mutation dict
    for file in dir_l:                                                     # Read and filter mutations
        df = read_vcf(file)
        f_df = filter_vcf(df)
        d = d_mut(f_df, d, indels)                                         # Add items to mutation dict
        
    pr_mat = pd.DataFrame(d).fillna(0).astype(int)                         # Mutation dict to 1/0 matrix
    
    return pr_mat

def d_mut(df, d, indels):
    
    """ From variant dataframe, updates double key mutation dict --> d[mutation][organism] = 1. """
    
    org = df.columns[-4]                                                          # Get sample id
    if not indels: df.drop(df[df.TYPE == 'INDEL'].index, inplace=True)            # Drop INDELS if required
    mut = df[['POS', 'REF', 'VAR']].apply(lambda x: 
                                        '_'.join(x.astype(str)), axis=1).tolist() # List of mutation ids (POS_REF_VAR) in organism
    for m in mut: d[m][org] = 1                                                   # Add items to dict
    
    return d


# In[ ]:


# ITERATIONS 4, 5, 6

def mut_dist_dnd(pr_mat, metric='jaccard', method='average', path=False):
    
    """ Calculates pairwise distance between samples 
        and plots dendrogram of the results. """
    
    linkage_matrix = linkage(pr_mat, metric=metric, method=method) # Pairwise distance according to specified metric
    
    plt.figure(figsize=(10, 7))
    plt.grid(axis='x')
    dnd = dendrogram(linkage_matrix,                               # Plot dendrogram
               labels = pr_mat.index,
               orientation='left')
    plt.title('Distance dendrogram - ' + metric, fontsize=15)
    
    if path : plt.savefig(path, bbox_inches='tight')               # Save plot if specified
   
    plt.show()
    
    return linkage_matrix
