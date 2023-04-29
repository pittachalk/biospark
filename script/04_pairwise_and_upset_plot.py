#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: alymakhlouf
"""

# %%
from pathlib import Path
import pandas as pd
import numpy as np
import itertools
import tqdm
from upsetplot import from_memberships, plot
from matplotlib import pyplot as plt

directory = Path(__file__).resolve().parent.parent.joinpath("data/")
upset_directory = Path(__file__).resolve().parent.joinpath("figure/05_upset_plots/")

df = pd.read_csv(directory.joinpath('filtered_hit_list_summarised.csv'))

# %%
# store broad cell type-specific hit genes in a dictionary

broad_cell_types, counts_cell_types = np.unique(df['cell_type'].values, return_counts = True)
hit_genes_dict = {}

for i, cell_type in enumerate(broad_cell_types[0:]):
    
    hit_genes_dict[cell_type] = np.array(df.loc[df['cell_type'] == cell_type]['gene'])
    
# %%
# find intersections in hit gene sets between all pairs of broad cell types

pairs = list(itertools.permutations(broad_cell_types,2))

pairs_upset = [] # make pairs list format compatible with upsetplot library
counts_upset = [] # store count values for upset plot

# initialise dataframe to be written
df_intersect = pd.DataFrame(data=None, columns=["cell_type1", "cell_type2", "n_intersect", "n_genes1", "intersect_list"]) 

for j, pair in enumerate(pairs):

    pairs_upset.append(list(pair))
    
    genes_1 = np.array(df.loc[df['cell_type'] == pair[0]]['gene']) # gene list for broad cell type 1
    genes_2 = np.array(df.loc[df['cell_type'] == pair[1]]['gene']) # gene list for broad cell type 2
    
    # find intersection size of hit genes between the two broad cell types
    genes_intersect = np.intersect1d(genes_1, genes_2)
    counts_upset.append(len(genes_intersect))

    # concatenate to dataframe
    tmp_df = pd.DataFrame({'cell_type1': [pair[0]], 'cell_type2': [pair[1]],
                           'n_intersect': [len(genes_intersect)], 
                           'n_genes1': [len(genes_1)], 
                           'intersect_list': ','.join(genes_intersect)})
    df_intersect = pd.concat([df_intersect, tmp_df], axis = 0)

# %%    
# take subsets of the lists of all pairs and all intersections matching each broad cell type

count = 0
subset = len(broad_cell_types) - 1 # define the size of the subset

for i, cell_type in enumerate(tqdm.tqdm(broad_cell_types)):  
    
    upset = from_memberships(pairs_upset[count : count + subset], 
            (counts_upset[count : count + subset])/counts_cell_types[i]) # normalise intersection size by the number of hit genes for the given broad cell type
    
    # generate and save UpSet plots for the given broad cell type
    plot(upset) 
    
    file_name = cell_type + '.png'
    file_directory = upset_directory.joinpath(file_name)
    
    plt.savefig(file_directory) 
    plt.close()
    
    count += subset # update subset to move on to the next broad cell type


# %%

# generate csv file of pairs and intersection size

# matching to self
df_self = pd.DataFrame({'cell_type1': hit_genes_dict.keys(), 
                        'cell_type2': hit_genes_dict.keys(), 
                        'n_intersect': [len(x) for x in hit_genes_dict.values()],
                        'n_genes1': [len(x) for x in hit_genes_dict.values()],
                        'intersect_list': [','.join(x) for x in hit_genes_dict.values()] })

# concatenate and sort
df_out = (pd.concat([df_intersect, df_self], axis=0)
    .sort_values(by=["cell_type1", "cell_type2"])
)

df_out.to_csv(directory.joinpath("pairwise_intersects.txt"), sep='\t', index=False)

# %%
