#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %%
import pandas as pd
import numpy as np
import tqdm
import re
from pathlib import Path

directory = Path(__file__).resolve().parent.parent

# read the compiled hit list
# setting the first column containing the broad cell types as the index
df = pd.read_csv(directory.joinpath('data/filtered_hit_list.csv'), index_col=0)

# %%
# initialise empty master data frame
column_names_ordered = ['cell_type', 'gene', 'num_hits', 'num_screens']
df_master = pd.DataFrame(
    columns = column_names_ordered
)

# %%
# ITERATE THROUGH THE BROAD CELL TYPES
# COUNT THE NUMBER OF HITS FOR A GENE IN EACH CELL TYPE

# for each broad cell type...
for i in tqdm.tqdm(range(df.shape[0])):

    # get occurences of genes 
    current_hits = df.iloc[i, 1:]
    current_hits = current_hits.dropna()

    # retrieve cell type, remove spaces and substituting slashes
    current_cell = df.index[i]
    current_cell = re.sub("/", '-', current_cell)
    current_cell = re.sub(" ", '', current_cell)

    # retrieve the number of times each gene was a hit as a dataframe
    genes, counts = np.unique(current_hits, return_counts=True)
    df_subset = pd.DataFrame({'gene': genes, 'num_hits': counts})
    df_subset['num_hits'] = df_subset['num_hits'].astype("category")

    # add cell name column
    df_subset['cell_type'] = current_cell

    # add total number of screens done
    df_subset['num_screens'] = df.iloc[i, 0]

    # rearrange order of columns
    df_subset = df_subset[column_names_ordered]

    # concatenate to master dataframe
    df_master = pd.concat([df_master, df_subset])

# %%
df_master['num_screens'] = df_master['num_screens'].astype(int)

# %%
# add Ensembl ID to the genes
gene_id_mapping = (
    pd.read_table(directory.joinpath('data/homo_sapiens_gene_id_mapping.txt'))
      .rename(columns={'Symbol': 'gene', 'EnsemblID': 'gene_ensembl'})
)

df_master = (df_master
    .merge(gene_id_mapping[['gene', 'gene_ensembl']], how='left', on='gene')
    .sort_values(by=['cell_type', 'num_hits'])
    .reset_index(drop=True)
)

# %%
# write summarised data frame  
(df_master
    .to_csv(directory.joinpath('data/filtered_hit_list_summarised.csv'), index=False)
)
 
# %%
