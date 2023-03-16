#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: alymakhlouf
"""

import pandas as pd
import glob
import tqdm
from pathlib import Path

from BioGRID import hit_search, filter_essential_genes

directory = Path(__file__).resolve().parent.parent.joinpath("data/")

#directory = '/Users/alymakhlouf/Desktop/biospark/data/'

# SELECT ONLY TIMECOURSE, INHIBITION/KNOCKOUT EXPERIMENTS

df = pd.read_csv(directory + 'index_file_polished.csv')
df = df[df['EXPERIMENTAL_SETUP'] == 'Timecourse']
df = df[df['LIBRARY_METHODOLOGY'] != 'Activation']

# EXTRACT UNIQUE AND NON-UNIQUE, SINGLE AND MULTIPLE CELL TYPES


cell_type_broad = [i for i in df.CELL_TYPE_BROAD.unique()]
cell_type_broad_single = [i for i in set([i for j in [i.split(' / ') 
                    for i in df.CELL_TYPE_BROAD.unique()] for i in j])]
cell_type_broad_multiple = [i for i in cell_type_broad if '/' in i]
cell_type_broad_all = cell_type_broad_single + cell_type_broad_multiple

cell_type_broad_non_unique = [i for i in [j for j in cell_type_broad_single 
                          for k in cell_type_broad_multiple if j in k]]
cell_type_broad_non_unique_subset = [x for x in set(cell_type_broad_non_unique) 
                          if cell_type_broad_non_unique.count(x) == 1]
cell_type_broad_non_unique = [i for i in set([j for j in cell_type_broad_single 
                          for k in cell_type_broad_multiple if j in k])]    

screen_list = glob.glob(directory + "BIOGRID-ORCS-ALL-homo_sapiens-1.1.13.screens/*.txt")


# GROUP CELLS BASED ON UNIQUE (SINGLE & MULTIPLE) BROAD CELL TYPE


x = 10 # filter top x% of each hit list

cell_hit_dict = {}
cell_alias_dict = {}

number_screens_dict = {} # track number of screens per broad cell type - needed for normalisation


for cell_broad in tqdm.tqdm(cell_type_broad_all[0:]):
    
    cell_broad_key = cell_broad
    
    print(cell_broad_key)   
    
    group = df[df['CELL_TYPE_BROAD'] == cell_broad]
    
    if len(group) == 0:
        continue
    
    cell_hit_dict[cell_broad_key], cell_alias_dict[cell_broad_key] = hit_search(screen_list,x,group)
    number_screens_dict[cell_broad_key] = len(group)


# GROUP CELLS BASED ON NON-UNIQUE BROAD CELL TYPE


for cell_broad in tqdm.tqdm(cell_type_broad_non_unique[0:]):
    
    cell_broad_key = cell_broad + '*'
    
    print(cell_broad_key)   

    group = df[df['CELL_TYPE_BROAD'].str.contains(cell_broad)]
    
    
    cell_hit_dict[cell_broad_key], cell_alias_dict[cell_broad_key] = hit_search(screen_list,x,group)

    # delete redundant (duplicated) entries
    if cell_broad in cell_type_broad_non_unique_subset:
        cell_hit_dict.pop(cell_broad_key)
        cell_alias_dict.pop(cell_broad_key)
        continue
    
    number_screens_dict[cell_broad_key] = len(group)

df_hit = pd.DataFrame.from_dict(cell_hit_dict, orient='index')
df_alias = pd.DataFrame.from_dict(cell_alias_dict, orient='index')
df_number_screens = pd.DataFrame.from_dict(number_screens_dict, orient='index').rename(columns={0:'# Screens'})

df_hit = pd.concat([df_number_screens, df_hit], axis=1)
df_alias = pd.concat([df_number_screens, df_alias], axis=1)

df_hit.to_csv(directory + 'unfiltered_hit_list.csv')        
df_alias.to_csv(directory + 'unfiltered_alias_list.csv')


# FILTER OUT ESSENTIAL GENES

df_hit_filtered = filter_essential_genes(directory, df_hit)
df_hit_filtered.to_csv(directory + 'filtered_hit_list.csv')
