#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: alymakhlouf
"""

import pandas as pd
import math
import tqdm

def hit_search(screen_list, x, group):
    
    hit_list = []
    alias_list = []
    
    for screen_id in group['SCREEN_ID']:
        
        screen_match = [match for match in screen_list if 'SCREEN_' + str(screen_id) + '-' in match][0]
        df_screen = pd.read_csv(screen_match, sep='\t')

        num_hits = math.ceil(len(df_screen[df_screen['HIT'] == 'YES']) * (x/100))

        hit_list.append([i for i in df_screen[df_screen['HIT'] == 'YES']['OFFICIAL_SYMBOL']][0:num_hits])
        alias_list.append([i.split('|') for i in df_screen[df_screen['HIT'] == 'YES']['ALIASES']][0:num_hits])

    
    hit_list = [i for sublist in hit_list for i in sublist]
    alias_list = [i for sublist in alias_list for i in sublist]
    
    return hit_list, alias_list

def filter_essential_genes(directory, df_hit):
    
    # FIND LISTS OF ESSENTIAL AND NON-ESSENTIAL GENES FROM DepMap
    non_essential_genes = pd.read_csv(directory.joinpath('Achilles/AchillesNonessentialControls.csv'))['Gene'].tolist()
    non_essential_genes, non_essential_gene_ids = [x.split(' ')[0] for x in non_essential_genes], [x.split(' ')[1] for x in non_essential_genes]
    non_essential_gene_ids = [x[1:-1] for x in non_essential_gene_ids]
    
    df_essential_common = pd.read_csv(directory.joinpath('Achilles/AchillesCommonEssentialControls.csv'))
    df_essential_infer = pd.read_csv(directory.joinpath('Achilles/CRISPRInferredCommonEssentials.csv')).rename(columns={'Essentials':'Gene'})
    #df_essential_intersect = pd.merge(df_essential_common, df_essential_infer, how='inner', on='Gene')
    essential_genes = pd.concat([df_essential_common, df_essential_infer]).drop_duplicates().sort_values(by=['Gene'])['Gene'].tolist()
    essential_genes, essential_gene_ids = [x.split(' ')[0] for x in essential_genes], [x.split(' ')[1] for x in essential_genes]
    essential_gene_ids = [x[1:-1] for x in essential_gene_ids]
    
    filtered_cell_hit_dict = {}
    
    for i in tqdm.tqdm(list(range(len(df_hit)))[0:]):
    
        hit_list = df_hit.iloc[i][2:].unique()[0:-1] # get list of unique hits, excluding 'nan' value
        filtered_hit_list = [x for x in hit_list if x not in essential_genes]
        
        filtered_cell_hit_dict[df_hit.index[i]] = filtered_hit_list
    
    df_hit_filtered = pd.DataFrame.from_dict(filtered_cell_hit_dict, orient='index')
        
    return df_hit_filtered
