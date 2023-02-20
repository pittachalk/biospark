#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: alymakhlouf
"""

import pandas as pd
import math

def hit_search(screen_list, x, group):
    
    hit_list = []
    alias_list = []
    
    for screen_id in group['#SCREEN_ID']:
        
        screen_match = [match for match in screen_list if 'SCREEN_' + str(screen_id) + '-' in match][0]
        df_screen = pd.read_csv(screen_match, sep='\t')

        num_hits = math.ceil(len(df_screen[df_screen['HIT'] == 'YES']) * (x/100))

        hit_list.append([i for i in df_screen[df_screen['HIT'] == 'YES']['OFFICIAL_SYMBOL']][0:num_hits])
        alias_list.append([i.split('|') for i in df_screen[df_screen['HIT'] == 'YES']['ALIASES']][0:num_hits])

    
    hit_list = [i for sublist in hit_list for i in sublist]
    alias_list = [i for sublist in alias_list for i in sublist]
    
    return hit_list, alias_list