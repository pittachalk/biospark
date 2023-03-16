#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %%
import pandas as pd
import numpy as np
import plotly.express as px
import tqdm
import re
from pathlib import Path

directory = Path(__file__).resolve().parent.parent

# read the compiled hit list
# setting the first column containing the broad cell types as the index
df = pd.read_csv(directory.joinpath('data/filtered_hit_list.csv'), index_col=0)

# %%
# for each broad cell type...
for i in tqdm.tqdm(range(df.shape[0])):

    # get occurences of genes 
    current_hits = df.iloc[i]
    current_hits = current_hits.dropna()

    # retrieve cell type, remove spaces and substituting slashes
    current_cell = df.index[i]
    current_cell = re.sub("/", '-', current_cell)
    current_cell = re.sub(" ", '', current_cell)

    # retrieve the number of times each gene was a hit as a dataframe
    genes, counts = np.unique(current_hits, return_counts=True)
    df_subset = pd.DataFrame({'gene': genes, 'num_hits': counts})
    df_subset['num_hits'] = df_subset['num_hits'].astype("category")

    # for each category of num_hits, join the genes as a single string
    tmp_df = (df_subset
        .groupby('num_hits')
        .agg({'gene': [lambda x: ', '.join(x), 'count']})
    )
    df_summary = (tmp_df
        .loc[:, "gene"]
        .rename(columns={'<lambda_0>': "Genes", 'count': "Frequency"})
        .reset_index()
    )

    # create interactive bar chart to visualise this data 
    fig = px.bar(df_summary, 
        x="num_hits", y="Frequency", 
        hover_data={'Genes': True, 'num_hits': False},
        template="seaborn"
    )
    fig.update_layout(hovermode='x', 
        xaxis_title="Number of Hits", 
        yaxis_title="Frequency",
        title_text=f"Gene hits in {df.index[i]} screens",
        title_font_size=30
    )
    fig.update_traces(showlegend=False)
    fig.write_html(directory.joinpath(f"script/figure/03_count_hits/hitHistogram_{current_cell}.html"))


# %%
