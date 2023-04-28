#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dash import Dash, html, dash_table, dcc, callback, Output, Input
import pandas as pd
import numpy as np
import plotly.express as px
from pathlib import Path
import re

directory = Path(__file__).resolve().parent.parent

# function to wrap a comma-joined gene list by making a break after a few commas
def wrap_gene_names(text):
    replacement = "<br>"
    pattern = r"((?:[^,]*,){14}[^,]*),"  # match every 15th comma
    result = re.sub(pattern, lambda match: match.group(1) + replacement, text)
    return result

# read summarised data frame, renaming columns to be human friendly
df = pd.read_csv(directory.joinpath('data/filtered_hit_list_summarised.csv'))
df = df.rename(
    columns = {
        'cell_type': "Cell Type", 
        'gene': 'Gene', 
        'num_hits': 'Number of Hits',
        'num_screens': "Number of Screens",
        'gene_ensembl': "Ensembl ID"
    }
)

df_intersect = pd.read_table(directory.joinpath('data/pairwise_intersects.txt'))
df_intersect = df_intersect.rename(
    columns = {
        'cell_type1': "Cell Type",
        'cell_type2': "Compared Cell Type",
        'n_intersect': "Number of Intersects",
        'n_genes1': "Number of Genes",
        'intersect_list': "Genes in Common"
    }
)
df_intersect["Genes in Common"] = [wrap_gene_names(str(x)) for x in df_intersect["Genes in Common"]]

app = Dash(__name__)

app.layout = html.Div([
    html.H1(children='Top Gene Hits in BioGRID CRISPR Screens', style={'textAlign':'center'}),
    html.H2(children='Choose cell type', style={'textAlign':'left'}),
    dcc.Dropdown(df['Cell Type'].unique(), 'Breast', id='dropdown-selection'),
    dcc.Graph(figure={}, id='graph-content-genehits'),
    dash_table.DataTable(
        page_size=12, 
        id='bar_chart',
        sort_action="native"
    ),
    html.H2(children='Comparisons to other cell types (double click to zoom out)', style={'textAlign':'left'}),
    dcc.Graph(figure={}, id='graph-content-intersect'),
    # html.P(children='Here is some dummy text.', style={'textAlign':'left'}),
    html.Hr()
])

##########

@callback(
    Output('graph-content-genehits', 'figure'),
    Input('dropdown-selection', 'value')
)
def update_graph_genehits(value):
    # for each num_hits, join the genes as a single string and count genes
    df_agg = (df[df['Cell Type']==value]
        .groupby('Number of Hits')
        .agg({'Gene': [lambda x: ', '.join(x), 'count']})
        .loc[:, "Gene"]
        .rename(columns={'<lambda_0>': "Genes", 'count': "Frequency"})
        .reset_index()
    )

    # wrap "Genes" text by making a break after a few commas
    def wrap_gene_names(text):
        replacement = "<br>"
        pattern = r"((?:[^,]*,){9}[^,]*),"  # match every 10th comma
        result = re.sub(pattern, lambda match: match.group(1) + replacement, text)
        return result
    df_agg['Genes'] = [wrap_gene_names(x) for x in df_agg['Genes']]

    # create interactive bar chart to visualise this data 
    fig = px.bar(df_agg, x="Number of Hits", y="Frequency", 
                 hover_data={'Genes': True, 'Number of Hits': False},
                 template="seaborn")
    fig.update_layout(hovermode='x', xaxis_title="Number of Hits", yaxis_title="Frequency")
    fig.update_traces(showlegend=False)
    return fig

##########

@callback(
    Output('bar_chart', 'data'),
    Input('dropdown-selection', 'value')
)
def update_rows(value):
    dff = df[df['Cell Type']==value]
    return dff.sort_values("Number of Hits", ascending=False).to_dict('records')

##########

@callback(
    Output('graph-content-intersect', 'figure'),
    Input('dropdown-selection', 'value')
)
def update_graph_intersect(value):
    tmp_df = (df_intersect
        .query("`Cell Type`==@value & `Compared Cell Type`!=@value")
        .copy()
        .sort_values("Number of Intersects", ascending=False)      
    )

    # create interactive bar chart to visualise this data 
    fig = px.bar(tmp_df, #.iloc[0:20],  # show top hits
                x="Compared Cell Type", y="Number of Intersects", 
                hover_data={'Genes in Common': True, 'Number of Intersects': True},
                template="seaborn")
    fig.update_layout(hovermode='x', yaxis_title="Genes in Common")
    fig.update_traces(marker_color='#4e7846')
    fig.update_xaxes(range=[-0.5, 25.5])
    return fig

##########

if __name__ == '__main__':
    app.run_server(debug=True)

