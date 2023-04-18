#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dash import Dash, html, dash_table, dcc, callback, Output, Input
import pandas as pd
import numpy as np
import plotly.express as px
from pathlib import Path
import re

directory = Path(__file__).resolve().parent.parent

# read summarised data frame, renaming columns to be human friendly
df = pd.read_csv(directory.joinpath('data/filtered_hit_list_summarised.csv'))
df = df.rename(
    columns={
        'cell_type': "Cell Type", 
        'gene': 'Gene', 
        'num_hits': 'Number of Hits',
        'num_screens': "Number of Screens",
        'gene_ensembl': "Ensembl ID"
    }
)

app = Dash(__name__)

app.layout = html.Div([
    html.H1(children='Top Gene Hits in BioGRID CRISPR Screens', style={'textAlign':'center'}),
    html.H2(children='Choose cell type', style={'textAlign':'left'}),
    dcc.Dropdown(df['Cell Type'].unique(), 'Breast', id='dropdown-selection'),
    dcc.Graph(figure={}, id='graph-content'),
    dash_table.DataTable(
        page_size=12, 
        id='bar_chart',
        sort_action="native"
    ),
    # html.P(children='Here is some dummy text.', style={'textAlign':'left'}),
    html.Hr()
])

@callback(
    Output('graph-content', 'figure'),
    Input('dropdown-selection', 'value')
)
def update_graph(value):
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

@callback(
    Output('bar_chart', 'data'),
    Input('dropdown-selection', 'value')
)
def update_rows(value):
    dff = df[df['Cell Type']==value]
    return dff.sort_values("Number of Hits", ascending=False).to_dict('records')

if __name__ == '__main__':
    app.run_server(debug=True)