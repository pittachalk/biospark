# %%
import pandas as pd
import numpy as np
import plotly.express as px
import tqdm
from pathlib import Path

# %%
directory = Path(__file__).resolve().parent.parent.joinpath("data/")



# %%
df = pd.read_csv(directory.joinpath('hit_list.csv'), index_col=0)




# %%
for i in tqdm.tqdm(range(df.shape[0])):
    current_hits = df.iloc[i]
    current_hits = current_hits.dropna()

    current_cell = df.index[i]

    genes, counts = np.unique(current_hits, return_counts=True)

    df_subset = pd.DataFrame({'gene': genes, 'count': counts})

    df_subset['count'] = df_subset['count'].astype("category")

    tmp_df = (df_subset.groupby('count')
            .apply(lambda dataframe: ', '.join(dataframe.gene))
            .reset_index()
            .rename(columns={0: "Genes"})
    )

    df_subset = pd.merge(df_subset, tmp_df)


    fig = px.bar(df_subset, 
        x="count", color='Genes', hover_data={'count': False}
    )
    fig.update_layout(hovermode='x', 
        xaxis_title="Number of Hits", yaxis_title="Frequency")
    fig.update_traces(showlegend=False)
    fig.write_html(directory.joinpath(f"figure/hitHistogram_{current_cell}.html"))


# %%
df_subset = pd.DataFrame({'gene': genes, 'count': counts})

tmp_df = (df_subset.groupby('count')
          .apply(lambda dataframe: ', '.join(dataframe.gene))
          .reset_index()
          .rename(columns={0: "Genes"})
)

df_subset = pd.merge(df_subset, tmp_df)


# %%
fig = px.histogram(df_subset, 
    x="count", color='Genes', hover_data={'count': False}
)
fig.update_layout(hovermode='x', 
    xaxis_title="Number of Hits", yaxis_title="Frequency")
fig.update_traces(showlegend=False)
fig.write_html(directory.joinpath(f"figure/hitHistogram_{current_cell}.html"))

# %%
