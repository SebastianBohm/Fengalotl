import glasbey
import numpy as np

import plotly.graph_objects as go
from fengalotl.fct.load import get_data
from fengalotl._constants import GENES_LABEL

def plot_space(input):

    if not get_data(input):
        return None

    fig = go.Figure()
    fig.update_layout(
        template="plotly_dark",
        showlegend=True,
        autosize=True,
        yaxis_scaleanchor="x",
        xaxis=dict(title="x"),
        yaxis=dict(title="y")
    )
    
    return fig

def add_space_clusters(input, widget):

    if not (adata := get_data(input)):
        return None
    
    if not input.select_resolution():
        return None

    widget.data = [trace for trace in widget.data if not trace.name.isdigit() and trace.name != 'cluster']

    if input.switch_clusters():

        resolution = input.select_resolution()
        cluster_ids = np.unique(adata.obs[resolution].dropna())
        colors = glasbey.create_palette(palette_size=len(cluster_ids), lightness_bounds=(50, 100))
        
        # Get spatial coordinates
        x_coords = adata.obsm['spatial'][:, 0]
        y_coords = adata.obsm['spatial'][:, 1]

        for color, cluster_id in zip(colors, cluster_ids):
            mask = adata.obs[resolution] == cluster_id
            widget.add_trace(
                go.Scatter(
                    x=x_coords[mask],
                    y=y_coords[mask],
                    name=str(cluster_id),
                    mode='markers',
                    marker=dict(
                        color=color,
                        size=input.slider_dotsize_space()
                    )
                )
            )

def add_space_expression(input, widget):

    if not (adata := get_data(input)):
        return None

    widget.data = [trace for trace in widget.data if trace.name not in GENES_LABEL]

    if input.switch_expression() and input.select_gene():

        gene_name = input.select_gene()
        if gene_name in adata.var_names:
            gene_expression = np.array(adata[:,gene_name].X.flatten())
            x_coords = adata.obsm['spatial'][:, 0]
            y_coords = adata.obsm['spatial'][:, 1]
            
            widget.add_trace(
                go.Scatter(
                    name = gene_name,
                    x = x_coords,
                    y = y_coords,
                    mode='markers',
                    marker=dict(
                        color = gene_expression,
                        colorscale = 'Viridis',
                        showscale = True,
                        size = input.slider_dotsize_space(),
                        colorbar=dict(
                            orientation = 'h',
                            lenmode='fraction',
                            len=0.25,
                            thickness=10,
                            y = 0.05,
                            x = 0.15)
                    ),
                )
            )
