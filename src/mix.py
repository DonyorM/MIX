from random import sample
import scipy.sparse as sparse
import networkx as nx
import pandas as pd

def MIX(adata, diffusion_level=2, input_obs='sample_labels', output_obs='mix_labels'):
    conns = adata.obsp['connectivities']
    input_labels = adata.obs[input_obs]
    mix_labels = pd.Series(0, input_labels.index, name=input_labels.name)
    G = nx.from_scipy_sparse_array(conns)

    def diffuse(from_node, current_level, value):
        for to_node in G[from_node]:
            mix_labels[to_node] += (value / current_level) / (len(G[to_node]) + 1)
            if current_level < diffusion_level:
                diffuse(to_node, current_level + 1, value)

    i = 0
    for node in G.nodes():
        value = input_labels.iloc[node]
        mix_labels.iloc[node] += value/(len(G[node]) + 1)
        diffuse(node, 1, value)
        i += 1
    
    adata.obs[output_obs] = mix_labels