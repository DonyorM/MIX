from random import sample
import scipy.sparse as sparse
import networkx as nx
import pandas as pd
import numpy as np


def MIX_daniel(
    adata, diffusion_level=2, input_obs="sample_labels", output_obs="mix_labels_daniel"
):
    conns = adata.obsp["connectivities"]
    input_labels = adata.obs[input_obs]
    mix_labels = pd.Series(0, input_labels.index, name=input_labels.name)
    G = nx.from_scipy_sparse_array(conns)

    def diffuse(from_node, current_level, value):
        for to_node in G[from_node]:
            mix_labels[to_node] += (value / current_level) / (len(G[to_node]) + 1)
            if current_level < diffusion_level:
                diffuse(to_node, current_level + 1, value)

    for node in G.nodes():
        value = input_labels.iloc[node]
        mix_labels.iloc[node] += value / (len(G[node]) + 1)
        diffuse(node, 1, value)

    adata.obs[output_obs] = mix_labels


def MIX_monty(
    adata, num_steps, alpha=0.2, input_obs="sample_labels", output_obs="mix_labels"
):
    mix_labels = adata.obs[input_obs].astype(np.float32)
    G = nx.from_scipy_sparse_array(adata.obsp["connectivities"])

    steps_remaining = num_steps
    while steps_remaining > 0:
        if steps_remaining > len(G.nodes):
            this_steps = len(G.nodes)
        else:
            this_steps = steps_remaining

        steps_remaining -= this_steps

        for seed_node in sample(G.nodes, this_steps):
            for neighbour in G[seed_node]:
                flow = (mix_labels[seed_node] - mix_labels[neighbour]) * alpha
                mix_labels[neighbour] += flow

    adata.obs[output_obs] = mix_labels
