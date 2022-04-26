from random import sample
import networkx as nx
import numpy as np


def mix(
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
