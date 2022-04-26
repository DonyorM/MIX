# MIX

MIX is a tool to help determine prototypical cells in a single-cell dataset with multiple states (such as control and experimental groups)

## Installation

At the moment MIX isn't deployed anywhere, the simplest method of installation is probably to copy the MIX file to your project's workspace. It depends on numpy and networkx (see the requirements.txt file in this repo)

## Usage

MIX requires a knn graph defined in the 'connectivities' section of the observation of the adata passed to it. As such it can be used like this:

```
import scanpy as sc
import mix

sc.pp.neighbors(adata, n_neighbors= 10)
mix.mix(adata, num_steps = 10000, alpha=0.15, input_obs="genotype", output_obs="mix_labels")
print(adata.obs['mix_labels'])
```

## License

MIX is licensed under the MIT license
