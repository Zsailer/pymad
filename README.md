# Python implementation of MAD (minimal ancestor deviation)

This a re-implementation of the MAD (minimal ancestor deviation) algorithm
for rooting phylogenetic trees, published in Nature Ecology and Evolution.

See the paper [here](https://www.nature.com/articles/s41559-017-0193.epdf?shared_access_token=S62ZDIpEuBo7f8f8fO56xtRgN0jAjWel9jnR3ZoTv0Px0yVdsafzuduOQbkT4JkJOHFGG1kSo-AkPiJ94m3CK2Xm6hLEVZRd0qsbm2Zk_ZJCBMBoq9NRfBs3I65bR2aj3uJttkhIZL7CzhQnRstgJr_2jMXVxFFvQLZOXTo9zmw%3D).

The original implementations were done in R and Matlab. You can download that code
[here](https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen).

## Basic Usage

Using the command line interface:
```bash
pymad -i tree.nwk -o rooted-tree.nwk --schema "newick" # schema is optional
```

Using the API:
```python
import pymad

# Load the tree and quality control checks
tree = pymad.load_tree("tree.nwk")

# Root using the MAD algorithm
new_tree, RAindex = spymad.mad(tree)
```

## Install

Clone this repo and install in development:
```
pip install -e .
```

## Credit

All credit goes to the authors of the paper in Nature: Ecology and Evolution. This
is merely a translation of their work into Python. Please cite them if you use this
software!
