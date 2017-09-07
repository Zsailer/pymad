import itertools as it
import pandas as pd

def pairwise_distance_matrix(tree):
    """Construct a Pandas dataframe with pairwise distances between all
    nodes in tree.
    """
    dist = tree.node_distance_matrix()
    nodes = tree.nodes()
    node_labels = [n.label for n in nodes]

    matrix = pd.DataFrame(
        index=node_labels,
        columns=node_labels)

    for i in range(len(nodes)):
        for j in range(len(nodes)):
            node1 = nodes[i]
            node2 = nodes[j]
            matrix[node1.label][node2.label] = dist(node1, node2)
    return matrix

def split_node_pairs(branch, tree):
    """Generates two separate lists of node-tuple-pairs. The first list
    includes all pairs of nodes that straddle the branch. The second list includes
    pairs of nodes that are on the same side the given branch.
    """
    # Bundles nodes an each side of the branch
    dij = branch.length
    head = branch.head_node

    # Get all nodes on leftmost side of branch
    side_one = []
    # Preorder traverse the tree, and stop once we reach the edge.
    for n in head.preorder_iter(lambda x: x.is_leaf()):
        if n == branch.tail_node:
            break
        elif n.is_leaf():
            side_one.append(n)

    # Get nodes on right side
    side_two = list(set(tree.leaf_nodes()).difference(side_one))

    # Enumerate pairwise combinations for straddling nodes
    pairs_strad = [(n1, n2) for n1 in side_one for n2 in side_two]

    # Enumerate pairwise combinations for same side nodes
    pairs_side1 = list(it.combinations(side_one, 2))
    pairs_side2 = list(it.combinations(side_two, 2))
    pairs_same_side = pairs_side1 + pairs_side2

    return pairs_strad, pairs_same_side
