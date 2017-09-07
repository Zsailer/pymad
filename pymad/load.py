import dendropy as dp
from collections import Counter

def load(path, schema="newick"):
    """Read a tree from file and run quality controls.
    """
    tree = dp.Tree.get(path=path, schema=schema)

    # Remove any redundant leafs
    tree = remove_redundant_taxa(tree)

    # Label ancestor nodes
    tree = check_nodes_are_labelled(tree)
    return tree

def remove_redundant_taxa(tree):
    """Checks if tree contains redundant leaf nodes. If so, removes them.
    """
    leafs = tree.leaf_nodes()
    leaf_labels = tree.taxon_namespace.labels()

    # Check that all leafs are unique
    total = 0
    for label, count in Counter(leaf_labels).items():
        total += count-1
        for i in range(count-1):
            node = tree.find_node_with_taxon_label(label)
            tree.prune_nodes([node])

    if total > 0: print("{} taxa were redundant, so I removed them.".format(total))
    return tree

def check_nodes_are_labelled(tree):
    """Labels all nodes in the tree. Labels ancestral nodes.
    """
    internal = tree.internal_nodes()
    for i, node in enumerate(internal):
        node.label = "anc{}".format(i)
    for node in tree.leaf_nodes():
        node.label = node.taxon.label
    return tree
