import numpy as np
import pandas as pd
from . import distances

def mad(tree):
    """Run Minimal Ancestral Deviation Algorithm

    Parameters
    ----------
    tree : Dendropy Tree
        Tree to reroot

    Returns
    -------
    tree : Dendropy Tree
        Tree rooted at branch with minimal ancestral deviation
    RAI : float
        root ambiguity index
    """
    # -------------------------------------------
    # 1) Calculate pairwise distances
    # -------------------------------------------
    pair_matrix = distances.pairwise_distance_matrix(tree)
    dist_method = tree.node_distance_matrix()

    n_branches = len(tree.edges())
    branches = pd.Series(index=tree.edges(), dtype=float)
    rhos = pd.Series(index=tree.edges(), dtype=float)

    # ---------------------------------------------------------
    # 2) Calculate root-mean-square deviations for each branchs
    # ---------------------------------------------------------
    for branch in tree.edges():
        # -------------------------------------------
        # A) Split pairs of nodes for this branch
        # -------------------------------------------
        if branch.length == None:
            continue

        dij = branch.length
        head = branch.head_node

        # Split all pairs of nodes in to straddling or same-side pairs
        pairs_strad, pairs_same = distances.split_node_pairs(branch, tree)

        # Initialize the deviation score array for current branch
        n_strad = len(pairs_strad)
        n_same = len(pairs_same)
        rbca = np.empty(n_strad+n_same, dtype=float)

        # -------------------------------------------
        # B) Calculate deviation for straddling pairs
        # -------------------------------------------

        # Distance between two leaf nodes
        dbc = np.array([pair_matrix[n1.label][n2.label] for n1, n2 in pairs_strad])

        # Distance from leaf node to candidate branch
        dbi = np.array([dist_method(head,n1) for n1, n2 in pairs_strad])

        # Find optimal position for node on branch
        rho = np.sum((dbc-2*dbi)*dbc**-2)/(2*dij*np.sum(dbc**-2))
        rho = min(max(0, rho), 1)
        rhos[branch] = rho

        # Calculate distance from leaf to branch
        dab = dbi + dij*rho
        rbca[:n_strad] = np.abs(2*dab/dbc-1)

        # -------------------------------------------
        # C) Calculate deviation for same side pairs
        # -------------------------------------------

        # Calculate deviation for same-side node pairs
        dbc_same = np.array([pair_matrix[n1.label][n2.label] for n1, n2 in pairs_same])
        dab_same = np.array([dist_method(head,n1) for n1, n2 in pairs_same])
        rbca[n_strad:] = np.abs(2*dab_same/dbc_same-1)

        # -------------------------------------------
        # D) Calculate root-mean-square deviation for branch
        # -------------------------------------------

        # Calculate the root-mean-square.
        r = np.sqrt(np.mean(rbca**2))
        branches[branch] = r

    # ---------------------------------------------------------
    # 3) Find minimal ancestral deviation branch.
    # ---------------------------------------------------------
    sorted_branches = branches.sort_values(ascending=False)

    # Get best branch from MAD list
    ideal_branch = sorted_branches.index[0]
    ideal_branch2 = sorted_branches.index[1]

    # find position of new root node.
    dij = ideal_branch.length
    length1 = rhos[ideal_branch] *  dij
    length2 = dij - length1

    # ---------------------------------------------------------
    # 4) Reroot tree.
    # ---------------------------------------------------------
    tree.reroot_at_edge(ideal_branch, length1=length1, length2=length2, update_bipartitions=True)

    # ---------------------------------------------------------
    # 5) Calculate root ambiguity index
    # ---------------------------------------------------------
    RAI = branches[ideal_branch] / branches[ideal_branch2]

    return tree, RAI
