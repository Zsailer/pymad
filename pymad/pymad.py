edges = tree.edges()
M_dist = tree.node_distance_matrix()

# Test on single edge
edge = edges[9]

i = edge.head_node
j = edge.tail_node

I = set(i.child_nodes()).intersection(i.leaf_nodes())
J = set(j.child_nodes()).intersection(j.leaf_nodes())

dij = edge.length

dbc = np.empty((len(I), len(J)), dtype=float)
dbi = np.empty((len(I), len(J)), dtype=float)

for index0, b in enumerate(I):
    for index1, c in enumerate(J):
        dbc[index0, index1] = M_dist(b,c)
        dbi[index0, index1] = M_dist(b,i)

rho = np.sum((dbc - 2 * dbi) * dbc**-2) /  (2 * dij * np.sum(dbc**-2))
rho = min(max(0, rho), 1)

r = np.sum( ((2*(dbi + rho * dij))/dbc  - 1)**2)

print(r)
