#!/usr/bin/env python3
import networkx as nx


def groupAlignmentDict(D, namD):
    G = nx.Graph()
    for val in namD.values():
        G.add_node(val)
    for key in D:
        pairings = D[key].keys()
        for pairing in pairings:
            G.add_edge(key, pairing)
    return (list(nx.connected_components(G)))
