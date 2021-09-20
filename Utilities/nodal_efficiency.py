"""
=========================================
Author: Serafeim Loukas, May 2018
=========================================

"""
from numpy import *
from igraph import *
seterr(divide='ignore')


def nodal_eff(g):
    """
    This function calculates the nodal efficiency of a weighted graph object.
    
    Created by: Loukas Serafeim, Nov 2017

    Args:
     g: A igraph Graph() object.

    Returns:
     The nodal efficiency of each node of the graph
    """
    weights = g.es["weight"][:]
    weights = [1.0 / x for x in weights]
    sp = (1.0 / array(g.shortest_paths_dijkstra(weights=weights)))
    fill_diagonal(sp,0)
    N=sp.shape[0]
    ne= (1.0/(N-1)) * apply_along_axis(sum,0,sp)

    return ne