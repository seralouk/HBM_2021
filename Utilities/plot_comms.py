"""
=========================================
Author: Serafeim Loukas, May 2018
=========================================

"""
from numpy import *
from igraph import *

def plot_communities(X):
    X=X.values                                             
    X=X.tolist()
    g = Graph.Weighted_Adjacency(X, mode=ADJ_UNDIRECTED, attr="weight", loops = False)  
    N = g.vcount()
    E = g.ecount()
    #g.es["label"] = g.es["weight"]
    g.vs["name"] = range(N)

    layout = g.layout("kk")

    visual_style = {}
    visual_style["vertex_label"] = g.vs["name"]
   
    comms = g.community_fastgreedy(weights = g.es["weight"])
    clusters = comms.as_clustering()

    return clusters, visual_style, layout