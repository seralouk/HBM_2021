"""
=========================================
Author: Serafeim Loukas, Oct 2018
=========================================
"""
from numpy import *
from igraph import *
import string
from nodal_efficiency import nodal_eff
import time 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


###############################################################################################
###############################################################################################

def extract_global(subject_folder_names, path_FC_matrices, filename_of_FC):
  """
  This function creates an empty dataframe, builds the brain graphs,
  calculates the network GLOBAL features and writes them in the dataframe.
  Created by: Loukas Serafeim, Oct 2018

  Args:
    subject_folder_names: A list with subject numbers. e.g. ['Subject_01']
    path_FC_matrices: path to folder that contains the connectomes
    filename_of_FC: a list of string with the name of the connectomes e.g. 'PRE_pearson_{}.csv'

  Returns:
    3 Global network features: deg, cc, ne for all subjects in one dataframe
  """

  dt1 = pd.DataFrame(columns=['Deg', 'Cc', 'Ne'])
  for s in subject_folder_names:
    X=pd.read_csv(path_FC_matrices + filename_of_FC.format(s),header=None)
    X=X.values
    X=X.tolist()
    g = Graph.Weighted_Adjacency(X, mode=ADJ_UNDIRECTED, attr="weight", loops = False)
    N=g.vcount()
    E=g.ecount()
    g.es["label"] = g.es["weight"]
    lst = range(N)
    g.vs["name"] = lst

    d = g.strength(loops=False, weights=g.es["weight"])
    av_d = average(d)  
    cc = g.transitivity_local_undirected(vertices=None,mode="zero", weights=g.es["weight"])
    av_cc = average(cc)    
    ne = nodal_eff(g)
    av_ne = average(ne)

    dt1.loc[s] =[av_d,av_cc,av_ne]

  return dt1

###############################################################################################
###############################################################################################

def extract_nodal(subject_folder_names, path_FC_matrices, filename_of_FC):
  """
  This function creates an empty dataframe, builds the brain graphs,
  calculates the network NODAL features and writes them in the dataframe.
  Created by: Loukas Serafeim, Oct 2018

  Args:
    subject_folder_names: A list with subject numbers. e.g. ['Subject_01']
    path_FC_matrices: path to folder that contains the connectomes
    filename_of_FC: a list of string with the name of the connectomes e.g. 'PRE_pearson_{}.csv'

  Returns:
    3  network features per node: deg, cc, ne for all subjects in one dataframe
  """
  counter = 0
  for s in subject_folder_names:
    X=pd.read_csv(path_FC_matrices + filename_of_FC.format(s),header=None)
    X=X.values
    np.fill_diagonal(X,0)
    X=X.tolist()
    g = Graph.Weighted_Adjacency(X, mode=ADJ_UNDIRECTED, attr="weight", loops = False)
    N=g.vcount()
    E=g.ecount()
    g.es["label"] = g.es["weight"]
    lst = range(N)
    g.vs["name"] = lst

    d = g.strength(loops=False, weights=g.es["weight"]) 
    cc = g.transitivity_local_undirected(vertices=None,mode="zero", weights=g.es["weight"])   
    ne = nodal_eff(g)

    d=np.array(d)
    cc=np.array(cc)

    if counter == 0:
      dt1 = pd.DataFrame(columns=['Deg']*d.shape[0] + ['Cc']*d.shape[0] + ['Ne']*d.shape[0])

    all_together = np.concatenate((d,cc,ne), axis = 0)
    dt1.loc[s] = all_together[:]
    counter +=1

  return dt1

###############################################################################################
###############################################################################################


def extract_connections(subject_folder_names, path_FC_matrices, filename_of_FC, number_of_ROIs):
  """
  This function creates an empty dataframe, builds the brain graphs,
  calculates the CONENCTION features and return them in a numpy array.
  Created by: Loukas Serafeim, Oct 2018

  Args:
    subject_folder_names: A list with subject numbers. e.g. ['Subject_01']
    path_FC_matrices: path to folder that contains the connectomes
    filename_of_FC: a list of string with the name of the connectomes e.g. 'PRE_pearson_{}.csv'

  Returns:
   The connectomes of the CTRL stored in a tensor X
  """
  X = np.zeros((number_of_ROIs, number_of_ROIs, len(subject_folder_names)))

  for s, n in enumerate(subject_folder_names):
    XX=pd.read_csv(path_FC_matrices + filename_of_FC.format(n),header=None)
    #print(XX.shape)
    for se in range(X.shape[0]):
      XX.iloc[se,se] = 0.0
    X[:,:,s-1] = XX[:]
  
  X_connections = X[np.triu_indices(X.shape[0], k = 1)]
  X_connections = X_connections.T

  return X_connections

###############################################################################################
###############################################################################################

def Build_mean_network(subject_folder_names, path_FC_matrices, filename_of_FC, number_of_ROIs):
  """
  This function calculates the Mean network. The decomposition of this network
  will be used as prior information for the decomposition of all groups/subjects.
  
  Created by: Loukas Serafeim, Oct 2018

  Args:
    subject_folder_names: A list with subject numbers. e.g. ['Subject_01']
    path_FC_matrices: path to folder that contains the connectomes
    filename_of_FC: a list of string with the name of the connectomes e.g. 'PRE_pearson_{}.csv'
    number_of_ROIs: number

  Returns:
    The mean network in a pandas dataframe
  """
  X = np.zeros((number_of_ROIs, number_of_ROIs, len(subject_folder_names)))
  for i, s in enumerate(subject_folder_names):
    XX=pd.read_csv(path_FC_matrices + filename_of_FC.format(s),header=None)
    for se in range(X.shape[0]):
      XX.iloc[se,se] = 0.0
    X[:,:,i-1] = XX
    Xout = np.sum(X, axis=2) / float(len(subject_folder_names))
    Xout = pd.DataFrame(Xout)

  return Xout

###############################################################################################
###############################################################################################
  
def decompose_fg(X, subject_folder_names, path_FC_matrices, filename_of_FC, verbose=0):
  """
  This function decomposes the inputed subjects using the fast greedy algorithm
  based on the prior information (decomposition of X, X= mean network).
  
  Created by: Loukas Serafeim, Oct 2018

  Args:
    X: mean netwrok to be used as prior information
    subject_folder_names: A list with subject numbers. e.g. ['Subject_01']
    path_FC_matrices: path to folder that contains the connectomes
    filename_of_FC: a list of string with the name of the connectomes e.g. 'PRE_pearson_{}.csv'


  Returns:
   The average degree, clustering coefficient and nodal efficiency of each community
   (3*comms features)
  """

  X=X.values 
  np.fill_diagonal(X,0)
  X=X.tolist()

  g = Graph.Weighted_Adjacency(X, mode=ADJ_UNDIRECTED, attr="weight", loops = False)  
  N = g.vcount()
  E = g.ecount()
  g.es["label"] = g.es["weight"]
  g.vs["name"] = range(N)
 
  comms = g.community_fastgreedy(weights = g.es["weight"])
  clusters = comms.as_clustering()
  if verbose:
    print("\nFast greedy: after calculating the average network the comms are:  ")
    print(clusters)

  c =['deg', 'cc','ne'] * len(clusters)
  dt_fg = pd.DataFrame(columns = c)

  for j , n in enumerate(subject_folder_names):
    Xc = pd.read_csv(path_FC_matrices + filename_of_FC.format(n),header=None)
    Xc = Xc.values
    np.fill_diagonal(Xc,0)
    Xc = Xc.tolist()
    gc = Graph.Weighted_Adjacency(Xc, mode=ADJ_UNDIRECTED,attr="weight", loops=False)
    Nc = gc.vcount()
    Ec = gc.ecount()
    gc.es["label"] = gc.es["weight"]
    lst_c = range(Nc)
    gc.vs["name"] = lst_c

    lists = []

    for i in range(len(clusters)):
      sub = gc.induced_subgraph(clusters[i])
      av_d = average(sub.strength(loops=False,weights=sub.es["weight"]))
      av_cc = average(sub.transitivity_local_undirected(vertices=None,mode="zero", weights=sub.es["weight"]))
      av_ne =average(nodal_eff(sub))

      lists.extend([av_d,av_cc,av_ne])
    dt_fg.loc[j] = lists

  return dt_fg, clusters

###############################################################################################
###############################################################################################
 
def decompose_lev(X, subject_folder_names, path_FC_matrices, filename_of_FC, verbose=0):
  """
  This function decomposes the inputed subjects using the leading eigenvector algorithm
  based on the prior information (decomposition of X, X= mean network).
  
  Created by: Loukas Serafeim, Oct 2018

  Args:
    X: mean netwrok to be used as prior information
    subject_folder_names: A list with subject numbers. e.g. ['Subject_01']
    path_FC_matrices: path to folder that contains the connectomes
    filename_of_FC: a list of string with the name of the connectomes e.g. 'PRE_pearson_{}.csv'

  Returns:
   The average degree, clustering coefficient and nodal efficiency of each community
   (3*comms features)
  """
  X=X.values
  np.fill_diagonal(X,0)                                             
  X=X.tolist()

  g = Graph.Weighted_Adjacency(X,mode=ADJ_UNDIRECTED,attr="weight", loops = False)  
  N = g.vcount()
  E = g.ecount()
  g.es["label"] = g.es["weight"]
  g.vs["name"] = range(N)
 
  comms = g.community_leading_eigenvector(weights = g.es["weight"])
  clusters = comms
  if verbose:
    print("\nLEV: after calculating the average network the comms are:  ")
    print(clusters)

  c =['deg', 'cc','ne'] * len(clusters)
  dt_lev = pd.DataFrame(columns = c)

  for j , n in enumerate(subject_folder_names):
    Xc = pd.read_csv(path_FC_matrices + filename_of_FC.format(n),header=None)
    Xc = Xc.values
    np.fill_diagonal(Xc,0)
    Xc = Xc.tolist()
    gc = Graph.Weighted_Adjacency(Xc, mode=ADJ_UNDIRECTED,attr="weight", loops=False)
    Nc = gc.vcount()
    Ec = gc.ecount()
    gc.es["label"] = gc.es["weight"]
    lst_c = range(Nc)
    gc.vs["name"] = lst_c

    lists = []

    for i in range(len(clusters)):
      sub = gc.induced_subgraph(clusters[i])
      av_d = average(sub.strength(loops=False,weights=sub.es["weight"]))
      av_cc = average(sub.transitivity_local_undirected(vertices=None,mode="zero", weights=sub.es["weight"]))
      av_ne =average(nodal_eff(sub))

      lists.extend([av_d,av_cc,av_ne])
    dt_lev.loc[j] = lists

  return dt_lev, clusters



