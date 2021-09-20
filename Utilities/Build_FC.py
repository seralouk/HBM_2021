"""
=========================================
Author: Serafeim Loukas, Oct 2018
=========================================
"""
# Get rid of warnings
def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

import pandas as pd
import numpy as np
import os

#########################################################################################################
#########################################################################################################

def standardize_TCS(X):
	mean = X.mean(axis =1)
	std = X.std(axis= 1, ddof=1)
	X= X - mean[:, np.newaxis]
	X= X / std[:, np.newaxis]
	return X

#########################################################################################################
#########################################################################################################

def acc_disc(X,u,l): #define main function
	X=np.array(X) #type of X array
	N=X.shape[0] #N:regions (in X,the rows represent the regions)
	T=X.shape[1] #T:time (in X,the columns represent the time)

	Xu = np.zeros(X.shape) # create zero arrays
	Xl = np.zeros(X.shape)

	Xu[X >= u] = 1 # apply thresholds on X and put 1 in Xu 
	Xl[X <= l] = -1 #for the corresponding positions 

	Txu=np.transpose(Xu) #transpose Xu and Xl
	Txl=np.transpose(Xl) #we use them in the next lines

	Scav =(float(1) / float(T))*(np.dot(Xu,Txu))
	Scdv =(float(1) / float(T))*(np.dot(Xl,Txl))
	Sa = Scav + Scdv
	E =np.diag(Sa) ** (- 0.5) #Calculate energy
	d_E=np.diag(E)
	Acc =d_E.dot(Sa).dot(d_E) #Calculate Accordance array

	Sd =(float(1) /float(T))*(np.dot(Xu,Txl) + np.dot(Xl,Txu))
	Dis=d_E.dot(Sd).dot(d_E) #Calculate Discordance array

	return Acc, Dis

#########################################################################################################
#########################################################################################################

def Build_using_Acc_Disc(subject_folder_names, path_to_TCS, filename_TCS_pre, filename_TCS_post, save_results_to):

	if not os.path.exists(os.path.join(save_results_to, 'FC/', 'Acc/')):
		os.makedirs(os.path.join(save_results_to, 'FC/', 'Acc/'))
	if not os.path.exists(os.path.join(save_results_to, 'FC/', 'Disc/')):
		os.makedirs(os.path.join(save_results_to, 'FC/', 'Disc/'))

	# build pre matrices
	for s in subject_folder_names:
		X=pd.read_csv(path_to_TCS + filename_TCS_pre.format(s),header=None)
		X = X.T
		#print(X.shape)
		X=X.values
		X = standardize_TCS(X)
		Acc, Dis = acc_disc(X, 0.84, -0.84)
		np.fill_diagonal(Acc,0.0)
		Dis = 1.0-np.abs(Dis)
		np.fill_diagonal(Dis,0.0)

		np.savetxt(os.path.join(save_results_to, 'FC/', 'Acc/') + 'PRE_acc_{}.csv'.format(s), Acc, delimiter=',')
		np.savetxt(os.path.join(save_results_to, 'FC/', 'Disc/') + 'PRE_disc_{}.csv'.format(s), Dis, delimiter=',')

	del(Acc,Dis,X,s)

	# build post matrices
	for s in subject_folder_names:
		X=pd.read_csv(path_to_TCS + filename_TCS_post.format(s),header=None)
		X = X.T
		X=X.values
		X = standardize_TCS(X)
		Acc, Dis = acc_disc(X, 0.84, -0.84)
		np.fill_diagonal(Acc,0.0)
		Dis = 1.0-np.abs(Dis)
		np.fill_diagonal(Dis,0.0)

		np.savetxt(os.path.join(save_results_to, 'FC/', 'Acc/') + 'POST_acc_{}.csv'.format(s), Acc, delimiter=',')
		np.savetxt(os.path.join(save_results_to, 'FC/', 'Disc/') + 'POST_disc_{}.csv'.format(s), Dis, delimiter=',')
	del(X,s)
	
	return

#########################################################################################################
#########################################################################################################


def Build_using_Pearson(subject_folder_names, path_to_TCS, filename_TCS_pre, filename_TCS_post, save_results_to):


	if not os.path.exists(os.path.join(save_results_to, 'FC/', 'Pearson/')):
		os.makedirs(os.path.join(save_results_to, 'FC/', 'Pearson/'))

	# build pre matrices
	for s in subject_folder_names:
		X=pd.read_csv(path_to_TCS + filename_TCS_pre.format(s),header=None)
		N_vol_pre = X.shape[0]
		X = X.T
		X=X.values
		X = standardize_TCS(X)
		X = np.corrcoef(X, rowvar=True)
		np.fill_diagonal(X,0)
		threshld_pre = (1.0 / np.sqrt(N_vol_pre))
		X[X <= threshld_pre] = 0.0

		np.savetxt(os.path.join(save_results_to, 'FC/', 'Pearson/') + 'PRE_pearson_{}.csv'.format(s), X, delimiter=',')

	del(X,s)

	# build post matrices
	for s in subject_folder_names:
		X=pd.read_csv(path_to_TCS + filename_TCS_post.format(s),header=None)
		N_vol_post = X.shape[0]
		X = X.T
		X=X.values
		X = standardize_TCS(X)
		X = np.corrcoef(X, rowvar=True)
		np.fill_diagonal(X,0)
		threshld_post = (1.0 / np.sqrt(N_vol_post))
		X[X <= threshld_post] = 0.0

		np.savetxt(os.path.join(save_results_to, 'FC/', 'Pearson/') + 'POST_pearson_{}.csv'.format(s), X, delimiter=',')

	return

#########################################################################################################
#########################################################################################################










