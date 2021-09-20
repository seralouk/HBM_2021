"""
=========================================
Author: Serafeim Loukas, Oct 2018
# updated: 2021
=========================================
"""
#%load_ext autoreload
#%autoreload 2
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

import os
import sys
import cairocffi as cairo
from colorama import Fore, Back, Style
import pandas as pd
from scipy import stats
np.random.seed(0)


# Define the paths to utilities, data and results directories.
path = '/Users/loukas/Desktop/RS_DATA/'
path_to_utilities = '/Users/loukas/Desktop/RS_DATA/Pipeline/Utilities/'
path_to_TCS = '/Users/loukas/Desktop/RS_DATA/Pipeline/Dataset/initial_dataset/'

save_results_to = '/Users/loukas/Desktop/RS_DATA/Pipeline/Results/'

sys.path.append(path_to_utilities)
#sys.path.append(path_to_TCS)

from Build_FC import Build_using_Acc_Disc, Build_using_Pearson
from decomposition_functions import *

subject_folder_names = [f for f in sorted(os.listdir(path)) if f.startswith('Subject_')]
n_subjects = len(subject_folder_names)

MP = np.array(['Subject_001',
 'Subject_002',
 'Subject_006',
 'Subject_007',
 'Subject_008',
 'Subject_012',
 'Subject_016',
 'Subject_028',
 'Subject_032',
 'Subject_044',
 'Subject_047',
 'Subject_048',
 'Subject_051',
 'Subject_052',
 'Subject_054'])

FT = np.array(['Subject_003',
 'Subject_013',
 'Subject_014',
 'Subject_015',
 'Subject_020',
 'Subject_024',
 'Subject_025',
 'Subject_033',
 'Subject_034',
 'Subject_037',
 'Subject_038',
 'Subject_039',
 'Subject_040',
 'Subject_041',
 'Subject_045',
 'Subject_056'])


cases = ['PM', 'PN', 'FT']
for _case_ in range(len(cases)):
	subject_folder_names = [f for f in sorted(os.listdir(path)) if f.startswith('Subject_')]
	currect_case = cases[_case_]
	print("Working on {} group".format(currect_case))
	if currect_case == 'PM':
		subject_folder_names = [e for e in subject_folder_names if e in MP] # Perform the analysis only for the MP group
	elif currect_case == 'PN':
		subject_folder_names = [e for e in subject_folder_names if e not in np.concatenate((MP,FT))]
	elif currect_case == 'FT':
		subject_folder_names = [e for e in subject_folder_names if e in FT]
	#print(subject_folder_names)
	case = currect_case

	n_subjects = len(subject_folder_names)

	#########################################################
	############- Build the FC matrices -####################
	#########################################################
	filename_TCS_pre='pre_masked_{}.csv'
	filename_TCS_post='post_masked_{}.csv'

	# ACC / DISC method
	print("Constructing the connectomes using Accordance-Discordance")
	Build_using_Acc_Disc(subject_folder_names, path_to_TCS, filename_TCS_pre, filename_TCS_post, save_results_to)

	# Pearson R method
	print("Constructing the connectomes using Pearson\n")
	Build_using_Pearson(subject_folder_names, path_to_TCS, filename_TCS_pre, filename_TCS_post, save_results_to)

	#########################################################
	###################   - ACCORDANCE      -################
	#################### --------------- ####################

	path_FC = os.path.join(save_results_to, 'FC/', 'Acc/')
	filename_of_FC_pre = 'PRE_acc_{}.csv'
	filename_of_FC_post = 'POST_acc_{}.csv'

	if not os.path.exists(os.path.join(save_results_to, 'Based_on_Acc/')):
		os.makedirs(os.path.join(save_results_to, 'Based_on_Acc/'))

	###################- Global Approach -###################
	print("Extracting Global features - Accordance case")
	pre_acc_global = extract_global(subject_folder_names, path_FC, filename_of_FC_pre)
	pre_acc_global.round(6).to_csv(os.path.join(save_results_to, 'Based_on_Acc/') + 'pre_acc_global.csv')

	post_acc_global = extract_global(subject_folder_names, path_FC, filename_of_FC_post)
	post_acc_global.round(6).to_csv(os.path.join(save_results_to, 'Based_on_Acc/') + 'post_acc_global.csv')


	###################- Nodal Approach -###################
	print("Extracting Nodal features - Accordance case")
	pre_acc_nodal = extract_nodal(subject_folder_names, path_FC
		, filename_of_FC_pre)
	pre_acc_nodal.round(6).to_csv(os.path.join(save_results_to, 'Based_on_Acc/') + 'pre_acc_nodal.csv')

	post_acc_nodal = extract_nodal(subject_folder_names, path_FC, filename_of_FC_post)
	post_acc_nodal.round(6).to_csv(os.path.join(save_results_to, 'Based_on_Acc/') + 'post_acc_nodal.csv')


	###################- Connections Approach -###################
	number_of_ROIs = 90
	print("Extracting Connection features - Accordance case")
	pre_acc_connections = extract_connections(subject_folder_names, path_FC, filename_of_FC_pre, number_of_ROIs)
	# write data to csv
	np.savetxt(os.path.join(save_results_to, 'Based_on_Acc/') + 'pre_acc_connections.csv', np.round(pre_acc_connections,6), delimiter=",")

	post_acc_connections = extract_connections(subject_folder_names, path_FC, filename_of_FC_post, number_of_ROIs)
	np.savetxt(os.path.join(save_results_to, 'Based_on_Acc/') + 'post_acc_connections.csv', np.round(post_acc_connections,6), delimiter=",")


	###################- Module Approach -###################
	## Fast greedy
	print("Extracting Module features with fast greedy- Accordance case")
	Average_network_pre = Build_mean_network(subject_folder_names, path_FC, filename_of_FC_pre, number_of_ROIs)
	pre_acc_module_fg, clusters_pre = decompose_fg(Average_network_pre, subject_folder_names, path_FC,filename_of_FC_pre, verbose=0)
	pre_acc_module_fg.round(6).to_csv(os.path.join(save_results_to, 'Based_on_Acc/') + 'pre_acc_module_fg.csv')

	post_acc_module_fg, clusters_post = decompose_fg(Average_network_pre, subject_folder_names, path_FC,filename_of_FC_post, verbose=0)
	post_acc_module_fg.round(6).to_csv(os.path.join(save_results_to, 'Based_on_Acc/') + 'post_acc_module_fg.csv')

	## Leading Eigenvector
	print("Extracting Module features with leading eigenvector- Accordance case\n")
	pre_acc_module_lev, clusters_pre = decompose_lev(Average_network_pre, subject_folder_names, path_FC,filename_of_FC_pre, verbose=0)
	pre_acc_module_lev.round(6).to_csv(os.path.join(save_results_to, 'Based_on_Acc/') + 'pre_acc_module_lev.csv')

	post_acc_module_lev, clusters_post = decompose_lev(Average_network_pre, subject_folder_names, path_FC,filename_of_FC_post, verbose=0)
	post_acc_module_lev.round(6).to_csv(os.path.join(save_results_to, 'Based_on_Acc/') + 'post_acc_module_lev.csv')

	del(path_FC,filename_of_FC_pre, filename_of_FC_post, Average_network_pre)
	print("All features were extracted and saved to {} using Accordance\n".format(os.path.join(save_results_to, 'Based_on_Acc/')))


	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	############-   STATISTICAL TESTS    -###################
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	print("Starting the statistical analysis")

	# function to calculate Cohen's d for independent samples
	def cohend(d1, d2):
		# calculate the size of samples
		n1, n2 = len(d1), len(d2)
		# calculate the variance of the samples
		s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
		# calculate the pooled standard deviation
		s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
		# calculate the means of the samples
		u1, u2 = np.mean(d1), np.mean(d2)
		# calculate the effect size
		return (u1 - u2) / s

	#Project back the flattened upper triangle to (i,j) matrix indices
	# for flattened UT without the diagonal
	def flattened2ij(k, n):
	    rows = 0
	    for t, cols in enumerate(range(n - 1, -1, -1)):
	        rows += cols
	        if k in range(rows):
	            return (t, n - (rows - k))
	    return None

	# load the atlas regional names
	AAL_names = pd.read_csv('/Users/loukas/Desktop/RS_DATA/Pipeline/AAL.csv', header = None, index_col=None)

	# where to store the paired t-tests results and p values
	store_to = os.path.join(save_results_to, 'T_tests/')
	if not os.path.exists(store_to):
		os.makedirs(store_to)

	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	###############   - Based on ACCORDANCE      -###########
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	print("\nACCORDANCE CASE\n")
	print("Removing connections that are 0 for more than 8 subjects - ACCORDANCE Case")
	# remove columns that contain a lot of 0s
	remove_0s_connections = 1
	remove_pre = []
	remove_post = []
	n = (number_of_ROIs * (number_of_ROIs-1)) / 2.0
	if remove_0s_connections:
		for column_id in range(int(n)):
			if sum(pre_acc_connections[:,column_id] == 0) > 8:
				remove_pre.append(column_id)
			if sum(post_acc_connections[:,column_id] == 0) > 8:
				remove_post.append(column_id)
	merged_idx = remove_post + remove_pre
	merged_idx = np.sort(merged_idx)
	merged_idx = np.unique(merged_idx)

	if sum(merged_idx)!=0:
		pre_acc_connections = np.delete(pre_acc_connections, merged_idx, axis = 1)
		post_acc_connections = np.delete(post_acc_connections, merged_idx, axis = 1)


	print("Performing the paired t-tests - ACCORDANCE Case")
	# perform the paired t-tests on global, nodal, connection, module extracted network features
	t_acc_global, p_acc_global = stats.ttest_rel(pre_acc_global.values, post_acc_global.values)
	t_acc_nodal, p_acc_nodal = stats.ttest_rel(pre_acc_nodal.values, post_acc_nodal.values)
	t_acc_connections, p_acc_connections = stats.ttest_rel(pre_acc_connections, post_acc_connections)
	t_acc_module_fg, p_acc_module_fg = stats.ttest_rel(pre_acc_module_fg.values, post_acc_module_fg.values)
	t_acc_module_lev, p_acc_module_lev = stats.ttest_rel(pre_acc_module_lev.values, post_acc_module_lev.values)

	np.savetxt(store_to + 'p_acc_global.csv', np.concatenate((t_acc_global.reshape(1,-1).round(6) ,p_acc_global.reshape(1,-1).round(6)),axis=0), delimiter=',',fmt='%.6f')
	np.savetxt(store_to + 'p_acc_nodal.csv', np.concatenate((t_acc_nodal.reshape(1,-1).round(6) ,p_acc_nodal.reshape(1,-1).round(6)),axis=0), delimiter=',',fmt='%.6f')
	np.savetxt(store_to + 'p_acc_connections.csv', np.concatenate((t_acc_connections.reshape(1,-1).round(6) ,p_acc_connections.reshape(1,-1).round(6)),axis=0), delimiter=',',fmt='%.6f')
	np.savetxt(store_to + 'p_acc_module_fg.csv', np.concatenate((t_acc_module_fg.reshape(1,-1).round(6) ,p_acc_module_fg.reshape(1,-1).round(6)),axis=0), delimiter=',',fmt='%.6f')
	np.savetxt(store_to + 'p_acc_module_lev.csv', np.concatenate((t_acc_module_lev.reshape(1,-1).round(6) ,p_acc_module_lev.reshape(1,-1).round(6)),axis=0), delimiter=',',fmt='%.6f')
	print("ACCORDANCE Case: The paired t-test results were saved to {}".format(store_to))

	print(sum(p_acc_global < 0.05))
	print(sum(p_acc_nodal < 0.05))
	print(sum(p_acc_connections < 0.05))
	print(sum(p_acc_module_fg < 0.05))
	print(sum(p_acc_module_lev < 0.05))
	print('\n')

	print("Multiple comparison FDR correction plots")
	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
	#FDR multiple comparison correction for NODAL
	q = 0.05
	N = float(p_acc_connections.shape[0])
	sorted_p_acc_connections = np.sort(p_acc_connections)
	sorted_p_acc_connections = sorted_p_acc_connections / 2 # one-sided
	i = np.arange(1,N+1)
	plt.plot(i, sorted_p_acc_connections, 'b.', label = 'p(i)')
	plt.plot(i, q * i / N, 'r', label = 'q i / N')
	plt.xlabel('i (index)')
	plt.ylabel('p-value')
	plt.legend()
	plt.title('FDR plot for {} group'.format(case))
	plt.savefig('/Users/loukas/Desktop/' + 'FDR_plot_{}'.format(case))
	#plt.show()
	below = sorted_p_acc_connections <= (q * i / N)
	if sum(below) != 0:
		max_below = np.max(np.where(below)[0])
		print('p value:',sorted_p_acc_connections[max_below])
		print('i:', max_below + 1)
	else:
		print("No results survived")

	# verify also with statsmodels module
	from statsmodels.stats.multitest import fdrcorrection
	rejected, stat_p_adj = fdrcorrection(sorted_p_acc_connections, alpha=0.1, method='indep', is_sorted=True)
	print("statsmodels toolbox: n={} sign. results \n".format(sum(rejected)))
	print(stat_p_adj[np.where(rejected)[0]])

	print("Accordance Case: Find significant results based on the paired t-tests for CONNECTIONS features")
	####################################
	## Find the differences: CONNECTIONS
	connections_column_names = np.array(range(4005))
	if sum(merged_idx)!=0:
		connections_column_names = np.delete(connections_column_names, merged_idx)

	# find = ((p_acc_connections<0.001) & (t_acc_connections > 0)) # pre>post
	# connections_table= pd.DataFrame(columns = ['region1','region2','p_value','Cohen_d'])
	# connections_table.region1 = connections_table.region1.astype(object)
	# connections_table.Cohen_d = connections_table.Cohen_d.astype(object)
	# connections_table.region2 = connections_table.region2.astype(object)
	# connections_table.p_value = connections_table.p_value.astype(object)
	# for_brain_net_accordance_pre_post = []
	# for idx, con in enumerate(connections_column_names[find]):
	# 	tmp = flattened2ij(con,90)
	# 	co = cohend(pre_acc_connections[:,find][:,idx], post_acc_connections[:,find][:,idx])
	# 	reg1 = str(AAL_names.values[tmp[0]][1])
	# 	reg2 = str(AAL_names.values[tmp[1]][1])
	# 	p_val = p_acc_connections[find][idx]
	# 	for_brain_net_accordance_pre_post.append(tmp)
	# 	for_brain_net_accordance_pre_post.append(p_val)
	# 	connections_table.set_value(idx, 'region1', reg1)
	# 	connections_table.set_value(idx, 'region2', reg2)
	# 	connections_table.set_value(idx, 'p_value', p_val)
	# 	connections_table.set_value(idx, 'Cohen_d', co)
	# connections_table = connections_table.sort_values('p_value')
	# connections_table.to_csv(save_results_to + "Acc_connections_table_pre>post.csv")


	find = ((p_acc_connections<0.001) & (t_acc_connections < 0)) # pre<post
	connections_table= pd.DataFrame(columns = ['region1', 'region2', 'p_value', 'Cohen_d'])
	connections_table.region1 = connections_table.region1.astype(object)
	connections_table.Cohen_d = connections_table.Cohen_d.astype(object)
	connections_table.region2 = connections_table.region2.astype(object)
	connections_table.p_value = connections_table.p_value.astype(object)
	for_brain_net_accordance_post_pre = []
	for idx, con in enumerate(connections_column_names[find]):
		tmp = flattened2ij(con,90)
		co = cohend(pre_acc_connections[:,find][:,idx], post_acc_connections[:,find][:,idx])
		reg1 = str(AAL_names.values[tmp[0]][1])
		reg2 = str(AAL_names.values[tmp[1]][1])
		p_val = p_acc_connections[find][idx]
		for_brain_net_accordance_post_pre.append(tmp)
		for_brain_net_accordance_post_pre.append(p_val)
		#connections_table.set_value(idx, 'region1', reg1)
		connections_table.loc[idx,'region1'] = reg1
		#connections_table.set_value(idx, 'region2', reg2)
		connections_table.loc[idx,'region2'] = reg2
		#connections_table.set_value(idx, 'p_value', p_val)
		connections_table.loc[idx,'p_value'] = p_val
		#connections_table.set_value(idx, 'Cohen_d', co)	
		connections_table.loc[idx,'Cohen_d'] = co
	connections_table = connections_table.sort_values('p_value')
	connections_table.to_csv(save_results_to + "{}_Acc_connections_table_pre<post.csv".format(case))

	#print(connections_table)
