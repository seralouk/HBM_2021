"""Step 1: Remove all imaging data except the realigned images (data cleaning step 1)

This script allows the user to remove all the imaging data except images of a specific type
# v1.0 Serafeim Loukas, 22 Oct 2018

I M P O R T A N T: Be sure that you know what are you doing. This code will delete files with no recovery option available.

"""
import os
from subprocess import call
import scipy.io

def remove_all_expept_arf(path):
	"""This is the main function
	
	Parameters
	----------
	path : str
		The path to the location of the folder containg the subjects e.g. path = '/Users/Analysis/'. Do not forget the backslash
	"""
	path = path[:]
	# Get the subject names
	subject_folder_names = [f for f in sorted(os.listdir(path)) if f.startswith('Subject_')]
	n_subjects = len(subject_folder_names)
	print("found {} subjects in total".format(n_subjects))

	# new_names = list()
	# for haha in subject_folder_names:
	# 	if haha.startswith(('Subject_06', 'Subject_28','Subject_44')):
	# 		new_names.append(haha)

	# subject_folder_names = new_names[:]
	# n_subjects = len(subject_folder_names)

	# Iterate within each subject's folder, inside the 'func', 'anat', etc subfolders and delete files
	for i in range(n_subjects):
		subject_folder_names = subject_folder_names
		# Reading files and defining paths
		go_back = os.getcwd()
		subj_dir=os.path.join(path, subject_folder_names[i])
		session_name = [f for f in os.listdir(os.path.join(path, subject_folder_names[i])) if f.startswith('ses-')]
		session_dir = os.path.join(path, subject_folder_names[i],session_name[0])
		func_dir_name = [f for f in os.listdir(session_dir) if f.startswith('func_')]
		anat_dir_name = [f for f in os.listdir(session_dir) if f.startswith('anat')]

		# Iterate to all functional runs (if you have more than 1 fMRI run)
		for jj in range(len(func_dir_name)):
			run = os.path.join(path, subject_folder_names[i],session_name[0],func_dir_name[jj] )
			print("deleting files in {}".format(run))

			# I M P O R T A N T : remove all files except anything that starts with 'arf' or 'good_volumes_'
			remove_name = [f for f in os.listdir(run) if not f.startswith(('arf','good_volumes_'))]
			for file in remove_name:
				call(['rm', '-rf', run + '/' + file], shell= False)

		# Iterate to the anatomical folder and remove unwanted file types (everything that does not start with 's')
		anat_run = os.path.join(path, subject_folder_names[i],session_name[0],anat_dir_name[0] )
		# I M P O R T A N T : remove all files except anything that starts with 's'
		remove_anat_names = [f for f in os.listdir(anat_run) if not f.startswith('s')]

		for files in remove_anat_names:
			call(['rm', '-rf', anat_run + '/' + files], shell= False)
		print( subject_folder_names[i] + ' ' +'Done !')

	return


## Run the function
path = '/Users/miplab/Desktop/RS_DATA/'
remove_all_expept_arf(path)


