"""Step 2: Keep only imaging data that have a FD below 0.5mm (data cleaning step 2)

This script allows the user to keep only the imaging data that have a FD below 0.5mm. 
This information is stored in a matlab matrix format and should be inside the functional folders (fMRI runs folders)
# v1.0 Serafeim Loukas, 22 Oct 2018

I M P O R T A N T: Be sure that you know what are you doing. This code will delete files with no recovery option available.

"""
import os
from subprocess import call
import scipy.io

# The path to the location of the folder containg the subjects e.g. path = '/Users/Analysis/'. Do not forget the backslash
path = '/Users/loukas/Desktop/RS_DATA/'

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


# Iterate within each subject's folder, inside the 'func', 'anat', etc subfolders and keep only the 'GOOD' volumes
for i in range(n_subjects):
	# Reading files and defining paths
	go_back = os.getcwd()
	subj_dir=os.path.join(path, subject_folder_names[i])
	session_name = [f for f in os.listdir(os.path.join(path, subject_folder_names[i])) if f.startswith('ses-')]
	session_dir = os.path.join(path, subject_folder_names[i],session_name[0])

	func_dir_name = [f for f in sorted(os.listdir(session_dir)) if f.startswith('func_')]
	#anat_dir_name = [f for f in os.listdir(session_dir) if f.startswith('anat')]

	# Iterate to all functional runs (if you have more than 1 fMRI run)
	for jj in range(len(func_dir_name)):
		run = os.path.join(path, subject_folder_names[i],session_name[0],func_dir_name[jj] )
		print("deleting files in {}".format(run))

		# Read the .mat to get the indices of the good images (e.g. of indices, good = ['1','2'], list type)
		mat = scipy.io.loadmat(run + '/' + "good_volumes_FD0.5_DVARS30.mat")
		good_files =  mat['good']

		for general_file in os.listdir(run):
			try:
				file_numbers = int(general_file.split('-')[-2])
				if file_numbers not in good_files:
					call(['rm', run + '/' + general_file], shell= False)
			except:
				continue

	print( subject_folder_names[i] + ' ' +'Done !')




