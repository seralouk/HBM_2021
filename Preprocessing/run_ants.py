"""Step 4: Bring the atlas to subject space. This works for each subject separately 
and we do not normalize/ sample the functional data which is good.

This script allows the user to run the ANTs and bring the atlas to the subject space before performing atlasing of the functional data.
# v1.0 Serafeim Loukas, 22 Oct 2018

I M P O R T A N T: This is computationally very expensive. Run it on a Server. ETA: 3H / subject

"""
import os
import nibabel as nib
import numpy as np
import ants

# Paths on server
# Path to the folder containing the Subject folders
path = '/media/miplab-nas3/Data/Serafeim_RS_music/'

# Path to save the transformed atlases
save_global_data_to = '/home/loukas/rs_music_test/Results/atlases/'

# Path to save the transformed TPMs
save_registered_TPMs_to = '/home/loukas/rs_music_test/Results/registered_atlas/'

# Path to the UNC atlas (to be used as input in this algorithm)
path_to_global_mask = '/media/miplab-nas3/Data/Serafeim_dHCP_data/DHCP/derivatives/UNC_atlas/'


# Local paths in case I run them locally
## path = '/Users/miplab/Desktop/RS_DATA/'
## save_global_data_to = '/Users/miplab/Desktop/RS_DATA/RESULTS/atlases/'
## save_registered_TPMs_to = '/Users/miplab/Desktop/RS_DATA/RESULTS/registered_atlas/'
## path_to_global_mask = '/Users/miplab/Desktop/RS_DATA/UNC_atlas/'

# List of Subjects that belong to MP group.
MP= list(['Subject_001',
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

# Get the subject names
subject_folder_names = [f for f in sorted(os.listdir(path)) if f.startswith('Subject_')]

# Remove the subjects that have been already preprocessed (it does not always apply, so comment it out)
subject_folder_names = [e for e in subject_folder_names if e not in MP]
n_subjects = len(subject_folder_names)

# adding the new 3 subjects
# new_names = list()
# for file in subject_folder_names:
# 	if file.startswith(('Subject_006', 'Subject_028','Subject_044')):
# 		new_names.append(file)
# subject_folder_names = new_names[:]
# n_subjects = len(subject_folder_names)

# Create the folders where you wish to save the results
if not os.path.exists(save_global_data_to):
    os.mkdir(save_global_data_to)
if not os.path.exists(save_registered_TPMs_to):
    os.mkdir(save_registered_TPMs_to)


def Build_Connectomes(subject_folder_names):
	"""This is the main function

	Parameters
	----------
	subject_folder_names : list
		A list containing the subject folders names (it is automatically created in line 53)
	
	Returns
	----------
	empty: results are automatically saved to the defined paths / folders

	"""
	n_subjects = len(subject_folder_names)
	for i in range(n_subjects):
		subject_folder_names = subject_folder_names
		
		# Reading files and defining paths
		go_back = os.getcwd()
		subj_dir=os.path.join(path, subject_folder_names[i])
		session_name = [f for f in os.listdir(os.path.join(path, subject_folder_names[i])) if f.startswith('ses-')]
		session_dir = os.path.join(path, subject_folder_names[i],session_name[0])

		# Get the anatomical image of the subject 'i'
		anat_dir_name = [f for f in os.listdir(session_dir) if f.startswith('anat')]
		T2_struct_image = [f for f in os.listdir(os.path.join(session_dir, anat_dir_name[0])) if (f.startswith('s') & f.endswith('.nii'))]
		path_to_T2 = os.path.join(session_dir, anat_dir_name[0] , T2_struct_image[0])

		# Get the atlas image and the TPMs in atlas space
		global_mask_name = [f for f in os.listdir(path_to_global_mask) if f.startswith('infant-neo.nii')]
		TMP_GM_name = [f for f in os.listdir(path_to_global_mask) if f.startswith('infant-neo-seg-gm')]
		TMP_WM_name = [f for f in os.listdir(path_to_global_mask) if f.startswith('infant-neo-seg-wm')]
		TMP_CSF_name = [f for f in os.listdir(path_to_global_mask) if f.startswith('infant-neo-seg-csf')]

		path_global_mask = os.path.join(path_to_global_mask , global_mask_name[0])
		path_TMP_GM_name = os.path.join(path_to_global_mask , TMP_GM_name[0])
		path_TMP_WM_name = os.path.join(path_to_global_mask , TMP_WM_name[0])
		path_TMP_CSF_name = os.path.join(path_to_global_mask , TMP_CSF_name[0])


		# ANTs functions to get the deformation fields and then apply them on the atlas to bring it to subject space
		# fixed image to which we register the moving image
		# os.chdir(path_to_T2)
		fixed = ants.image_read(path_to_T2)
		if not fixed:
			raise ValueError('T2 was not found')
		else:
			print('\nT2 was selected as {}'.format(T2_struct_image))

		# Moving image to be mapped to fixed space.
		moving = ants.image_read(path_global_mask)
		if not moving:
			raise ValueError('Atlas was not found')
		else:
			print('\nAtlas was selected as {}'.format(global_mask_name))

		print("Registration of the atlas to the T2 subject space")
		mytx = ants.registration(fixed=fixed , moving=moving, type_of_transform='SyNCC', syn_sampling=8)

		# Some information:
		# dict containing follow key/value pairs:
		#         `warpedmovout`: Moving image warped to space of fixed image.
		#         `warpedfixout`: Fixed image warped to space of moving image.
		#         `fwdtransforms`: Transforms to move from moving to fixed image.
		#         `invtransforms`: Transforms to move from fixed to moving image

		# warped_moving = mytx['warpedmovout']
		# fixed.plot(overlay=warped_moving,title='After Registration')

		# You can also use the transforms output from registration and apply them directly to the image:
		moving_atlas_name = [f for f in os.listdir(path_to_global_mask) if f.startswith('infant-neo-aal')]
		path_moving_atlas = os.path.join(path_to_global_mask , moving_atlas_name[0])

		moving_atlas = ants.image_read(path_moving_atlas)
		if not moving_atlas:
			raise ValueError('Label image was not found')
		else:
			print('\nLabel image was selected as {}'.format(moving_atlas_name))

		print("Applying the deformation to the Label image")
		mywarpedimage = ants.apply_transforms(fixed=fixed, moving=moving_atlas,transformlist=mytx['fwdtransforms'], interpolator='genericLabel')

		print("Applying the deformation to the TPM images")
		C1 = ants.image_read(path_TMP_GM_name)
		C2 = ants.image_read(path_TMP_WM_name)
		C3 = ants.image_read(path_TMP_CSF_name)
		c1warped = ants.apply_transforms(fixed=fixed, moving=C1,transformlist=mytx['fwdtransforms'], interpolator='bSpline')
		c2warped = ants.apply_transforms(fixed=fixed, moving=C2,transformlist=mytx['fwdtransforms'], interpolator='bSpline')
		c3warped = ants.apply_transforms(fixed=fixed, moving=C3,transformlist=mytx['fwdtransforms'], interpolator='bSpline')

		# Some information about ants.apply_transforms function
			## ants.apply_transforms function:
		    # Arguments
		    # ---------
		    # fixed : ANTsImage
		    #     fixed image defining domain into which the moving image is transformed.

		    # moving : AntsImage
		    #     moving image to be mapped to fixed space.

		    # transformlist : list of strings
		    #     list of transforms generated by ants.registration where each transform is a filename.

		os.chdir(save_global_data_to)
		mywarpedimage.to_file('atlas_ML_{}.nii'.format(subject_folder_names[i]))
		os.chdir(save_registered_TPMs_to)
		c1warped.to_file('C1_{}.nii'.format(subject_folder_names[i]))
		c2warped.to_file('C2_{}.nii'.format(subject_folder_names[i]))
		c3warped.to_file('C3_{}.nii'.format(subject_folder_names[i]))
		print("\nRegistered atlas was saved in {}".format(save_global_data_to))
		print("\nRegistered TPMs were saved in {}".format(save_registered_TPMs_to))

		print ("Subject {} is done !".format(subject_folder_names[i]))
		print ("Subject number {} out of 40 is done !".format(i+1))

	return

###########################################################################################
###########################################################################################

# Run it
# subject_folder_names is defined in line 53
Build_Connectomes(subject_folder_names)

