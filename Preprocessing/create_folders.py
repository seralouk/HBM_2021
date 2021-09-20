"""
=========================================
Author: Serafeim Loukas, August 2018
=========================================

"""
import os

def create_directories(location, root):

	try: 
		location == [] or root == []
	except ValueError:
		print("Oops!  Specify the location of the root folder..Try again...")

	main_dir=[]
	for i in range(len(root)):
		main_dir.append(root[i]) 

	common_dir1 = ["ses-1"]
	common_dir2 = ["anat", "func_01", "func_02"]

	for dir1 in main_dir:
		for dir2 in common_dir1:
			for dir3 in common_dir2:
				try: os.makedirs(location+ '/' + os.path.join(dir1,dir2,dir3))
				except OSError: pass


if __name__ == "__main__":
	root = ['Subject_004', 'Subject_005', 'Subject_009','Subject_010','Subject_011','Subject_022', 'Subject_035', 'Subject_036','Subject_049','Subject_050','Subject_053','Subject_055','Subject_060','Subject_061','Subject_062']
	create_directories(os.getcwd(), root)