%% Define the path to the main directory and get the subject folder names
% path: the path to the directory that contains the subject folders (e.g.
% subjext_001)

%addpath(genpath('/home/loukas/spm12'));
%path = '/media/loukas-nas3/Data/Serafeim_RS_music/';
path = '/Users/loukas/Desktop/RS_DATA/';
%path_to_functions = '/media/loukas-nas3/Data/Serafeim_RS_music/functions/';
%addpath(path_to_functions);
path_to_functions = '/Users/loukas/Desktop/RS_DATA/functions/';
addpath(path_to_functions);
savepath;

cd(path);
S = dir('*Subject_*');
dirFlags = [S.isdir];
subFolders = S(dirFlags);

%% Create cell with the subject folder names
n_subjects= size(S,1);
subject_name= {};
for k = 1 : n_subjects
	subject_name{k} = subFolders(k).name;
end

%% run the normalization and altasing job
RS_Music_Connectomes_SUBJ_server_NEW(path, subject_name, n_subjects) % with nuisance reg
