%%img2nii.m--------------------------------------------

function Convert_img_to_nii()
%% Script to convert hdr/img files to 4D nii.
% v1.0 Serafeim Loukas, 22 Oct 2018
path = '/Users/miplab/Desktop/RS_DATA/';

cd(path);
S = dir('*Subject_*');
dirFlags = [S.isdir];
subFolders = S(dirFlags);
prefix_of_files_to_be_selected ='arf';
prefix_anat = 's';

%% Create cell with the subject folder names
n_subjects= size(S,1);
subject_name= {};
for k = 1 : n_subjects
	subject_name{k} = subFolders(k).name;
end

%% the next 3 lines of code are only for processing the additional 3 new
% subjects
%filtering = {'Subject_47'};
% filtering = {'Subject_06','Subject_28','Subject_44'};
%subject_name(~ismember(subject_name,filtering)) = [];
%n_subjects=size(subject_name,2);

for i =1:n_subjects
    fprintf('Working on subject %s\n',subject_name{i});
    go_back = pwd;
    path_to_subject = [path subject_name{i} filesep];
    cd(path_to_subject)
    Session_name = dir('*ses-*');
    if isempty(Session_name.name)
      error('Session name could not be defined')
    else
      disp('Session name was found !')
    end

    path_to_sess = [path subject_name{i} filesep Session_name.name];
    cd(path_to_sess);
    run_names = dir('*func_*');
    run_names = {run_names.name};
    if size(run_names,2) == 0
      error('Run names could not be defined')
    else
      disp('Run names were found !')
    end
    
    for ii = 1:size(run_names,2)
        fprintf('Working on run %s for subject %s\n',run_names{ii},subject_name{i} );
        step=0;
        path_to_run = [path subject_name{i} filesep Session_name.name filesep run_names{ii} ];
        %cd(path_to_run);
        
        %select files
        f = spm_select('FPList' , path_to_run  , ['^' prefix_of_files_to_be_selected '.*\.img$']);
        %f = spm_select(Inf,'img$','Select img files to be converted');
        fprintf('%d images with the specified prefix were found\n',size(f,1));

        %convert img files to nii
        for j=1:size(f,1)
          input = deblank(f(j,:));
          [pathstr,fname,ext] = fileparts(input);
           output = [pathstr filesep fname,'.nii'];
           V=spm_vol(input);
           ima=spm_read_vols(V);
           V.fname=output;
           spm_write_vol(V,ima);
        end
        
        %convert functional nii files to 4D nii file
        f = spm_select('FPList' , path_to_run  , ['^' prefix_of_files_to_be_selected '.*\.nii$']);
        %fprintf('%d images with the specified prefix were found\n',size(f,1));
        matlabbatch{step+1}.spm.util.cat.vols = cellstr(f);
        matlabbatch{step+1}.spm.util.cat.name = ['4D_' subject_name{i} '_' run_names{ii} '.nii'];
        matlabbatch{step+1}.spm.util.cat.dtype = 0;
        
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch );   
        clearvars matlabbatch;
        
    end
    
    path_to_anat = [path subject_name{i} filesep Session_name.name filesep 'anat' ];
    f_anat= spm_select('FPList' , path_to_anat  , ['^' prefix_anat '.*\.img$']);
    if size(f_anat,1) ~= 1
        error("Less than 1 T2 images were found")
    end

    %convert img files to nii
    input = deblank(f_anat);
    [pathstr,fname,ext] = fileparts(input);
     output = [pathstr filesep 's_' subject_name{i},'.nii'];
     V=spm_vol(input);
     ima=spm_read_vols(V);
     V.fname=output;
     spm_write_vol(V,ima);
end
    
    
return