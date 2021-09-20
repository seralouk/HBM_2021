function RS_Music_Connectomes_SUBJ_server_NEW(path, subject_name, n_subjects  )
%%%
% This function uses the UNC atlas to extract ROI TCS. It does segmentation
% of the T2, atlas masking with GM mask in the native space (T2) so that
% only GM voxels are considered, but also global whole brain masking is
% applied. 2 connectomes are outputed for each subject and each run.
% Cleaning (detrending, regression of motion, CSF, WM signals) is done and 
% the clean TCS are put back into brain images. The clean data are used to extract
% atlas-based TCS (90 regions) after warping the atlas into the functional
% native space.
%
% path: the path to the directory that contains the subject folders (e.g.
% subjext_001)
%
% v1.0 Oct 2018 Serafeim Loukas
% v2.0 Jul 2019 Revised by Serafeim Loukas (added data cleaning part)
% Compatible with SPM12 and SPM8. Defaults are based on SPM8


%save_global_connectomes_to = '/home/loukas/rs_music_test/Results/CONNECTOMES/';
%save_masked_connectomes_to = '/home/loukas/rs_music_test/Results/CONNECTOMES/';
%path_to_registered_global_atlas = '/home/loukas/rs_music_test/Results/atlases/';
%save_registered_TPMs_to = '/home/loukas/rs_music_test/Results/registered_atlas/';

%path_to_global_mask = '/media/miplab-nas3/Data/Serafeim_dHCP_data/DHCP/derivatives/UNC_atlas/';
%path_to_registered_TPMs= '/home/loukas/rs_music_test/Results/registered_atlas/';

save_global_connectomes_to = '/Users/loukas/Desktop/RS_DATA/RESULTS/CONNECTOMES/may_2021/';
save_masked_connectomes_to = '/Users/loukas/Desktop/RS_DATA/RESULTS/CONNECTOMES/may_2021/';
path_to_registered_global_atlas = '/Users/loukas/Desktop/RS_DATA/RESULTS/atlases/';
path_to_registered_TPMs= '/Users/loukas/Desktop/RS_DATA/RESULTS/registered_atlas/';


%%
for i =1:n_subjects
        step = 0;
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
        
    
    
        %% load the T2 image to do the segmentation
        path_to_T2 = [path  subject_name{i} filesep Session_name.name filesep 'anat' filesep];
        filepartition =regexp(path_to_T2,filesep,'split');
        T2_img = spm_select('FPList' , path_to_T2  , ['^' 's_' filepartition{end-3} '.*\.nii$']);
        if isempty(T2_img)
            error('No T2 was found')
        else
            disp(['The T2 was selected:' T2_img] )
        end
        
        %% Create SPM 8 segmentation batch

        %[filepath,name,ext] = fileparts(T2_img);
        % image to-be-segmented (the structural)
        matlabbatch{step+1}.spm.tools.oldseg.data = {T2_img};
        matlabbatch{step+1}.spm.tools.oldseg.output.GM = [0 0 1];
        matlabbatch{step+1}.spm.tools.oldseg.output.WM = [0 0 1];
        matlabbatch{step+1}.spm.tools.oldseg.output.CSF = [0 0 1];
        matlabbatch{step+1}.spm.tools.oldseg.output.biascor = 0;
        matlabbatch{step+1}.spm.tools.oldseg.output.cleanup = 0;

        % the TPMs
        
        segmfiles_cell = cellstr([path_to_registered_TPMs 'C1_' filepartition{end-3} '.nii']);
        segmfiles_cell(end+1) = {[path_to_registered_TPMs 'C2_' filepartition{end-3} '.nii']};
        segmfiles_cell(end+1) = {[path_to_registered_TPMs 'C3_' filepartition{end-3} '.nii']};

        matlabbatch{step+1}.spm.tools.oldseg.opts.tpm = reshape(segmfiles_cell,[3 1]);
        matlabbatch{step+1}.spm.tools.oldseg.opts.ngaus = [2
            2
            2
            4];
        matlabbatch{step+1}.spm.tools.oldseg.opts.regtype = 'subj';
        matlabbatch{step+1}.spm.tools.oldseg.opts.warpreg = 1;
        matlabbatch{step+1}.spm.tools.oldseg.opts.warpco = 25;
        matlabbatch{step+1}.spm.tools.oldseg.opts.biasreg = 0.0001;
        matlabbatch{step+1}.spm.tools.oldseg.opts.biasfwhm = 60;
        matlabbatch{step+1}.spm.tools.oldseg.opts.samp = 3;
        matlabbatch{step+1}.spm.tools.oldseg.opts.msk = {''};

        %% Run the batch

        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);

        %[filepath2,name2,ext2] = fileparts(path_to_T2);
        save_batch_to = [path  subject_name{i} filesep Session_name.name filesep];
        save([save_batch_to 'subject_space_segm_preprocessing.mat'],'matlabbatch');
    
        clearvars matlabbatch;        
        
        
        %% Mask the atlas
        [filepath3,name3,ext3] = fileparts(path_to_registered_global_atlas);
        atlasFile= [filepath3 filesep 'atlas_ML_' filepartition{end-3} '.nii'];
        atlasFolder= filepath3;
        segFolder =  path_to_T2;
        volExt ='nii';

        deform = fullfile(segFolder, ['s_' filepartition{end-3} '_' 'seg_sn' '.mat']);
        if isempty(deform)
            error('No deformation file was found')
        else
            disp(['The deformation file that was selected is:' deform] )
        end
        
        c1_file = ['c1s_' filepartition{end-3} '.nii']; % c1s_Subject_01.nii
        c2_file = ['c2s_' filepartition{end-3} '.nii'];
        c3_file = ['c3s_' filepartition{end-3} '.nii'];
        list = [];
        list=[c1_file;c2_file;c3_file];
        
        % segmentFiles will be used for the GM based masking of the atlas
        segmentFiles = cell(size(list,1),1);
        for jj = 1:size(list,1)
            segmentFiles{jj}= fullfile(segFolder, list(jj,:));
        end

        %% create an individual GM-masked AAL atlas in the subject space (C1 = GM segmentation file in needed as input) !!!
        new_Auto_Labelling(segmentFiles, atlasFile, deform, atlasFolder,[]);
        
        
        %% Load the GM-masked and the global atlas (they are in the T2 space
        %now)
        cd(atlasFolder);
        atlas_file_masked = dir(['s_' subject_name{i} '_AtlasGM-masked' '.nii']);
        atlas_file_global = dir(['atlas_ML_' subject_name{i} '.nii']);
        
        if isempty(atlas_file_global.name) || isempty(atlas_file_masked.name)
          error('Registered atlases not found')
        else
          disp('Registered atlases found !')
        end
        
        %% load BOLD of first run
        path_to_BOLD = [path  subject_name{i} filesep Session_name.name filesep run_names{1} filesep];
        cd(path_to_BOLD);
        BOLD_file = dir(['4D_' subject_name{i} '_' run_names{1} '.nii']);
        if isempty(BOLD_file.name)
          error('BOLD could not be defined')
        else
          disp('BOLD was found !')
        end
        BOLD_path = [path_to_BOLD BOLD_file.name];

        %% Reslice the created atlases to bring them from the T2 space to the
        % functional space and save the new images
        atlas_masked_path = [path_to_registered_global_atlas atlas_file_masked.name];
        atlas_global_path = [path_to_registered_global_atlas atlas_file_global.name];        
        
        Vin_fn_masked= atlas_masked_path;
        Vin_fn_global= atlas_global_path;
        Vout_fn= BOLD_path;

        %reslice the masked atlas to functional space
        [Vout, params1]=mapVolumeToVolume(Vin_fn_masked,Vout_fn);
        params1.fname = [path_to_registered_global_atlas 'final_masked_atlas_' subject_name{i} '.nii'];
        spm_write_vol(params1, Vout);
        clearvars Vout;
        
        %reslice the global atlas to functional space
        [Vout, params2]=mapVolumeToVolume(Vin_fn_global,Vout_fn);
        params2.fname = [path_to_registered_global_atlas 'final_global_atlas_' subject_name{i} '.nii'];
        spm_write_vol(params2, Vout);
        clearvars Vout;        
        
        %% NEW: Build masks for nuisance regression
        
        c2vol = spm_vol(segmfiles_cell{2});
        c3vol = spm_vol(segmfiles_cell{3});
            
        c2voldata = spm_read_vols(c2vol);
        c3voldata = spm_read_vols(c3vol);
            
        c2voldata = c2voldata / max(max(max(c2voldata)));
        c3voldata = c3voldata / max(max(max(c3voldata)));
            
            
        c2voldata(c2voldata< 0.90) = 0;
        c3voldata(c3voldata< 0.90) = 0;
            
        [filepath4,name4,ext4] = fileparts(c2vol.fname);
        [filepath5,name5,ext5] = fileparts(c3vol.fname);
            
        c2vol.fname = [path_to_T2, 'w',name4,'_0.9.nii'];
        c3vol.fname = [path_to_T2, 'w',name5,'_0.9.nii'];
            
        spm_write_vol(c2vol, c2voldata);
        spm_write_vol(c3vol, c3voldata);
    
        for r = 1:size(run_names,2)
            
            fprintf('Working on run %s for subject %s\n',run_names{r},subject_name{i} );
            path_to_run = [path subject_name{i} filesep Session_name.name filesep run_names{r} filesep ]; 
            
            %% NEW: Regress out nuisance variates
            paths.f = path_to_run;
            paths.sseg = [path_to_T2];
            
            Detrend = 1;
            removeFirstNscans = 0;
            
            pre_f='arf';
            
            %% Run the regression
            RegressOutNuissance(paths,pre_f,Detrend,removeFirstNscans)
            
            %% Convert the resulted regressed 3D volumes into a 4D nii
            prefix_of_files_to_be_selected = ['Det_' , pre_f];
            f = spm_select('FPList' , path_to_run  , ['^' prefix_of_files_to_be_selected '.*\.nii$']);
            matlabbatch{1}.spm.util.cat.vols = cellstr(f);
            matlabbatch{1}.spm.util.cat.name = ['Det_4D_' subject_name{i} '_' run_names{r} '.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 0;
        
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch );   
            clearvars matlabbatch;
            
            %% Extract the timecourse using these 2 new masks (there are in Bold space now)
            final_masked_atlas_path = params1.fname;
            final_global_atlas_path = params2.fname;

            % read the CLEAN functional volumes 
            cd(path_to_run);
            BOLD_file = dir(['Det_4D_' subject_name{i} '_' run_names{r} '.nii']);
            BOLD_path = [path_to_run BOLD_file.name];
            
            Volumes = spm_vol(BOLD_path);
            s1dim = Volumes(1).dim(1) * Volumes(1).dim(2) * Volumes(1).dim(3);

            %% Read the new atlases (masked and global)
            atlas_masked = spm_read_vols(spm_vol(final_masked_atlas_path));
            atlas_global = spm_read_vols(spm_vol(final_global_atlas_path));

            brainmask_masked = logical(atlas_masked);
            brainmask_global = logical(atlas_global);

            brainmask_masked = reshape(brainmask_masked, [1,s1dim]);
            brainmask_global = reshape(brainmask_global, [1,s1dim]);

            V_masked = zeros(length(Volumes),nnz(brainmask_masked));
            V_global = zeros(length(Volumes),nnz(brainmask_global));

            % iterate over all volumes and extract BOLD TCS based on the
            % atlases
            for ii = 1:length(Volumes)
                % vectorize the functional volume into a vector
                Vprime =  reshape(spm_read_vols(Volumes(ii)), [1,s1dim]);
                Vprime(isnan(Vprime)) = 0;
                Vprime(Vprime < 0) = 0;

                V_masked(ii,:) = Vprime(brainmask_masked);
                V_global(ii,:) = Vprime(brainmask_global);

                if mod(ii,100) == 0
                    disp(ii)
                end
            end
            disp(ii);

            %% Iterate for each ROI in the atlas and compute the mean CLEAN BOLD activity
            Global_Connectome = zeros(length(Volumes), 90);
            Masked_Connectome = zeros(length(Volumes), 90);

            flatten_global_atlas = reshape(atlas_global, [1, s1dim]);
            flatten_masked_atlas = reshape(atlas_masked, [1, s1dim]);

            for j = 1:90
              Global_Connectome(:,j) = mean(V_global(:,flatten_global_atlas(brainmask_global)==j),2);
              Masked_Connectome(:,j) = mean(V_masked(:,flatten_masked_atlas(brainmask_masked)==j),2);
            end


            fprintf('Subject %s , run %s  is finished.\n', subject_name{i},run_names{r});
            % Save the functional connectomes
            if run_names{r} == 'func_01'
                csvwrite( [save_global_connectomes_to 'pre' filesep 'global' filesep 'pre_global' '_' subject_name{i} '.csv'] , Global_Connectome)
                csvwrite( [save_masked_connectomes_to 'pre' filesep 'masked' filesep 'pre_masked' '_' subject_name{i} '.csv'] , Masked_Connectome)
            elseif run_names{r} == 'func_02'
                csvwrite( [save_global_connectomes_to 'post' filesep 'global' filesep 'post_global' '_' subject_name{i} '.csv'] , Global_Connectome)
                csvwrite( [save_masked_connectomes_to 'post' filesep 'masked' filesep 'post_masked' '_' subject_name{i} '.csv'] , Masked_Connectome)                
            end
            
        end
   
        
end
cd(go_back);
disp("ALL DONE!!");
return
    

