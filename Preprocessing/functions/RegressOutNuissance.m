function [] = RegressOutNuissance(paths,pre_f,Detrend,removeFirstNscans)
%% this part is useful only in case of WM-CSF regression [in fMRI individual space (usual case), and for masking the CSF-WM masks with realigned functional Vol]
% path to mean realigned functional image - takes only header so use
    % any functional image
    
    %HERE realigned image (not smoothed)
    if pre_f(1)=='s' %smoothing was done
        %pre_r=pre_f(3:end);
        pre_r='mean';
    else
        pre_r=pre_f;
    end
        
    [fMeanFN_noPath, foosubdirs] = spm_select('List',paths.f,['^' pre_r '.*\.nii$']);
    if isempty(fMeanFN_noPath)
        [fMeanFN_noPath, foosubdirs] = spm_select('List',paths.f,['^' pre_r '.*\.img$']);
    end
    fMeanFN_noPath=fMeanFN_noPath(1,:); % take first file
    fMean_fn=fullfile(paths.f,fMeanFN_noPath);
    rVol=spm_read_vols(spm_vol(fMean_fn));%read volume
    
    %% Regressing out CSF
    fprintf('Extracting CSF...');
    % path to CSF segmentation file
    [segFN_noPath, foosubdirs] = spm_select('List',paths.sseg,'^wC3_.*\.nii$'); % new file with overlap with mask - in individual space case
    % see eg preprocROItimeCourse for generation
    sSeg_fn=fullfile(paths.sseg,segFN_noPath);
    % downsample stuctural-res segmented CSF to functional res
    CSFV=mapVolumeToVolume(sSeg_fn,fMean_fn);
    CSFmaskLidx=CSFV;%>0.9; %previously thresholded as no mask used
    clear CSFV;
    fprintf('\n');
    
    
    %% Regressing out WM
    fprintf('Extracting WM...')
    %individual space
    [segFN_noPath, foosubdirs] = spm_select('List',paths.sseg,'^wC2_.*\.nii$');
    sSeg_fn=fullfile(paths.sseg,segFN_noPath);
    % downsample stuctural-res segmented WM to functional res
    WMV = mapVolumeToVolume(sSeg_fn,fMean_fn); % map seg to mean voxel space
    WMmaskLidx=WMV;%>0.9;
    clear WMV;
    fprintf('\n');
    
    
    %% Exract the time courses
    try
        %fVolsFNlist=getImageFNinAcqOrder(paths.f,pre_f,'nii',[]);
        fVolsFNlist=dir(fullfile(paths.f,[pre_f,'*.nii']));%,'tryParsing',tryParsingFileOrder);
        fVolsFNlist = struct2cell(fVolsFNlist); 
        fVolsFNlist = fVolsFNlist(1,:);
        
    catch
        fVolsFNlist=getImageFNinAcqOrder2(paths.f,pre_f,'img',[]);
    end
    
    %fVolsFNlist_fullpath=cellfun(@(x) fullfile(paths.f,x), fVolsFNlist,'UniformOutput',false);
    % read all headers and files in temporal order into a 4-D array
    fVolsFNlist_fullpath=cellfun(@(x) fullfile(paths.f,x), fVolsFNlist,'UniformOutput',false);
    V0i=spm_vol(fVolsFNlist_fullpath);
    
    V0idx=1:length(V0i);
    % clear and preallocate to avoid fragmenting
    clear V0;
    V0=zeros(V0i{1}.dim(1),V0i{1}.dim(2),V0i{1}.dim(3),length(V0i),'single');
    
    %% read all volumes
    fprintf('Reading all volumes\n');
    for iter=1:length(V0i)
        V0(:,:,:,V0idx(iter)) = spm_read_vols(V0i{iter});
    end
    
    %% Detrending
    AllVolume=reshape(V0,[],size(V0,4))'; % time x space
    CUTNUMBER=10; % cut space in parts
    
    % Detrend
    if Detrend==1
        %AllVolume=detrend(AllVolume);
        if Detrend, fprintf('Detrending'); end
        SegmentLength = ceil(size(AllVolume,2) / CUTNUMBER);
        for iCut=1:CUTNUMBER
            if iCut~=CUTNUMBER
                Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
            else
                Segment = (iCut-1)*SegmentLength+1 : size(AllVolume,2);
            end
            if Detrend, AllVolume(:,Segment) = detrend(AllVolume(:,Segment)); end
            fprintf('.');
        end
        fprintf('\n');
    end
            
    %% demeaning of each voxel (Ns space).
%    % for j = 1:size(AllVolume,2)
%        % AllVolume(:,j) = AllVolume(:,j) - mean(AllVolume(:,j));
%    % end
    
    % CHANGES FROM V0 TO V0detrended
    V0detrended = reshape(AllVolume',size(V0,1),size(V0,2),size(V0,3),size(V0,4)); % reshape back to 4D shape

    disp(['Computing average CSF signal...']);
    CSFmaskLidx(~logical(rVol))=0;
    % average signal of all CSF voxels corresponding to this block
    [Ci,Cj,Ck]=ind2sub(size(V0(:,:,:,1)),find(CSFmaskLidx));
    %CSFsignal=V0(Ci,Cj,Ck,V0idx);
    nCSFvoxels=numel(Ci); % also, sum(CSFmaskLidx)
    tcCSF=zeros(1,1,1,size(V0,4));
    %[csfx1,csfx2,csfx3]=ndgrid(1:szCSFvoxels(1),1:szCSFvoxels(2),1:szCSFvoxels(3));
    for idx=1:numel(Ci)
        tcCSF=tcCSF+V0detrended(Ci(idx),Cj(idx),Ck(idx),:);
    end
    CSFavg=squeeze(tcCSF/nCSFvoxels)';
    CSFavg_demean=CSFavg-mean(CSFavg);
    clear tcCSF;

    %%
    disp('Computing average WM signal...');
    % average signal of all WM voxels corresponding to this block
    WMmaskLidx(~logical(rVol))=0;
    [Ci,Cj,Ck]=ind2sub(size(V0(:,:,:,1)),find(WMmaskLidx));
    %CSFsignal=V0(Ci,Cj,Ck,V0idx);
    nWMvoxels=numel(Ci); % also, sum(CSFmaskLidx)
    tcWM=zeros(1,1,1,size(V0,4));
    %[csfx1,csfx2,csfx3]=ndgrid(1:szCSFvoxels(1),1:szCSFvoxels(2),1:szCSFvoxels(3));
    for idx=1:numel(Ci)
        tcWM=tcWM+V0detrended(Ci(idx),Cj(idx),Ck(idx),:);
    end
    WMavg=squeeze(tcWM/nWMvoxels)';
    WMavg_demean=WMavg-mean(WMavg);
    clear tcWM;
    
    %% Regress out CSF and WM using GLM
    CSFavg=ones(size(V0idx)); % need length to initialize X
    X=[ones(numel(CSFavg),1) [1:numel(CSFavg)]'/numel(CSFavg) ...
        [1:numel(CSFavg)].^2'/(numel(CSFavg)^2)]; % [const lin quadratic]

    X=[X,WMavg_demean'];
    X=[X,CSFavg_demean'];
    
    % LOAD MOTION PARAMS
    mot = dir(fullfile(paths.f,'rp*.txt'));
    goodvols = dir(fullfile(paths.f,'good_volumes_FD0.5_DVARS30.mat'));
    goodvols = load([goodvols.folder ,'/', goodvols.name]);
    goodvols = goodvols.good;
    
    Cov=load(fullfile(paths.f,mot(1).name));
    Cov = Cov(goodvols,:);
    X=[X,Cov(removeFirstNscans+1:end,1:6)];
%     
    V0Cov=zeros(size(V0));
    means = zeros(size(V0));
    for iii=1:size(V0,1)
        fprintf('.');
        for jjj=1:size(V0,2)
            for kkk=1:size(V0,3)
                [beta,res,SSE,SSR,T] = y_regress_ss(squeeze(V0detrended(iii,jjj,kkk,V0idx)),X);
                beta4(iii,jjj,kkk)=beta(4);
                beta5(iii,jjj,kkk)=beta(5);
                %beta6(i,j,k)=beta(6);
                %beta7(i,j,k)=beta(7);
                T4(iii,jjj,kkk)=T(4);
                T5(iii,jjj,kkk)=T(5);
                %T6(i,j,k)=T(6);
                %T7(i,j,k)=T(7);
%               V0Cov(i,j,k,V0idx)=res;%without mean
%               means(i,j,k) = mean(squeeze(V0(i,j,k,V0idx)));
                V0Cov(iii,jjj,kkk,V0idx)= res + mean(squeeze(V0(iii,jjj,kkk,V0idx))); %put the mean signal back
            end
        end
    end
    V0Cov(isnan(V0Cov))=0;
    V0=V0Cov;

    %% write V0Cov
    mV0i=V0i;

    fprintf('\nWriting new Clean volumes...');
    
    %voxel_size = diag(mV0i{1}.mat);
    %voxel_size = voxel_size(1:end-1)';

    %voxel_shift = mV0i{1}.mat(:,4);
    %voxel_shift = voxel_shift(1:end-1)';
        
    for y=1:size(V0Cov,4)
        %change hdr names
        s=strfind(V0i{y}.fname,'/');
        last=s(size(s,2));
        mV0i{y}.fname=strcat(V0i{y}.fname(1:last),'Det_',V0i{y}.fname(last+1:end));
        
        %write volumes
        spm_write_vol(mV0i{y}, V0(:,:,:,y));

        %Vres = double(reshape(V0Cov(:,:,:,y), mV0i{y}.dim));

        % Creates a new NIFTI with the right parameters
        %tmp_NIFTI = make_nii(Vres,voxel_size,-voxel_shift./voxel_size);
        %tmp_NIFTI.hdr.dime.datatype=64;
        %tmp_NIFTI.hdr.dime.bitpix=64;
            
        %save_nii(tmp_NIFTI,mV0i{y}.fname);
        fprintf('.');
    end
     
    fprintf('Done.')
end