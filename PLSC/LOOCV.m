% The script includes: 
%   1. Call of a script with PLS inputs and their description
%   2. Call of the functions to run PLS and plot the results
%%
clear all; clc;
addpath('./myPLS_functions')
addpath('./RotmanBaycrest')
addpath('./misc')
addpath('./Inputs')
% ZSCORING: replaces the measurement unit with "number of standard deviations" away from the mean.

% SETUP: Define all the inputs
% Modify this script to setup your PLS analysis
myPLS_inputs

% Check all inputs for validity
% check setup before running PLS
[input,pls_opts,save_opts] = myPLS_initialize(input,pls_opts,save_opts);

%% Get constants
% Number of subjects
nSubj = size(input.brain_data,1); 
% Number of behavior/design variables
nBehav = size(input.behav_data,2);
% Number of imaging variables
nImg = size(input.brain_data,2);  

%% Get initial imaging matrix X0
X0_all = input.brain_data;

% Generate the behavior/contrast/interaction matrix Y0
[Y0_all,design_names,nDesignScores] = myPLS_getY(pls_opts.behav_type,...
    input.behav_data,input.grouping,input.group_names,input.behav_names);

disp('... Input data information ...')
disp(['Number of observations (subjects): ' num2str(nSubj)]);
disp(['Number of brain measures (voxels/connections): ' num2str(nImg)]);
disp(['Number of design measures (behavior/contrasts): ' num2str(nDesignScores)]);
disp(' ')


%% Full run of PLSC
X0 = X0_all(:,:);
Y0 = Y0_all(:);

%% Normalize input data matrices X and Y 
% (if there is a group contrast in Y, groups won't be considered in any case)
% standardization of each column individually

current_mean_X0 = mean(X0,1);
current_std_X0 = std(X0,[],1);

current_mean_Y0 = mean(Y0,1);
current_std_Y0 = std(Y0,[],1);

X = zeros(size(X0));
for line=1:size(X0,1)
    X(line,:) = (X0(line,:) - current_mean_X0) ./ current_std_X0;
end

Y = zeros(size(Y0));
for line=1:size(Y0,1)
    Y(line,:) = (Y0(line,:) - current_mean_Y0) ./ current_std_Y0;
end

% Cross-covariance matrix
if pls_opts.grouped_PLS==0; disp('empty pls_opts.grouped_PLS'); end
% Cross-covariance matrix using the z-scores across subjects data
R = myPLS_cov(X,Y,input.grouping,pls_opts.grouped_PLS); % computes R = Y'*X;
assert(sum(Y'*X - R) ==0)

% Singular value decomposition
[U,S,V] = svd(R,'econ');
% U: behavior/design saliences
% S: singular values
% V: Imaging saliences

% Number of latent (LCs)
nLC = min(size(S));
% Amount of covariance explained by each LC
explCovLC = (diag(S).^2) / sum(diag(S.^2));

% Compute PLS scores & loadings
% see notes for estimation details
% Lx = X * V / Ly = Y * U / corr_Lx_X = corr(Lx,X)' / corr_Ly_Y = corr(Ly,Y)'
[Lx,Ly,corr_Lx_X,corr_Ly_Y,corr_Lx_Y,corr_Ly_X] = ...
    myPLS_get_PLS_scores_loadings(X,Y,V,U,input.grouping,pls_opts);
assert(sum(X*V - Lx) == 0)
assert(sum(corr(Lx,X)'-corr_Lx_X) ==0)
assert(sum(corr(Ly,Y)'-corr_Ly_Y) ==0)

% Permutation testing for LC statistical significance (p-value)
% (only rows (subjects) of Y are permuted in each iteration)
Sp_vect = myPLS_permutations(X,Y,U,input.grouping,pls_opts);

% Compute the p-values from the permutation null distribution
LC_pvals = myPLS_get_LC_pvals(Sp_vect,S,pls_opts);


% Bootstrapping to test stability of PLS loadings
%boot_results = myPLS_bootstrapping(X0,Y0,U,V,input.grouping,pls_opts);

% TODO: compute bootstrapping stats and maybe remove the original sampling
% data, depending on the type (or size?) of the data (for voxelwise data,
% saving the original sampling data would take up too much space)

% Save all result variables in a structure
%res.R_full = R; % Cross-covariance matrix using the z-scores across subjects data
res.U_full = U; % behavior Latent variables
res.S_full = S; % singular values
res.V_full = V; % imaging Latent variables
%res.explCovLC_full = explCovLC; % explained co-variance by LC
res.LC_pvals_full = LC_pvals; % p-values of LCs after permutations testing
res.Lx_full = Lx; % PLS imaging scores / projections
res.Ly_full = Ly; % PLS behavior scores / projections
%res.LC_img_loadings_full = corr_Lx_X; %%% can be changed to corr_Ly_X
%res.LC_behav_loadings_full = corr_Ly_Y; %%% can be changed to corr_Lx_Y
%res.Sp_vect_full = Sp_vect; % singular values obtained by permutation testing (by permuting rows of Y)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOCV loop
for LOO = 1:nSubj
    disp(LOO);
    % training samples
    X0 = X0_all(1:end ~= LOO , :);
    Y0 = Y0_all(1:end ~= LOO);
    % left-out sample
    X0_left_out = X0_all(LOO , :);
    Y0_left_out = Y0_all(LOO);

    %% Normalize input data matrices X and Y 
    % (if there is a group contrast in Y, groups won't be considered in any case)
    % standardization of each column individually

    training_mean_X0 = mean(X0,1);
    training_std_X0 = std(X0,[],1);

    training_mean_Y0 = mean(Y0,1);
    training_std_Y0 = std(Y0,[],1);

    X = zeros(size(X0));
    for line=1:size(X0,1)
        X(line,:) = (X0(line,:) - training_mean_X0) ./ training_std_X0;
    end

    Y = zeros(size(Y0));
    for line=1:size(Y0,1)
        Y(line,:) = (Y0(line,:) - training_mean_Y0) ./ training_std_Y0;
    end

    % Cross-covariance matrix
    if pls_opts.grouped_PLS==0; disp('empty pls_opts.grouped_PLS'); end
    % Cross-covariance matrix using the z-scores across subjects data
    input.grouping_ = input.grouping(1:end ~= LOO);
    R = myPLS_cov(X,Y,input.grouping_,pls_opts.grouped_PLS); % computes R = Y'*X;
    assert(sum(Y'*X - R) ==0)

    % Singular value decomposition
    [U,S,V] = svd(R,'econ');
    % U: behavior/design saliences
    % S: singular values
    % V: Imaging saliences

    % Number of latent (LCs)
    nLC = min(size(S));
    % Amount of covariance explained by each LC
    explCovLC = (diag(S).^2) / sum(diag(S.^2));

    % Compute PLS scores & loadings
    % see notes for estimation details
    % Lx = X * V / Ly = Y * U / corr_Lx_X = corr(Lx,X)' / corr_Ly_Y = corr(Ly,Y)'
    [Lx,Ly,corr_Lx_X,corr_Ly_Y,corr_Lx_Y,corr_Ly_X] = ...
        myPLS_get_PLS_scores_loadings(X,Y,V,U,input.grouping_,pls_opts);
    assert(sum(X*V - Lx) == 0)
    assert(sum(corr(Lx,X)'-corr_Lx_X) ==0)
    assert(sum(corr(Ly,Y)'-corr_Ly_Y) ==0)

    % Permutation testing for LC statistical significance (p-value)
    % (only rows (subjects) of Y are permuted in each iteration)
    Sp_vect = myPLS_permutations(X,Y,U,input.grouping_,pls_opts);

    % Compute the p-values from the permutation null distribution
    LC_pvals = myPLS_get_LC_pvals(Sp_vect,S,pls_opts);


    %% Predict using the left-out sample
    % apply the training normalization + PLSC component to it
    X0_left_out = ( X0_left_out - training_mean_X0 ) ./ training_std_X0;
    Y0_left_out = ( Y0_left_out - training_mean_Y0 ) ./ training_std_Y0;

    %% Left-out sample scores
    Lx_out_sample = X0_left_out * V; 
    Ly_out_sample = Y0_left_out * U; 

    % Save all result variables in a structure
    %res.R_cv = R; % Cross-covariance matrix using the z-scores across subjects data
    res.U_cv = U; % behavior Latent variables
    res.S_cv = S; % singular values
    res.V_cv = V; % imaging Latent variables
    %res.explCovLC_cv = explCovLC; % explained co-variance by LC
    res.LC_pvals_cv = LC_pvals; % p-values of LCs after permutations testing
    
    myfield = (strcat('Lx_CV_fold_',num2str(LOO)));
    res.(myfield) = Lx_out_sample; % PLS imaging scores / projections
    Lxs_cv(LOO) = Lx_out_sample;
    
    myfield = (strcat('Ly_CV_fold_',num2str(LOO)));
    res.(myfield) = Ly_out_sample;  % PLS behavior scores / projections
    Lys_cv(LOO) = Ly_out_sample;
end

% Save results to CSV
csvwrite('../LOOCV/Lx_CV.csv', Lxs_cv);
csvwrite('../LOOCV/Ly_CV.csv', Lys_cv);
csvwrite('../LOOCV/Lx_full.csv', res.Lx_full);
csvwrite('../LOOCV/Ly_full.csv', res.Ly_full);
