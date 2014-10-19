

% Case-insensitive
% After entering the parameters, use icatb_batch_file_run(inputFile); 

%========================================================================%
modalityType = 'fMRI';

%- Type of analysis
% Options are 1 and 2.
% 1 - Regular Group ICA
% 2 - Group ICA using icasso
which_analysis = 1;

%% ICASSO options.
% This variable will be used only when which_analysis variable is set to 2.
icasso_opts.sel_mode = 'randinit';  % Options are 'randinit', 'bootstrap' and 'both'
icasso_opts.num_ica_runs = 5; % Number of times ICA will be run
% Most stable run estimate is based on these settings. 
icasso_opts.min_cluster_size = 2; % Minimum cluster size
icasso_opts.max_cluster_size = 15; % Max cluster size. Max is the no. of components

%- Group PCA performance settings.
% Options are:
% 1 - Maximize Performance
% 2 - Less Memory Usage
% 3 - User Specified Settings
perfType = 1;

%- Enter location (full file path) of the image file to use as mask or use Default mask which is []
p=which('rest'); p=fileparts(p); restmask=[p,filesep,'mask',filesep,'BrainMask_05_61x73x61.img'];
maskFile = restmask;


%========================================================================%
%- Design matrix selection
% 1. 'no' - means no design matrix.
% 2. 'same_sub_same_sess' - same design over subjects and sessions
% 3. 'same_sub_diff_sess' - same design matrix for subjects but different
% over sessions
% 4. 'diff_sub_diff_sess' - means one design matrix per subject.
keyword_designMatrix = 'no';

% specify location of design matrix here if 'same_sub_same_sess' or 'same_sub_diff_sess'
OnedesignMat = 'C:\MATLAB6p5p2\work\Example Subjects\Visuomotor_data\SPM.mat';


%========================================================================%
%- Enter directory to put results of analysis
outputDir = 'F:\Users\dengzhou\modified\op\ICAs';

%- Enter Name (Prefix) Of Output Files
prefix = 'test';


%========================================================================%
%- Group PCA Type. Used for analysis on multiple subjects and sessions.
% Options are 'subject specific' and 'grand mean'. 
%   a. Subject specific - Individual PCA is done on each data-set before group
%   PCA is done.
%   b. Grand Mean - PCA is done on the mean over all data-sets. Each data-set is
%   projected on to the eigen space of the mean before doing group PCA.
%
% NOTE: Grand mean implemented is from FSL Melodic. Make sure that there are
% equal no. of timepoints between data-sets.
%
group_pca_type = 'subject specific';

%- Back reconstruction type. Options are 'str(Spatial-Temporal Correlation)','GICA', 'GICA 3'
backReconType = 'gica';

%- Data Pre-processing options
% 1 - Remove mean per time point
% 2 - Remove mean per voxel
% 3 - Intensity normalization
% 4 - Variance normalization
preproc_type = 1;


%- PCA Type. Also see options associated with the selected pca option. EM
% PCA options and SVD PCA are commented.
% Options are 1, 2 and 3
% 1 - Standard 
% 2 - Expectation Maximization
% 3 - SVD
pcaType = 1;

%------------------------------------------------------------------------%
%- PCA options (Standard)
% a. Options are yes or no
% 1a. yes - Datasets are stacked. Memory-demaned
% 2a. no - A pair of datasets are loaded at a time. Lower
pca_opts.stack_data = 'yes';

% b. Options are full or packed of covariance matrix.
% 1b. full
% 2b. packed
pca_opts.storage = 'full';

% c. Options are double or single precison.
% 1c. double
% 2c. single
pca_opts.precision = 'double';

% d. Type of eigen solver. Options are selective or all
% 1d. selective - Selective eigen solver is used. If there are convergence issues, use option all.
% 2d. all - All eigen values are computed. This might run very slow if you are using packed storage. Use this only when selective option doesn't converge.
pca_opts.eig_solver = 'selective';

%------------------------------------------------------------------------%
% %% PCA Options (Expectation Maximization)
% % a. Options are yes or no
% % 1a. yes - Datasets are stacked. This option uses lot of memory depending
% % on datasets, voxels and components.
% % 2a. no - A pair of datasets are loaded at a time. This option uses least
% % amount of memory and can run very slower if you have very large datasets.
% pca_opts.stack_data = 'yes';
% 
% % b. Options are double or single.
% % 1b. double - Double precision is used
% % 2b. single - Floating point precision is used.
% pca_opts.precision = 'single';
% 
% % c. Stopping tolerance 
% pca_opts.tolerance = 1e-4;
% 
% % d. Maximum no. of iterations
% pca_opts.max_iter = 1000;


% %% PCA Options (SVD)
% % a. Options are double or single.
% % 1a. double - Double precision is used
% % 2a. single - Floating point precision is used.
% pca_opts.precision = 'single';

% % b. Type of eigen solver. Options are selective or all
% % 1b. selective - svds function is used.
% % 2b. all - Economy size decomposition is used.
% pca_opts.solver = 'selective';


%- Maximum reduction steps you can select is 3
% Note: This number will be changed depending upon the number of subjects
% and sessions
numReductionSteps = 3;

%- Batch Estimation. If 1, Automatic estimation of  the components and the corresponding PC numbers are associated
% Options are 1 or 0
doEstimation = 1; 

%- MDL Estimation options. This variable will be used only if doEstimation is set to 1.
% Options are 'mean', 'median' and 'max' for each reduction step. The length of cell is equal to the no. of data reductions used.
estimation_opts.PC1 = 'mean';
estimation_opts.PC2 = 'mean';
estimation_opts.PC3 = 'mean';

%- Number of pc to reduce each subject down to at each reduction step  
numOfPC1 = 30;
numOfPC2 = 25;
numOfPC3 = 20;

%========================================================================%
%- Scale the Results. Options are 0, 1, 2, 3 and 4
% 0 - Don't scale
% 1 - Scale to Percent signal change
% 2 - Scale to Z scores
% 3 - Normalize spatial maps using the maximum intensity value and multiply timecourses using the maximum intensity value
% 4 - Scale timecourses using the maximum intensity value and spatial maps using the standard deviation of timecourses
scaleType = 2;


%========================================================================%
%- 'Which ICA Algorithm Do You Want To Use';
% Note: Use only one subject and one session for Semi-blind ICA. Also specify atmost two reference function names
% 1 means infomax, 2 means fastICA, etc.
algoType = 1;


%? Specify atmost two reference function names if you select Semi-blind ICA algorithm.
% Reference function names can be acessed by loading SPM.mat in MATLAB and accessing structure SPM.xX.name.
refFunNames = {'Sn(1) right*bf(1)', 'Sn(1) left*bf(1)'};


%? Specify spatial reference files for multi-fixed ICA
refFiles = {which('ref_default_mode.nii'), which('ref_left_visuomotor.nii'), which('ref_right_visuomotor.nii')};

%- ICA Options - Name by value pairs in a cell array. Options will vary depending on the algorithm. See icatb_icaOptions for more details. Some options are shown below.
% Infomax -  {'posact', 'off', 'sphering', 'on', 'bias', 'on', 'extended', 0}
% FastICA - {'approach', 'symm', 'g', 'tanh', 'stabilization', 'on'}
icaOptions = {'posact', 'off', 'sphering', 'on', 'bias', 'on', 'extended', 1};


%=======================================================================%
%% There are three ways to enter the subject data
% options are 1, 2, 3 or 4
dataSelectionMethod = 4;

%- Method 3 (Uses Regular expressions)

% Input data directory name
input_directory_name = '\nfs\data';

% Subject directory regular expression. This variable can have nested paths
% like Sub01_vis\Study1. To match this Sub\w+; Study\w+ regular expression can be used where semi-colon
% is used as a path separator. If there are no subject directories inside the input directory, leave it as empty like ''
subject_dir_regexp = 'Sub\w+';

% Session directory regular expression. This variable cannot have nested
% paths. If there are no session directories inside subject directories, leave it as empty.
session_dir_regexp = '';

% Data file pattern. Use wild card for this and not regular expression.
data_file_pattern = 'ns*.img';

% File numbers to include. Leave it as empty if you want to include all of
% them.
file_numbers_to_include = [];

% SPM stats directory name relative to subject or session directories. Use this only when you specify
% 'diff_sub_diff_sess' as the value for keyword_designMatrix variable. GIFT
% will first search in the subject directories and later session
% directories to find SPM.mat files
spm_stats_dir = '';

%%%%%%%% End for Method 3 %%%%%%%%%%%%

%- Method 4
% Input data file pattern for data-sets must be in a cell array. The no. of rows of cell array correspond to no. of subjects
% and columns correspond to sessions. In the below example, there are 3
% subjects and 1 session. If you have multiple sessions, please see
% Input_data_subjects_2.m file.
input_data_file_patterns = {
                            'F:\MultimodeAnalysisFC\FC_EPI\FunImgNormalizedSmoothedDetrendedFilteredCovremoved\rest_s12009\*.img';...
                            'F:\MultimodeAnalysisFC\FC_EPI\FunImgNormalizedSmoothedDetrendedFilteredCovremoved\rest_s12010\*.img';...
                            'F:\MultimodeAnalysisFC\FC_EPI\FunImgNormalizedSmoothedDetrendedFilteredCovremoved\rest_s12011\*.img';...
                            'F:\MultimodeAnalysisFC\FC_EPI\FunImgNormalizedSmoothedDetrendedFilteredCovremoved\rest_s12012\*.img';...
                            'F:\MultimodeAnalysisFC\FC_EPI\FunImgNormalizedSmoothedDetrendedFilteredCovremoved\rest_s12013\*.img';...
                            'F:\MultimodeAnalysisFC\FC_EPI\FunImgNormalizedSmoothedDetrendedFilteredCovremoved\rest_s12014\*.img';
                            };

% Input for design matrices will be used only if you have a design matrix
% for each subject i.e., if you have selected 'diff_sub_diff_sess' for
% variable keyword_designMatrix.
input_design_matrices = {'C:\MATLAB6p5p2\work\Example Subjects\Visuomotor_data\SPM.mat';
    'C:\MATLAB6p5p2\work\Example Subjects\Visuomotor_data\SPM.mat';
    'C:\MATLAB6p5p2\work\Example Subjects\Visuomotor_data\SPM.mat'};

% Enter no. of dummy scans to exclude from the group ICA analysis. If you have no dummy scans leave it as 0.
dummy_scans = 0;

%%%%%%%% End for Method 4 %%%%%%%%%%%%

