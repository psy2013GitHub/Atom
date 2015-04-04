


RawDir='/home/dengzhou/WM_fMRI/TaskImg';
GLMDir='/home/dengzhou/WM_fMRI/GLM_1stlevel';
nii_filter='^swr.*\.(img|nii)$';
PPIDir='/home/dengzhou/WM_fMRI/PPI_1stlevel';
MultiRuns_flag=0;
MergeRuns_flag=1;


% voi
clear voi;
voi{1}.ajust_con=1;
voi{1}.sesseion=1;
voi{1}.name='VC';
voi{1}.roi{1}.spm.con=4;
voi{1}.roi{1}.spm.thresdesc='none';
voi{1}.roi{1}.spm.thresh=0.05;
voi{1}.roi{1}.spm.extend=0;
% voi{1}.roi{2}.diy_mask_peak.mask_file='/home/dengzhou/WM_fMRI/Mask/WorkingMemory_JN2012/3mm/3mm_VC.nii';
% voi{1}.roi{2}.diy_mask_peak.spm_con=4;
% voi{1}.roi{2}.diy_mask_peak.sphere.radius=6;
% voi{1}.expression='i1&i2';
voi{1}.roi{2}.sphere.center=[48,-81,27];
voi{1}.roi{2}.sphere.radius=6;
voi{1}.roi{2}.sphere.move.fixed= 1;
voi{1}.expression='i1&i2';

% ppi variable construction
ppi_specify.u=[1,1,0.5;2,1,0.5;3,1,-1];
ppi_specify.name='VC(2back-Fixation)';
ppi_specify.disp=1;

% Design .e. ppi_specify;
ppi_specify.timing.units='scans';
ppi_specify.timing.RT=2;
ppi_specify.timing.fmri_t=32;
ppi_specify.timing.fmri_t0=16;
% regress motion ?
ppi_specify.regress_motion=1;





% -------------------------------------------------------------------------
p=fileparts(pwd); Batch_PPI_Dir=[p,filesep,'PPI']; Batch_1stlevel_Dir=[p,filesep,'1stlevel']; Batch_Contrast_Dir=[p,filesep,'Contrast'];
PPI_batch_func(RawDir, GLMDir, PPIDir, nii_filter, Batch_PPI_Dir, Batch_1stlevel_Dir, Batch_Contrast_Dir, voi, MultiRuns_flag, MergeRuns_flag, ppi_specify);