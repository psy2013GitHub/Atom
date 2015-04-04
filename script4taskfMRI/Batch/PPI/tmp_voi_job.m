

%----------------------------------------------------
% Created at 2014Apr05 12:27:12
%----------------------------------------------------




clear matlabbatch;
matlabbatch{1}.spm.util.voi.spmmat = {'/home/dengzhou/WM_fMRI/GLM_1stlevel/sub9_task_suc/SPM.mat'};
matlabbatch{1}.spm.util.voi.adjust = 1;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name = 'VC';
matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {'/home/dengzhou/WM_fMRI/GLM_1stlevel/sub9_task_suc/SPM.mat'};
matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 4;
matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 5.000000e-02;
matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = [48.0 -81.0 27.0];
matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 6;
matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.roi{3}.mask.image = {'/home/dengzhou/WM_fMRI/GLM_1stlevel/sub9_task_suc/mask.img'};
matlabbatch{1}.spm.util.voi.roi{3}.mask.threshold = 0.5;
matlabbatch{1}.spm.util.voi.expression = 'i1&i2&i3';
