clear spm_pipeline

spm_pipeline.Batch_Dir=pwd;

spm_pipeline.Batch_SliceTiming_Dir=[spm_pipeline.Batch_Dir,filesep,'SliceTiming'];
spm_pipeline.Batch_Realign_Dir=[spm_pipeline.Batch_Dir,filesep,'Realign'];
spm_pipeline.Batch_Normalise_Dir=[spm_pipeline.Batch_Dir,filesep,'Normalise'];
spm_pipeline.Batch_1stlevel_Dir=[spm_pipeline.Batch_Dir,filesep,'1stlevel'];
spm_pipeline.Batch_Contrast_Dir=[spm_pipeline.Batch_Dir,filesep,'Contrast'];
spm_pipeline.Batch_Smooth_Dir=[spm_pipeline.Batch_Dir,filesep,'Smooth'];


spm_pipeline.Utility_Dir=[spm_pipeline.Batch_Dir,filesep,'Utility'];
spm_pipeline.Utility_Del_files=[spm_pipeline.Utility_Dir,filesep,'Del_files'];
spm_pipeline.Utility_Move_con=[spm_pipeline.Utility_Dir,filesep,'Move_con'];
