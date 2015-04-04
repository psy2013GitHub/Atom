

DataDir='/home/dengzhou/WM_fMRI/PPI_1stlevel';
DestDir='/home/dengzhou/WM_fMRI/PS/PPI_1stlevel';
MultiRuns_flag=0;

conspec.title='';
conspec.contrast=1;
conspec.thresdesc='none';
conspec.thresh=0.001;
conspec.extend=0;
conspec.print=1;

%--------------------------------------------------------------------------
Batch_ResultDisp_Dir='.';
Batch_ResultDisp_func(DataDir,DestDir,Batch_ResultDisp_Dir,MultiRuns_flag, conspec)