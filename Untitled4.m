

clc;clear LadyFirst;

%--------------------------------------------------------------------------
%- Modality
LadyFirst.Modality='fMRI';

%- Steps
%-> LadyFirst.Steps={'Dcm2Nii','RemoveTimePoints','SliceTiming','Realign','T1Normalize','Smooth',''};
LadyFirst.Steps={'CovRegress','Filter','Detrend'};

%- Input
LadyFirst.InputDir.Functionals={'E:\LadyFirst\Rest\SmoothNormalizeRealignSliceTimingFunImg'};
LadyFirst.InputDir.Anatomicals={'E:\LadyFirst\3D\oldSegment'}; % if postpreprocess step, then specify it with aCompCor masks direct
LadyFirst.InputDir.Motions={'E:\LadyFirst\Rest\MotionParameter'};
% filter for files & directorys
LadyFirst.dirFilter.Anatomicals={'*'}; % if not advanced user, not to change !
LadyFirst.fileFilter.Anatomicals=[{'wc1*'},{'wc2*'},{'wc3*'}];
LadyFirst.dirFilter.Functionals={'*'}; % if not advanced user, not to change !
LadyFirst.fileFilter.Functionals={'swra*'};

%- Output
LadyFirst.OutputDir.Functionals={'E:\LadyFirst\Rest'};
LadyFirst.OutputDir.Anatomicals={'E:\LadyFirst\3D\AnaImg\'};

%- Remove Time Points
LadyFirst.RemoveTimePoints.nFirstImgDiscard=5;

%- SliceTiming
LadyFirst.SliceTiming.SliceNumber=32;
LadyFirst.SliceTiming.TR=2;
LadyFirst.SliceTiming.SliceOrder=[2:2:32,1:2:31];
LadyFirst.SliceTiming.ReferenceSlice=32;

%- Realign
LadyFirst.Realign.Registaer2Mean=0;

%- Smooth
LadyFirst.Smooth.FWHM=[6,6,6];


%- Covariables
LadyFirst.Covariables.Names              ={'WM','CSF','Motion'};
LadyFirst.Covariables.Dim                =[5,   5,    6]; % 0 for mean, >1 for PCA dimensions
LadyFirst.Covariables.Derivative         =[0,   0,    1]; % derivative   
LadyFirst.Covariables.Thresh             =[0.99,0.99, 1];
LadyFirst.Covariables.Erode              =[1,   0,    1]; % usually, for 'CSF' no need to erode
LadyFirst.Covariables.Path               ={}; % leave empty uless other wise stated

%- Filter
LadyFirst.Filter.Type='bandpass'; % 'bandpass'/'butter'
LadyFirst.Filter.FrequencyBand=[0.01,0.1];

%- Derend
LadyFirst.Detrend.PolyTrend=1;

%--------------------------------------------------------------------------
%- Final Run
LadyFirst=LadyFirst_Run(LadyFirst);