clc;clear Atom;

%--------------------------------------------------------------------------
%- Modality
Atom.Modality='fMRI';

%- Steps
%-> Atom.Steps={'Dcm2Nii','RemoveTimePoints','SliceTiming','Realign','T1Normalize','Smooth',''};
Atom.Steps={'Dcm2Nii','RemoveTimePoints','SliceTiming','Realign','T1Normalize','Detrend','CovRegress','Smooth','Filter'};

%- I/O
Atom.InputDir.Functionals={'E:\Rest'};
Atom.InputDir.Anatomicals={'E:\3DRaw'};
Atom.OutputDir.Functionals={'E:\Atom\Rest'};
Atom.OutputDir.Anatomicals={'E:\Atom\3D'};
% filter for files & directorys
% Atom.dirFilter.Anatomicals={'*.dcm'}; % if not advanced user, not to change !
Atom.fileFilter.Anatomicals={'*.dcm'};
% Atom.dirFilter.Functionals={'*'}; % if not advanced user, not to change !
Atom.fileFilter.Functionals={'*.dcm'};

%- Remove Time Points
Atom.RemoveTimePoints.nFirstImgDiscard=5;

%- SliceTiming
Atom.SliceTiming.SliceNumber=32;
Atom.SliceTiming.TR=2;
Atom.SliceTiming.SliceOrder=[2:2:32,1:2:31];
Atom.SliceTiming.ReferenceSlice=32;

%- Realign
Atom.Realign.Registaer2Mean=0;

%- Smooth
Atom.Smooth.FWHM=[6,6,6];


%- Covariables
Atom.Covariables.Names              ={'WM','CSF','Motion'};
Atom.Covariables.Dim                =[5,   5,    6]; % 0 for mean, >1 for PCA dimensions
Atom.Covariables.Derivative         =[0,   0,    1]; % derivative   
Atom.Covariables.Thresh             =[0.99,0.99, 1];
Atom.Covariables.Erode              =[1,   0,    1]; % for CSF usually, erode=0
Atom.Covariables.Path               ={};

%- Filter
Atom.Filter.Type='bandpass'; % 'bandpass'/'butter'
Atom.Filter.FrequencyBand=[0.01,0.1];

%- Derend
Atom.Detrend.PolyTrend=1;

%--------------------------------------------------------------------------
%- Final Run
Atom=Atom_Run(Atom);


