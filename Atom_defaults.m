

function Atom_defaults

global DEFAULTS;


%- All Supported Steps
DEFAULTS.SupportedSteps.FullNames ={'Dcm2Nii','SliceTiming','Realign','EPINormalize','T1Normalize','Smooth','CovRegress','Filter','Detrend'};
DEFAULTS.SupportedSteps.ShortNames={'',       'A',          'R',      'W',           'W',          'S',     'Cr',        'Fl',    'Dt'     };
DEFAULTS.SupportedSteps.Prefix    ={'',       'a',          'r',      'w',           'w',          's',     'c',         'f',     'd'      };

%- DeleteImediateFiles=1;
DEFAULTS.DeleteImediateFiles=0;

%- File Arrange Style, DPARSF or SPM ?
DEFAULTS.FileArrangeStyle=1; % 0 for SPM, 1 for DPARSF

%- nBlock
DEFAULTS.nBlock=10;

%- DirName for Dcm2Nii & MotionParameter
DEFAULTS.NiftiFunDirName='FunImg';
DEFAULTS.NiftiAnaDirName='AnaImg';
DEFAULTS.MotionDirName='MotionParameter';

%- Dcm2Nii default format
DEFAULTS.Dcm2NiiFormat='.nii';

%- TR (used in all frequency analysis, so list here)
DEFAULTS.TR=2;

%- Remove First Time Points
%------------------------------------------------------------------------%
DEFAULTS.RemoveTimePoints.nFirstImgDiscard=0;


%- SliceTiming
%------------------------------------------------------------------------%
DEFAULTS.SliceTiming.SliceNumber=32;
DEFAULTS.SliceTiming.SliceOrder=[2:2:32,1:2:31];
DEFAULTS.SliceTiming.ReferenceSlice=32;
DEFAULTS.SliceTiming.TR=DEFAULTS.TR;

%- Realign
%------------------------------------------------------------------------%
DEFAULTS.Realign.Register2Mean=0;


%- Normalize
%------------------------------------------------------------------------%
DEFAULTS.Normalize.BoundingBox=[-90 -126 -72;90 90 108];
DEFAULTS.Normalize.VoxSize=[3,3,3];
%---T1 Normalize---%
DEFAULTS.T1Normalize.Segment.Method='oldSegment';  % spm12 supplys new segment
DEFAULTS.T1Normalize.Segment.Output.GM   =[0,1,1]; % Native + Modulated + Unmodulated
DEFAULTS.T1Normalize.Segment.Output.WM   =[0,1,1];
DEFAULTS.T1Normalize.Segment.Output.CSF  =[0,1,1];
DEFAULTS.T1Normalize.Segment.Output.CleanUp =1; % Light Clean
DEFAULTS.T1Normalize.Segment.AffineRegularisationInSegmentation='eastern'; % /mni/eastern

%- Smooth
%------------------------------------------------------------------------%
DEFAULTS.Smooth.FWHM=[6,6,6];

%- Regress out Covariables
DEFAULTS.PolyTrend=1;

%- Filter
DEFAULTS.FrequencyBand=[0.01,0.1];


end