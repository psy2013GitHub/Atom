

% Init LadyFirst

% LadyFirst.Steps

% LadyFirst.RT/nslices

% LadyFirst.nSubj

% - LadyFirst.Functionals
%   {Subject}{Session}

% - SliceTiming
%             .SliceNumber
%             .SliceOrder
%             .TR
%             .TA
%             .ReferenceSlice
% - Realign
%             .
%             .
%             .
%             .
%             .

% LadyFirst.Anatomincals

function [Atom]=Atom_Run(Atom)

%---------------%
%- Environment Init
%  Path for software
[SPMPath, fileN, extn] = fileparts(which('spm.m'));
[SPMv, SPMr] = spm('Ver'); SPMv=str2num(SPMv(4:end));
ProgramPath            = fileparts(which('Atom.m'));
ImgExt='.nii'; % forced to be nii
Atom.ImgExt=ImgExt;

%  Modality
Modality=Atom.Modality;
sMRI_flag=0; fMRI_flag=0; dMRI_flag=0;
switch Modality
    case 'sMRI', sMRI_flag=1;
    case 'fMRI', fMRI_flag=1;
    case 'dMRI', dMRI_flag=1;
    otherwise,   fprintf('Error: Unrecognized Modality: %s\n',Modality); return;
end
Atom.sMRI_flag=sMRI_flag; Atom.fMRI_flag=fMRI_flag; Atom.dMRI_flag=dMRI_flag;
%  Output

%- Parameter Init
fprintf('\n>>>>>>>>>>>>\t Init...>>>>>>>>>>\n\n');Atom=Atom_init(Atom);


%----------------------------PreProcess-----------------------------------%
nAlreadySteps=0; DirPrefix=cell(3,1); ImgPrefix=cell(3,1); % Fun - T1 -
iii=1;
while iii<=length(Atom.Steps)
    tmpStep=lower(Atom.Steps{iii}); post_preprocessFlag=0;
    %---------------%
    %- dcm2nii / to speed, directory only
    if ~isempty(strfind(tmpStep,'dcm2nii'))
        % Routine Check
        fprintf('\n>>>>>>>>>>>>\t Dcm2Nii...>>>>>>>>>>\n\n');
        if ~nAlreadySteps,
            [status,Atom]=dzIOCheck(Atom,'Dcm2Nii');
            if ~status, return; end;
            try, tmp=Atom.fileFilter.Functionals; catch, Atom.fileFilter.Functionals={'*.dcm'}; end
            try, tmp=Atom.fileFilter.Anatomicals; catch, Atom.fileFilter.Anatomicals={'*.dcm'}; end
            Atom=Atom_SetupData(Atom); % fix me
        end

        %   % Dcm2Nii
        Options.Dcm2NiiFormat=Atom.Dcm2NiiFormat;
        if Atom.sMRI_flag, % fix me infuture for multiple anatomical files
            Atom.Anatomicals=cell(Atom.nSubj);
            Options.fun.names={'dzDcm2Nii'}; Options.fun.methods={'firstfile'};
            Atom.Anatomicals=cell(Atom.nSubj,1);
            for s=1:Atom.nSubj
                Options.OutAna=[Atom.OutputDir.Anatomicals{s},filesep,Atom.NiftiAnaDirName,filesep,Atom.SubjID{s}];
                Options.tmpDir=Atom.InputDir.Anatomicals{s};
                if exist(Options.OutAna,'dir')==7, rmdir(Options.OutAna,'s'); end; mkdir(Options.OutAna); % fix me, add an option in 'Atom_defaults.m'
                Options.cmd{1}=['[tmpp,tmpf,tmpe]=fileparts(inFile);',...
                    'tmpidx=strfind(tmpp,Options.tmpDir);',...
                    'tmpoutputdir=[Options.OutAna,filesep,tmpp(tmpidx(1)+length(Options.tmpDir)+1:end)];',...
                    'dzDcm2Nii(tmpoutputdir,inFile);',... % dcm3nii will mkdir automatically
                    ]; % specific arg: 'inFiles', see below & only the first file will be converted
                Options.method='concise'; dzRecurDir(Options.tmpDir,Options); % Dcm2Nii
                Options.method='spm'; Options.fun.names={}; Options.fun.methods={}; Options.cmd={};
                tmpFileList=dzRecurDir(Options.OutAna,Options); if ~DirWithin(Options.OutAna), tmpFileList={tmpFileList}; end % wrap only only one session
                Atom.Anatomicals{s}=tmpFileList;
            end
            DirPrefix{2}=[DirPrefix{2},Atom.NiftiAnaDirName];
        end
        if Atom.fMRI_flag
            Options.fun.names={'dzDcm2Nii'}; Options.fun.methods={'firstfile'};
            Atom.Functionals=cell(Atom.nSubj,1);
            for s=1:Atom.nSubj
                Atom.Functionals=cell(Atom.nSubj,1);
                Options.OutFun=[Atom.OutputDir.Functionals{s},filesep,Atom.NiftiFunDirName,filesep,Atom.SubjID{s}];
                Options.tmpDir=Atom.InputDir.Functionals{s};
                if exist(Options.OutFun,'dir')==7, rmdir(Options.OutFun,'s'); end; mkdir(Options.OutFun); % fix me, add an option in 'Atom_defaults.m'
                Options.cmd{1}=['[tmpp,tmpf,tmpe]=fileparts(inFile);',...
                    'tmpidx=strfind(tmpp,Options.tmpDir);',...
                    'tmpoutputdir=[Options.OutFun,filesep,tmpp(tmpidx(1)+length(Options.tmpDir)+1:end)];',...
                    'dzDcm2Nii(tmpoutputdir,inFile);',... % dcm3nii will mkdir automatically
                    ]; % specific arg: 'inFiles', see below & only the first file will be converted
                Options.method='concise'; dzRecurDir(Options.tmpDir,Options);
                Options.method='spm'; Options.fun.names={}; Options.fun.methods={}; Options.cmd={};
                tmpFileList=dzRecurDir(Options.OutFun,Options); if ~DirWithin(Options.OutFun), tmpFileList={tmpFileList}; end % wrap only one session
                Atom.Functionals{s}=tmpFileList;
            end
            DirPrefix{1}=[DirPrefix{1},Atom.NiftiFunDirName];
        end
        nAlreadySteps=nAlreadySteps+1;
    end
    
    %---------------%
    %- Remove The First Few TimePoints
    if Atom.fMRI_flag&&~isempty(strfind(tmpStep,'removetimepoints'))&&Atom.RemoveTimePoints.nFirstImgDiscard>0
        if ~nAlreadySteps,
            [status,Atom]=dzIOCheck(Atom,'Dcm2Nii');
            if ~status, return; end;
            try, tmp=Atom.fileFilter.Functionals; catch, Atom.fileFilter.Functionals={['*.',ImgExt]}; end
            try, tmp=Atom.fileFilter.Anatomicals; catch, Atom.fileFilter.Anatomicals={['*.',ImgExt]}; end
            Atom=Atom_SetupData(Atom); % fix me
        end
        fprintf('\n>>>>>>>>>>>>\t Remove TimePoints...>>>>>>>>>>\n\n');
        % Remove Time Points
        nFirstImgDiscard=Atom.RemoveTimePoints.nFirstImgDiscard;
        for subj=1:Atom.nSub
            fprintf('\t\t Subject %d\n\n',subj);
            Atom.Functionals{subj}{1}=Atom.Functionals{subj}{1}(nFirstImgDiscard+1:end); % temporally, only remove first timepoints in the first session! not delete!
        end
        nAlreadySteps=nAlreadySteps+1;
    end
    
    
    spm('defaults','fmri');
    %---------------%
    %- SliceTiming
    if Atom.fMRI_flag&&~isempty(strfind(tmpStep,'slicetiming'))
        if ~nAlreadySteps,
            [status,Atom]=dzIOCheck(Atom,'SliceTiming');
            if ~status, return; end;
            try, tmp=Atom.fileFilter.Functionals; catch, Atom.fileFilter.Functionals={['*.',ImgExt]}; end
            try, tmp=Atom.fileFilter.Anatomicals; catch, Atom.fileFilter.Anatomicals={['*.',ImgExt]}; end
            Atom=Atom_SetupData(Atom); % fix me
        end
        stepidx=strmatch('SliceTiming',Atom.SupportedSteps.FullNames);
        fprintf('\n>>>>>>>>>>>>\t Slice Timing...>>>>>>>>>>\n\n');
        clear matlabbatch;
        for s=1:Atom.nSubj
            subj=Atom.SubjID{s}; fprintf('\t\t Subject %d %s\n\n',s,subj);
            spm_jobman('initcfg');
            matlabbatch{1}.spm.temporal.st.scans=Atom.Functionals{s};
            matlabbatch{1}.spm.temporal.st.nslices=Atom.SliceTiming.SliceNumber;
            matlabbatch{1}.spm.temporal.st.tr=Atom.SliceTiming.TR;
            matlabbatch{1}.spm.temporal.st.ta=Atom.SliceTiming.TR-Atom.SliceTiming.TR/Atom.SliceTiming.SliceNumber;
            matlabbatch{1}.spm.temporal.st.so=Atom.SliceTiming.SliceOrder;
            matlabbatch{1}.spm.temporal.st.refslice=Atom.SliceTiming.ReferenceSlice;
            matlabbatch{1}.spm.temporal.st.prefix = Atom.SupportedSteps.Prefix{stepidx};
            spm_jobman('run',matlabbatch);
            Atom=subfun3(s,DirPrefix,ImgPrefix,'SliceTiming',nAlreadySteps,Atom);
        end
        ImgPrefix{1}=[Atom.SupportedSteps.Prefix{stepidx},ImgPrefix{1}]; DirPrefix{1}=[DirPrefix{1},Atom.SupportedSteps.ShortNames{stepidx}];
        %- remove to SliceTming
        prefix1=ImgPrefix{1}; % for normalize use fileout=conn_prepend(str1,file,varargin)
        nAlreadySteps=nAlreadySteps+1;
    end
    
    
    %---------------%
    %- Realign
    if Atom.fMRI_flag&&~isempty(strfind(tmpStep,'realign'))
        if ~nAlreadySteps,
            [status,Atom]=dzIOCheck(Atom,'Realign');
            if ~status, return; end;
            try, tmp=Atom.fileFilter.Functionals; catch, Atom.fileFilter.Functionals={['*.',ImgExt]}; end
            try, tmp=Atom.fileFilter.Anatomicals; catch, Atom.fileFilter.Anatomicals={['*.',ImgExt]}; end
            Atom=Atom_SetupData(Atom); % fix me
        end
        stepidx=strmatch('Realign',Atom.SupportedSteps.FullNames);
        fprintf('\n>>>>>>>>>>>>\t Realign...>>>>>>>>>>\n\n');
        clear matlabbatch;
        for s=1:Atom.nSubj
            Subj=Atom.SubjID{s};
            fprintf('\t\t Subject %d %s\n\n',s,Subj);
            spm_jobman('initcfg');
            matlabbatch{1}.spm.spatial.realign.estwrite.data =Atom.xFunctionals{s};
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = Atom.Realign.Register2Mean; % register to the first volume
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = Atom.SupportedSteps.Prefix{stepidx};
            spm_jobman('run',matlabbatch);
            Atom=subfun3(s,DirPrefix,ImgPrefix,'Realign',nAlreadySteps,Atom);
            if Atom.DeleteImediateFiles, rawFiles=Atom.Functionals{s} ;subfun2(ImgExt,ImgPrefix,rawFiles); end
        end
        ImgPrefix{1}=[Atom.SupportedSteps.Prefix{stepidx},ImgPrefix{1}]; DirPrefix{1}=[DirPrefix{1},Atom.SupportedSteps.ShortNames{stepidx}];
        nAlreadySteps=nAlreadySteps+1;
    end
    
    %---------------%
    %- Normalize
    if Atom.fMRI_flag&&~isempty(strfind(tmpStep,'normalize'))
        if ~nAlreadySteps,
            [status,Atom]=dzIOCheck(Atom,'Normalize');
            if ~status, return; end;
            try, tmp=Atom.fileFilter.Functionals; catch, Atom.fileFilter.Functionals={['*.',ImgExt]}; end
            try, tmp=Atom.fileFilter.Anatomicals; catch, Atom.fileFilter.Anatomicals={['*.',ImgExt]}; end
            Atom=Atom_SetupData(Atom); % fix me
        end
    end
    if Atom.fMRI_flag&&~isempty(strfind(tmpStep,'epinormalize'))     % EPINormalize
        clear matlabbatch;
        stepidx=strmatch('EPINormalize',Atom.SupportedSteps.FullNames);
        for s=1:Atom.nSubj
            Subj=Atom.SubjID{s};
            fprintf('\t\t Subject %d %s\n\n',subj,Subj);
            load([ProgramPath,filesep,'Jobmats',filesep,'Normalize.mat']);
            jobs{1,1}.spatial{1,1}.normalise{1,1}.estwrite.subj(1,1).source=Atom.meanFunctionals(s); % cell demanded ()
            jobs{1,1}.spatial{1,1}.normalise{1,1}.estwrite.subj(1,1).resample=Atom.xFunctionals{s};
            jobs{1,1}.spatial{1,1}.normalise{1,1}.estwrite.eoptions.template={[SPMPath,filesep,'templates',filesep,'EPI.nii,1']};
            jobs{1,1}.spatial{1,1}.normalise{1,1}.estwrite.roptions.bb=Atom.Normalize.BoundingBox;
            jobs{1,1}.spatial{1,1}.normalise{1,1}.estwrite.roptions.vox=Atom.Normalize.VoxSize.func; % fix me
            spmJobRun(jobs);
            Atom=subfun3(s,DirPrefix,ImgPrefix,'EPINormalize',nAlreadySteps,Atom);
            if Atom.DeleteImediateFiles, rawFiles=Atom.xFunctionals{s}; subfun2(ImgExt,ImgPrefix,rawFiles); end % fix me
        end
        ImgPrefix{1}=[Atom.SupportedSteps.Prefix{stepidx},ImgPrefix{1}]; DirPrefix{1}=[DirPrefix{1},Atom.SupportedSteps.ShortNames{stepidx}];
        nAlreadySteps=nAlreadySteps+1;
    elseif Atom.fMRI_flag&&~isempty(strfind(tmpStep,'t1normalize'))  % T1Normalize
        if ~Atom.sMRI_flag, fprintf('Error: No T1 Available for T1Normaliz\n'); return; end
        stepidx=strmatch('T1Normalize',Atom.SupportedSteps.FullNames);
        fprintf('\n>>>>>>>>>>>>\t Normalize Using Unified-Segmentation...>>>>>>>>>>\n\n');
        for s=1:Atom.nSubj
            Subj=Atom.SubjID{s};
            fprintf('\t\t Subject %d %s\n\n',s,Subj);
            clear matlabbatch;
            % Session by Session
            AnaImg=Atom.Anatomicals{s}{1}{2}; % select the cropted one
            RealignMean=Atom.Motion.meanImages{s};
            nSes=length(Atom.xFunctionals{s});
            for ses=1:nSes
                spm_jobman('initcfg');
                % Coregister
                RealignMean=RealignMean{ses};
                matlabbatch{1}.spm.spatial.coreg.estimate.ref = RealignMean(1,:); % fix me
                matlabbatch{1}.spm.spatial.coreg.estimate.source = {AnaImg};
                matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm =[7 7];
                % Segment
                if strcmpi(Atom.T1Normalize.Segment.Method,'oldsegment')
                    matlabbatch{end+1}.spm.spatial.preproc.data = {AnaImg};
                    matlabbatch{end}.spm.spatial.preproc.output.GM  = Atom.T1Normalize.Segment.Output.GM;
                    matlabbatch{end}.spm.spatial.preproc.output.WM  = Atom.T1Normalize.Segment.Output.WM;
                    matlabbatch{end}.spm.spatial.preproc.output.CSF = Atom.T1Normalize.Segment.Output.CSF;
                    matlabbatch{end}.spm.spatial.preproc.output.biascor = 1;
                    matlabbatch{end}.spm.spatial.preproc.output.cleanup = 0;
                    if SPMv==8 % version control
                        matlabbatch{end}.spm.spatial.preproc.opts.tpm = {[SPMPath,filesep,'tpm',filesep,'grey.nii'];[SPMPath,filesep,'tpm',filesep,'white.nii'];[SPMPath,filesep,'tpm',filesep,'csf.nii']};
                    elseif SPMv==12
                        oldPath=[SPMPath,filesep,'toolbox',filesep,'OldSeg'];
                        matlabbatch{end}.spm.spatial.preproc.opts.tpm = {[oldPath,filesep,'grey.nii'];[oldPath,filesep,'white.nii'];[oldPath,filesep,'csf.nii']};
                    end
                    matlabbatch{end}.spm.spatial.preproc.opts.ngaus = [2;2;2;4];
                    matlabbatch{end}.spm.spatial.preproc.opts.regtype = Atom.T1Normalize.Segment.AffineRegularisationInSegmentation; % or eastern
                    matlabbatch{end}.spm.spatial.preproc.opts.warpreg = 1;
                    matlabbatch{end}.spm.spatial.preproc.opts.warpco = 25;
                    matlabbatch{end}.spm.spatial.preproc.opts.biasreg = 0.0001;
                    matlabbatch{end}.spm.spatial.preproc.opts.biasfwhm = 60;
                    matlabbatch{end}.spm.spatial.preproc.opts.samp = 3;
                    matlabbatch{end}.spm.spatial.preproc.opts.msk = {''};
                else
                    dzdummy(); % fix me
                end
                % Normalize_Write
                if strcmpi(Atom.T1Normalize.Segment.Method,'oldsegment')
                    matlabbatch{end+1}.spm.spatial.normalise.write.subj.matname = {[AnaImg(1:end-4),'_seg_sn.mat']};
                else
                    dzdummy(); % fix me
                end
                matlabbatch{end}.spm.spatial.normalise.write.subj.resample  = Atom.xFunctionals{s}{ses}; % no need to wrap here
                matlabbatch{end}.spm.spatial.normalise.write.roptions.preserve = 0;
                matlabbatch{end}.spm.spatial.normalise.write.roptions.bb = Atom.Normalize.BoundingBox;
                matlabbatch{end}.spm.spatial.normalise.write.roptions.vox = Atom.Normalize.VoxSize;
                matlabbatch{end}.spm.spatial.normalise.write.roptions.interp = 1;
                matlabbatch{end}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
                matlabbatch{end}.spm.spatial.normalise.write.roptions.prefix = Atom.SupportedSteps.Prefix{stepidx};
                spm_jobman('run',matlabbatch);
                Atom=subfun3(s,DirPrefix,ImgPrefix,'T1Normalize',nAlreadySteps,Atom);
                if Atom.DeleteImediateFiles, rawFiles=Atom.xFunctionals{s}; subfun2(ImgExt,ImgPrefix,rawFiles); end
            end
        end
        ImgPrefix{1}=[Atom.SupportedSteps.Prefix{stepidx},ImgPrefix{1}]; DirPrefix{1}=[DirPrefix{1},Atom.SupportedSteps.ShortNames{stepidx}];
        ImgPrefix{2}=['m',ImgPrefix{2}]; ImgPrefix{2}=['w',ImgPrefix{2}]; % segment + normalize % fix me
        ImgPrefix{3}=['w',ImgPrefix{3}];
        nAlreadySteps=nAlreadySteps+1;
    elseif ~isempty(strfind(tmpStep,'dartelnormalize'))  % Normalize by Dartel
        clear matlabbatch;
        disp('coming soon');
    end
    
    %---------------%
    %- Smooth
    if ~isempty(strfind(tmpStep,'smooth'))     % Smooth
        if ~nAlreadySteps,
            [status,Atom]=dzIOCheck(Atom,'Smooth');
            if ~status, return; end;
            try, tmp=Atom.fileFilter.Functionals; catch, Atom.fileFilter.Functionals={['*.',ImgExt]}; end
            try, tmp=Atom.fileFilter.Anatomicals; catch, Atom.fileFilter.Anatomicals={['*.',ImgExt]}; end
        end
        stepidx=strmatch('Smooth',Atom.SupportedSteps.FullNames);
        fprintf('\n>>>>>>>>>>>>\t Smooth...>>>>>>>>>>\n\n');
        clear matlabbatch;
        for s=1:length(Atom.nSubj)
            Subj=Atom.SubjID{s};
            fprintf('\t\t Subject %d %s\n\n',s,Subj);
            spm_jobman('initcfg');
            matlabbatch{1}.spm.spatial.smooth.data=dzFileCat(Atom.xFunctionals{s});
            matlabbatch{1}.spm.spatial.smooth.fwhm=Atom.Smooth.FWHM;
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = Atom.SupportedSteps.Prefix{stepidx};
            spm_jobman('run',matlabbatch);
            Atom=subfun3(s,DirPrefix,ImgPrefix,'Smooth',nAlreadySteps,Atom);
            if Atom.DeleteImediateFiles, rawFiles=Atom.xFunctionals{s}; subfun2(ImgExt,ImgPrefix,rawFiles); end
        end
        ImgPrefix{1}=[Atom.SupportedSteps.Prefix{stepidx},ImgPrefix{1}]; DirPrefix{1}=[DirPrefix{1},Atom.SupportedSteps.ShortNames{stepidx}];
        nAlreadySteps=nAlreadySteps+1;
    end
    
    
    %---------------%
    %- Post-Preprocess
    if Atom.fMRI_flag
        % setup Data
        while iii<=length(Atom.Steps)&&~isempty(strfind(tmpStep,'covregress'))||~isempty(strfind(tmpStep,'filter'))||...
                ~isempty(strfind(tmpStep,'detrend'))
            post_preprocessFlag=post_preprocessFlag+1;
            if post_preprocessFlag==1
                try, mask=Atom.Mask; catch, mask=''; end
                if ~nAlreadySteps,
                    Atom=Atom_SetupData(Atom);
                    DirPrefix={}; ImgPrefix=cell(3,1);
                end
                % prepare WM & GM & CSF
                if ~nAlreadySteps&&~isempty(strmatch('CovRegress',Atom.Steps,'exact'))
                    Atom=subfun4(Atom);
                end
                try, Options.GMs=Atom.GMs;                           end
                try, Options.WMs=Atom.WMs;                           end
                try, Options.CSFs=Atom.CSFs;                         end
                try, Options.Motion.rpTxts=Atom.Motion.rpTxts;       end % mean image is not necessary
                Options.MotionDirName=Atom.MotionDirName;
                
                % subjID
                Options.SubjID=Atom.SubjID;
                % 4D prefix
                Options.D4Prefix='4D';
                % output
                Options.OutputDir=Atom.OutputDir;
                
                Options.fun.names={}; Options.cmd={}; Options.fun.methods={};
                Options.nBlock=Atom.nBlock;
                Options.rawImgPrefix=ImgPrefix; Options.rawDirPrefix=DirPrefix;
                % !!! Attention!
                %               specific arg: 'inFiles';
                % !!!
                %               be caution of twisted 'Options', do not overlap
            end
            
            if strcmpi('CovRegress',tmpStep)
                stepidx=strmatch('CovRegress',Atom.SupportedSteps.FullNames);
                ImgPrefix{1}=[Atom.SupportedSteps.Prefix{stepidx},ImgPrefix{1}]; DirPrefix{1}=[DirPrefix{1},Atom.SupportedSteps.ShortNames{stepidx}];
                Options.DirPrefix=DirPrefix; Options.ImgPrefix=ImgPrefix;
                Options.fun.names{end+1}  ='dzCovRegress';
                Options.fun.methods{end+1}='inherit';
                Options.Covariables=Atom.Covariables;
                try, tmp=Atom.Covariables.Path{1}; catch, fprintf('!!!warning: currently Covariables.Path will be reset to empty'); end % fix me
                Options.PolyTrend=Atom.PolyTrend;
                Options.cmd{end+1}=['for ii=1:length(Options.Covariables.Names),',...                           % Ensure that masks is ordered in subject order
                                         'tmpcov=Options.Covariables.Names{ii};',...
                                         'if strcmpi(tmpcov,''gm''),',...
                                             'Options.Covariables.Path=[Options.Covariables.Path,Options.GMs{Options.subjidx}];',...
                                         'elseif strcmpi(tmpcov,''wm''),',...
                                             'Options.Covariables.Path=[Options.Covariables.Path,Options.WMs{Options.subjidx}];',...
                                         'elseif strcmpi(tmpcov,''csf''),',...
                                             'Options.Covariables.Path=[Options.Covariables.Path,Options.CSFs{Options.subjidx}];',...
                                         'elseif strcmpi(tmpcov,''motion''),',...
                                             'tmpMotionPath=fileparts(Options.dataHdr.fname);',...
                                             'if ~isempty(Options.rawDirPrefix{1}),',...
                                                 'idx=strfind(tmpMotionPath,[filesep,Options.rawDirPrefix{1}]);',...
                                                 'if isempty(idx), error(''No %s in %s'',Options.MotionDirName,tmpMotionPath); end; idx=idx(end);',...
                                                 'tmpMotionPath=[tmpMotionPath(1:idx),Options.MotionDirName,tmpMotionPath(idx+length([filesep,Options.rawDirPrefix{1}]):end)];',...
                                                 'rp_txt=dir([tmpMotionPath,filesep,''rp_*.txt'']);',...
                                                 'if length(rp_txt)~=1, error(''%d rp_*.txt in %s'',length(rp_txt),tmpMotionPath); end;',...
                                                 'Options.Covariables.Path=[Options.Covariables.Path,{[tmpMotionPath,filesep,rp_txt.name]}];',...
                                             'else,',...
                                                 '',...
                                             'end;',...
                                         'end;',...
                                    'end;',...
                                    'AllVolume=dzCovRegress(AllVolume,Options.Covariables,Options);'
                               ];
            end
            if strcmpi('Filter',tmpStep)
                stepidx=strmatch('Filter',Atom.SupportedSteps.FullNames);
                ImgPrefix{1}=[Atom.SupportedSteps.Prefix{stepidx},ImgPrefix{1}]; DirPrefix{1}=[DirPrefix{1},Atom.SupportedSteps.ShortNames{stepidx}];
                Options.DirPrefix=DirPrefix; Options.ImgPrefix=ImgPrefix;
                Options.fun.names{end+1}  ='dzIdealFilter';
                Options.fun.methods{end+1}='inherit';
                Options.TR=Atom.TR; Options.FrequencyBand=Atom.FrequencyBand;
                Options.cmd{end+1}=['AllVolume=dzIdealFilter(AllVolume, Options.TR, Options.FrequencyBand);'];
            end
            if strcmpi('Detrend',tmpStep)
                stepidx=strmatch('Detrend',Atom.SupportedSteps.FullNames);
                ImgPrefix{1}=[Atom.SupportedSteps.Prefix{stepidx},ImgPrefix{1}]; DirPrefix{1}=[DirPrefix{1},Atom.SupportedSteps.ShortNames{stepidx}];
                Options.DirPrefix=DirPrefix; Options.ImgPrefix=ImgPrefix;
                Options.fun.names{end+1}  ='dzDetrend';
                Options.fun.methods{end+1}='inherit';
                Options.PolyTrend=Atom.Detrend.PolyTrend;
                Options.cmd{end+1}=['AllVolume=dzDetrend(AllVolume,Options);'];
            end
            iii=iii+1; if iii>length(Atom.Steps), break; end; tmpStep=lower(Atom.Steps{iii});
        end
        if post_preprocessFlag
           if ~isempty(Options.fun.names)
                % Head, read Data;
                Options.fun.names  =['dzReadAllVolume',Options.fun.names];
                Options.fun.methods=['everydirect',Options.fun.methods];
                Options.mask=mask; % if not specified, will use a default mask
                Options.cmd=[
                             [
                             '[AllVolume,Options.dataHdrAll,Options.maskHdr,maskDat]=dzReadAllVolume(inFiles{1},Options.mask);'... % 'Hdr'&'Dat' saved in Options
                             'Options.dataHdr=Options.dataHdrAll(1);',... % 'dataHdr' was used for writing 3D
                             'Options.VoxIdxInRawMatrix=find(maskDat); clear maskDat;',... % VoxIdxInRawMatrix
                             ],...
                             Options.cmd...
                            ]; % specific arg
        
                % Tail, write 4D;
                Options.fun.names=[Options.fun.names,{'dzWrite4DNifti'}];
                Options.fun.methods=[Options.fun.methods,{'inherit'}];
                Options.cmd{end+1}=[
                                   'Options.dataHdr.dt=[16,0];',... % 'float32' enough!
                                   'fprintf(''Writing *.nii....\n'');',...
                                   'nTp=size(AllVolume,1);',...
                                   'inFiles=cell(nTp,1);',... % reconfigure inFiles as cell, i.e. 'spm' style
                                   'for(ii=1:nTp),',...
                                       'tmpHdr=Options.dataHdrAll(ii);',...
                                       '[path,fname,ext]=fileparts(tmpHdr.fname);',...
                                       'idx=strfind(path,[filesep,Options.rawDirPrefix{1}]);',...
                                       'if isempty(idx), error(''No %s in %s'',Options.rawDirPrefix{1},path); end; idx=idx(end);',...
                                       'fname=[path(1:idx),Options.DirPrefix{1},path(idx+length([filesep,Options.rawDirPrefix{1}]):end),filesep,Options.ImgPrefix{1},fname,ext];',...
                                       'if(size(AllVolume,2)~=prod(Options.dataHdr.dim)),',...
                                          'if(~isfield(Options,''VoxIdxInRawMatrix'')||isempty(Options.VoxIdxInRawMatrix)),',...
                                             'error(''Please specify Voxel Idx'');',...
                                          'end;',...
                                          'Data=zeros(Options.dataHdr.dim);',...
                                          'Data(Options.VoxIdxInRawMatrix)=AllVolume(ii,:);',...
                                          'dzWrite(Data,Options.dataHdr,fname);',... % 'dzWrite4DNIfti(AllVolume,Options.dataHdr);',...
                                       'else,',...
                                          'dzWrite(reshape(AllVolume(ii,:),[Options.dataHdr.dim]),Options.dataHdr,fname);',...
                                       'end;',...
                                   'inFiles{ii}=fname;',...
                                   'end;'
                                   ];
            end
            Options.method='concise';
            Options.fileFilter={[Options.rawImgPrefix{1},'*',ImgExt]};
            % eval command
            for s=1:Atom.nSubj
                Subj=Atom.SubjID{s}; Options.subjidx=s;
                fprintf('\t\t Subject %d %s\n\n',s,Subj);
                if Atom.FileArrangeStyle&&nAlreadySteps
                    tmpDir=[Atom.OutputDir.Functionals{1},filesep,Options.rawDirPrefix{1},filesep,Atom.SubjID{s}];
                else
                    tmpDir=[Atom.InputDir.Functionals{1},Atom.SubjID{s}];
                end % always write in a new directory
                Options.DestDir=[Atom.OutputDir.Functionals{1},filesep,DirPrefix{1},filesep,Atom.SubjID{s}];
                if ~exist(Options.DestDir,'dir'), mkdir(Options.DestDir); end
                Atom.xFunctionals{s}=dzRecurDir(tmpDir,Options);
            end
        end
    end
    if ~post_preprocessFlag, iii=iii+1; end 
end

fprintf('\n');

return
end

function [status]=subfun1(DirLen,FileCell)
status=1;
if length(DirLen)==1&&~isempty(FileCell{1})&&ischar(FileCell{1})
    status=0;
end
end

function [status,Atom]=dzIOCheck(Atom,field)
status=1; Modality=Atom.Modality;
% if ~isfield(Atom,field),     fprintf('Error: %s isn''t a Field of Atom\n',field); status=0; return; end
if ~isfield(Atom,'InputDir'),fprintf('Error: InputDir isn''t a Field of Atom\n'); status=0; return; end
if Atom.fMRI_flag, if ~isfield(Atom.InputDir,'Functionals'),fprintf('Error: Functionals isn''t a Field of Atom.InputDir\n');   status=0; return; end; end
if Atom.sMRI_flag, if ~isfield(Atom.InputDir,'Anatomicals'),fprintf('Error: Anatomicals isn''t a Field of Atom.InputDir\n\n'); status=0; return; end; end
if Atom.dMRI_flag, if ~isfield(Atom.InputDir,'Anatomicals'),fprintf('Error: Anatomicals isn''t a Field of Atom.InputDir\n\n'); status=0; return; end; end
%- Modality specific
if Atom.fMRI_flag&&eval('isfield(Atom.InputDir,''Anatomicals'')&&~isempty(Atom.InputDir.Anatomicals)'), Atom.sMRI_flag=1; end
%- field specific

if ~isfield(Atom,'OutputDir'),fprintf('Error: OutputDir isn''t a Field of Atom\n'); status=0; return; end
if Atom.fMRI_flag&&~isfield(Atom.OutputDir,'Functionals'),fprintf('Error: Functionals isn''t a Field of Atom.OutputDir\n'); status=0; return; end
if Atom.sMRI_flag&&~isfield(Atom.OutputDir,'Anatomicals'),fprintf('Error: Anatomicals isn''t a Field of Atom.OutputDir\n'); status=0; return; end

end

function ndir=DirWithin(indir)

if ~isdir(indir), error('Not A Dir'); end
subdirs=dir(indir); subdirs=subdirs(3:end);
ndir=length(subdirs([subdirs(:).isdir]));

end

function subfun2(ImgExt,ImgPrefix,incell)
% delete immediate files
if strcmpi(ImgExt,'.img')||strcmpi(ImgExt,'.hdr')
    dzDeleteFile(conn_prepend(ImgPrefix{1},incell,'.*'));
else
    dzDeleteFile(conn_prepend(ImgPrefix{1},incell));
end
end

function atom=subfun3(subjidx,dirprefix,imgprefix,currstep,nalreadysteps,atom)
% --------
% !!! currently, only support 'fMRI', fix me in future
%
% --------
% remove files
if ~atom.FileArrangeStyle,
    atom.xFunctionals{subjidx}=conn_prepend(imgprefix{1},Atom.xFunctionals{subj});
    return;
end

% Step idx
rawDirPrefix=dirprefix;
stepidx=strmatch(currstep,atom.SupportedSteps.FullNames);
imgprefix{1}=[atom.SupportedSteps.Prefix{stepidx},imgprefix{1}]; % caaution, 'ra' not 'ar'
dirprefix{1}=[dirprefix{1},atom.SupportedSteps.ShortNames{stepidx}];

% 1st, move every file
Options=dzRecurDir_options('moveeveryfile'); % kept in Options, below refresh
% 2st, method & RawDir&DestDir&fileFilter
Options.method='spm'; Options.Dcm2NiiFormat=atom.Dcm2NiiFormat;

if nalreadysteps
    Options.RawDir=[atom.OutputDir.Functionals{1},filesep,rawDirPrefix{1},filesep,atom.SubjID{subjidx}];
else
    Options.RawDir=[atom.InputDir.Functionals{1}];
end
Options.DestDir=[atom.OutputDir.Functionals{1},filesep,dirprefix{1},filesep,atom.SubjID{subjidx}];
if exist(Options.DestDir,'dir')~=7, mkdir(Options.DestDir); end;
Options.fileFilter={[imgprefix{1},'*',atom.Dcm2NiiFormat]};
atom.xFunctionals{subjidx}=dzRecurDir(Options.RawDir,Options);
if ~DirWithin(Options.DestDir), atom.xFunctionals{subjidx}={atom.xFunctionals{subjidx}}; end % wrap only one session

% ++++++ Step specific
if strcmpi(currstep,'Realign') % if 'Realign', mv 'rp*/mean*' as well
    Options.DestDir=[atom.OutputDir.Functionals{1},filesep,atom.MotionDirName,filesep,atom.SubjID{subjidx}];
    % rp_*.txt
    Options.Subj=atom.SubjID{subjidx}; Options.ImgExt=atom.ImgExt; Options.rawDirPrefix=rawDirPrefix;
    Options.OutputDir=atom.OutputDir;  Options.MotionDir=atom.MotionDirName;
    Options.fileFilter={'rp_*.txt'}; % 'rp' on top
    Options.fun.names       =['dzRealignParaCalc';Options.fun.names]; % put this as 'Head'
    Options.fun.methods     =['everyfile';Options.fun.methods];
    Options.cmd=[[
                 'tmpfile=inFile;',...
                 '[tmpp,tmpf,tmpe]=fileparts(tmpfile);',...
                 'tmpidx=strfind(tmpp,Options.RawDir);',... % Options.RawDir
                 'tmpoutputdir=[Options.DestDir,filesep,tmpp(tmpidx(1)+length(Options.RawDir)+1:end)];',... % Options.DestDir
                 'if exist(tmpoutputdir,''dir'')~=7, mkdir(tmpoutputdir); end; ',... % mkdir, if exist, then ignore
                 'tmpmeanImage=[tmpp,filesep,''mean'',tmpf(4:end-4),''*''];',...
                 'tmpmeanImage=dir(tmpmeanImage);',...
                 'if (strcmpi(Options.ImgExt,''.nii'')&&length(tmpmeanImage)~=1)&&(strcmpi(Options.ImgExt,''.img'')||strcmpi(Options.ImgExt,''.hdr'')&&length(tmpmeanImage)~=2), error(''In correct number of mean images''); end;',...
                 'tmpmeanImage=[tmpp,filesep,tmpmeanImage(1).name];',...
                 'tmpidx=strfind(tmpp,[Options.rawDirPrefix{1},filesep]); tmpp=[tmpp(tmpidx+length(Options.rawDirPrefix{1}):end),filesep]; tmpidx=strfind(tmpp,filesep);',... % Notice here, tmpphad been changed
                 'ExcludeInfo=[]; start=1; for ii=1:length(tmpidx), ExcludeInfo=[ExcludeInfo,tmpp(start:tmpidx(ii)-1)]; start=start+tmpidx(1); end;',...
                 'dzRealignParaCalc(tmpfile,tmpmeanImage,tmpoutputdir,[Options.OutputDir.Functionals{1},filesep,Options.MotionDir],ExcludeInfo);',...
                 ];
                 Options.cmd...
                ];% specific arg: 'inFiles', see below & only the first file will be converted
    atom.Motion.rpTxts{subjidx}=dzRecurDir(Options.RawDir,Options);     if ~DirWithin(Options.RawDir), atom.Motion.rpTxts{subjidx}=atom.Motion.rpTxts(subjidx); end
    % mean*.nii
    Options.fun.names=Options.fun.names(2); Options.fun.methods=Options.fun.methods(2); Options.cmd=Options.cmd(2); % Caution!, remove 1st arg
    Options.fileFilter={['mean*',atom.Dcm2NiiFormat]}; % 'rp' on top
    atom.Motion.meanImages{subjidx}=dzRecurDir(Options.RawDir,Options); if ~DirWithin(Options.RawDir), atom.Motion.meanImages{subjidx}=atom.Motion.meanImages(subjidx); end
elseif strcmpi(currstep,'T1Normalize') % move 'Segments' as well
    if nalreadysteps
        Options.RawDir=[atom.OutputDir.Anatomicals{1},filesep,rawDirPrefix{2},filesep,atom.SubjID{subjidx}];
    else
        Options.RawDir=[atom.InputDir.Anatomicals{1},filesep,atom.SubjIdx{subjidx}]; % fix me
    end
    Options.DestDir=[atom.OutputDir.Anatomicals{1},filesep,atom.T1Normalize.Segment.Method,filesep,atom.SubjID{subjidx}];
    Options.fileFilter={'c1*','c2*','c3*'};
    if strcmpi(atom.T1Normalize.Segment.Method,'oldsegment')
        Options.fileFilter = [Options.fileFilter,{'*_seg_sn.mat'},{'*_seg_inv_sn.mat'}];
    else
        dzdummy(); % fix me
    end
    dzRecurDir(Options.RawDir,Options);
    % normalized GM
    Options.fileFilter={'wc1*'};
    atom.GMs{subjidx}=dzRecurDir(Options.RawDir,Options);
    % normalized WM
    Options.fileFilter={'wc2*'};
    atom.WMs{subjidx}=dzRecurDir(Options.RawDir,Options);
    % normalized CSF
    Options.fileFilter={'wc3*'};
    atom.CSFs{subjidx}=dzRecurDir(Options.RawDir,Options);
end
return
end

function atom=subfun4(atom)

%--- GM || WM || CSF
phyNoiseFlag=0;
if ~isempty(strmatch('gm', lower(atom.Covariables.Names)))||~isempty(strmatch('wm', lower(atom.Covariables.Names)))||...
        ~isempty(strmatch('csf',lower(atom.Covariables.Names)))
    phyNoiseFlag=1;
    fileFilter=atom.fileFilter.Anatomicals;
    AllInputAnaInAFolder=1;  if length(atom.InputDir.Anatomicals)>1,  AllInputAnaInAFolder=0; end
    % GM
    GMFlag=0;
    if ~isempty(strmatch('gm', lower(atom.Covariables.Names))),
        GMFlag=1;
        tmp=cellfun(@strfind,fileFilter,repmat({'c1'},size(fileFilter)),'UniformOutput',false);
        if length([tmp{:}])~=1, error('error with the programe, please fix'); end
        tmp=cellfun(@isempty,tmp,'UniformOutput',false); tmp=~[tmp{:}];
        atom.GM.fileFilter=fileFilter(tmp); % in cell format
        atom.GMs=cell(atom.nSubj,1);
    end
    WMFlag=0;
    % WM
    if ~isempty(strmatch('wm', lower(atom.Covariables.Names))),
        WMFlag=1;
        tmp=cellfun(@strfind,fileFilter,repmat({'c2'},size(fileFilter)),'UniformOutput',false);
        if length([tmp{:}])~=1, error('error with the programe, please fix'); end
        tmp=cellfun(@isempty,tmp,'UniformOutput',false); tmp=~[tmp{:}];
        atom.WM.fileFilter=fileFilter(tmp); % in cell format
        atom.WMs=cell(atom.nSubj,1);
    end
    CSFFlag=0;
    % CSF
    if ~isempty(strmatch('csf',lower(atom.Covariables.Names))),
        CSFFlag=1;
        tmp=cellfun(@strfind,fileFilter,repmat({'c3'},size(fileFilter)),'UniformOutput',false);
        if length([tmp{:}])~=1, error('error with the programe, please fix'); end
        tmp=cellfun(@isempty,tmp,'UniformOutput',false); tmp=~[tmp{:}];
        atom.CSF.fileFilter=fileFilter(tmp); % in cell format
        atom.CSFs=cell(atom.nSubj,1);
    end
end

%--- Motion
motNoiseFlag=0;
if ~isempty(strmatch('motion',lower(atom.Covariables.Names))),
    motNoiseFlag=1;
    atom.Motion.rpTxt.fileFilter={'rp_*.txt'}; % fix me
    atom.Motion.meanImage.fileFilter={'mean*.*'}; % fix me
    atom.Motion.rpTxts=cell(atom.nSubj,1);
    atom.Motion.meanImags=cell(atom.nSubj,1);
end

Options=dzRecurDir_options('subjmatch'); Options.method='concise'; Options.SubjID=atom.SubjID;
for s=1:atom.nSubj % in subj order!!!
    if phyNoiseFlag
        if AllInputAnaInAFolder
            tmpAnaDir=[atom.InputDir.Anatomicals{1},filesep,atom.SubjID{s}];
        else
            tmpAnaDir=[atom.InputDir.Anatomicals{s}];
        end
        % GM
        if GMFlag,
            Options.fileFilter=atom.GM.fileFilter;
            [results,Options2]=dzRecurDir(tmpAnaDir,Options);
            if isempty(results), error('No GMs found'); end
            atom.GMs(Options2.subjidx)=results;
        end
        % WM
        if WMFlag,
            Options.fileFilter=atom.WM.fileFilter;
            [results,Options2]=dzRecurDir(tmpAnaDir,Options);
            if isempty(results), error('No WMs found'); end
            atom.WMs(Options2.subjidx)=results;
        end
        % CSF
        if CSFFlag,
            Options.fileFilter=atom.CSF.fileFilter;
            [results,Options2]=dzRecurDir(tmpAnaDir,Options);
            if isempty(results), error('No CSFs found'); end
            atom.CSFs(Options2.subjidx)=results;
        end
    end
    % Motion
    if motNoiseFlag,
        AllInputMotionInAFolder=1;  if length(atom.InputDir.Motions)>1,  AllInputMotionInAFolder=0; dzdummy(); error('sorry'); end
        tmpMotionDir=[atom.InputDir.Motions{1},filesep,atom.SubjID{s}];
        % rpTxt
        Options.fileFilter=atom.Motion.rpTxt.fileFilter;
        [results,Options2]=dzRecurDir(tmpMotionDir,Options);
        if isempty(results), error('No rp_*.txt found'); end
        atom.Motion.rpTxts(Options2.subjidx)=results;
        % meanImage
        Options.fileFilter=atom.Motion.meanImage.fileFilter;
        [results,Options2]=dzRecurDir(tmpMotionDir,Options);
        if isempty(results), error('No mean*.* found'); end
        atom.Motion.meanImages(Options2.subjidx)=results;
    end
end

end



