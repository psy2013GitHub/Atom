

CompTimeSeriesFile='\\QIU&ZHANG\dengzhou\bigsample\rumination\ica_25_spm\Component\TimeSeries'; % directory of file
BehaveScore=rumination;
FrequencyBand=[0.01 0.1];
TR=2;
parCov=Demo;
pNoise=[{WM},{CSF}];

%-------------------------------------------------------------------------%
if exist(CompTimeSeriesFile,'dir')~=7&&exist(CompTimeSeriesFile,'file')~=2, error('%s doesnt exist, please'); end

% CompList Initialization
if ischar(CompTimeSeriesFile)
    if exist(CompTimeSeriesFile,'dir')==7
        CompTimeSeriesDir=CompTimeSeriesFile;
        CompTimeSeriesList=dir(CompTimeSeriesDir);
        CompTimeSeriesList=CompTimeSeriesList(3:end);
    elseif exist(CompTimeSeriesFile,'file')==2
        [CompTimeSeriesDir,f,e]=fileparts(CompTimeSeriesFile);
        CompTimeSeriesList.name=f;
    else
        error('Invalid %s\n',CompTimeSeriesFile);
    end
elseif iscell(CompTimeSeriesFile)
    ncomp=length(CompTimeSeriesFile);
    for cc=1:ncomp
        CompTimeSeriesList(cc).name=CompTimeSeriesFile{cc};
    end
end

ncomp=length(CompTimeSeriesList);


% ALFF
for cc=1:ncomp
    [p,f,e]=fileparts(CompTimeSeriesList(cc).name);
    Comp=load([CompTimeSeriesDir,filesep,CompTimeSeriesList(cc).name]);
    fdnames=fieldnames(Comp); Comp=Comp.(fdnames{1});
    % Remove Covariables
    for subj=1:size(Comp,1)
        pnoise=[];
        for ii=1:length(pNoise)
            pnoise=[pnoise,pNoise{ii}(subj,:)'];
        end
        [b,Comp(subj,:)]=y_regress_ss(Comp(subj,:)',pnoise);
    end
    % ALFF and Bandpass Filter
    [ALFF,mALFF]=TimeSeriesALFF(Comp',TR,FrequencyBand(1),FrequencyBand(2));
    % Remove Covariables
    [b,mALFF_r]=y_regress_ss(mALFF',parCov);
    % Correlation
    [r,p]=corr(BehaveScore,mALFF_r,'type','Spearman');
    fprintf('%s: r=%.3f, p=%.3f\n',f,r,p);
end