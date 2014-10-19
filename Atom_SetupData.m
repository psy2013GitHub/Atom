

function LadyFirst=LadyFirst_SetupData(LadyFirst)


% if ~isfield(LadyFirst,field), fprintf('Error: %s is not a field of LadyFirst\n',field); return; end

%- InputDir/OutputDir
if ~isfield(LadyFirst,'InputDir') ||isempty(LadyFirst.InputDir),   fprintf('Error: Check LadyFirst.InputDir\n');  return; end
if ~isfield(LadyFirst,'OutputDir')||isempty(LadyFirst.OutputDir),  fprintf('Error: Check LadyFirst.OutputDir\n'); return; end
%- LadyFirst.OutputDir, if output to one folder
if LadyFirst.fMRI_flag&&size(LadyFirst.OutputDir.Functionals,1)>1; fprintf('Error: Multiple Output Directory in LadyFirst.OutputDir.Functionals\n'); return; end;
if LadyFirst.sMRI_flag&&size(LadyFirst.OutputDir.Anatomicals,1)>1; fprintf('Error: Multiple Output Directory in LadyFirst.OutputDir.Anatomicals\n'); return; end;
if LadyFirst.dMRI_flag&&size(LadyFirst.OutputDir.Dtis,1)>1;        fprintf('Error: Multiple Output Directory in LadyFirst.OutputDir.Dtis\n');        return; end;


% Supported Input  Style
% 1) One Parent Directory               | {'/home/allsub'}
% 2) One Directory for each Subj        | {'/subj1'; '/subj2'; '/subj3'; ...}
% Supported Outpur Style
% 1) One Parent Directory               | {'/home/allsub'}
% 2) One Directory for each Subj        | {'/subj1'; '/subj2'; '/subj3'; ...}

InputDir=LadyFirst.InputDir; OutputDir=LadyFirst.OutputDir;
AllInputFunInAFolder=1;  if LadyFirst.fMRI_flag&&length(InputDir.Functionals)>1,  AllInputFunInAFolder=0; end
AllInputAnaInAFolder=1;  if LadyFirst.sMRI_flag&&length(InputDir.Anatomicals)>1,  AllInputAnaInAFolder=0; end
AllInputDtiInAFolder=1;  if LadyFirst.dMRI_flag&&length(InputDir.Dtis)>1,         AllInputDtiInAFolder=0; end

nSubjFun=0;nSubjAna=0;nSubjDti=0;
%- LadyFirst.Functionals/Anatomical/Dtis
Flag         =[LadyFirst.fMRI_flag, LadyFirst.sMRI_flag, LadyFirst.dMRI_flag ];
Dirs         ={'Functionals',       'Anatomicals',       'Dtis'              };
nSubj        =[nSubjFun,            nSubjAna,            nSubjDti            ];
SubjID       =[{},                  {},                  {}];
AInputFolder =[AllInputFunInAFolder,AllInputAnaInAFolder,AllInputDtiInAFolder];
Idx=find(Flag==1); % valid modality

%- LadyFirst.nSub
LadyFirst.nSub=max(nSubj);
%- LadyFirst.SubID & LadyFirst.OutputDir.Functionals/Anatomicals

if LadyFirst.fMRI_flag 
    try, tmpOptions.fileFilter=LadyFirst.fileFilter.Functionals; catch, tmpOptions.fileFilter={'*'}; end
    try, tmpOptions.dirFilter=LadyFirst.dirFilter.Functionals;   catch, tmpOptions.dirFilter={'*'};  end  
    [SubjID{1},nSubj(1),LadyFirst.InputDir.Functionals]=subfun1(AllInputFunInAFolder,InputDir.Functionals,tmpOptions);
end
if LadyFirst.sMRI_flag
    try, tmpOptions.fileFilter=LadyFirst.fileFilter.Anatomicals; catch, tmpOptions.fileFilter={'*'}; end
    try, tmpOptions.dirFilter=LadyFirst.dirFilter.Anatomicals;   catch, tmpOptions.dirFilter={'*'};  end
    [SubjID{2},nSubj(2),LadyFirst.InputDir.Anatomicals]=subfun1(AllInputAnaInAFolder,InputDir.Anatomicals,tmpOptions);
end
if LadyFirst.dMRI_flag
    try, tmpOptions.fileFilter=LadyFirst.fileFilter.Dtis; catch, tmpOptions.fileFilter={'*'}; end
    try, tmpOptions.dirFilter=LadyFirst.dirFilter.Dtis;   catch, tmpOptions.dirFilter={'*'};  end
   [SubjID{3},nSubj(3),LadyFirst.InputDir.Dti]=subfun1(AllInputDtiInAFolder,InputDir.Anatomicals,tmpOptions);
end

if any(diff(nSubj(Flag==1))), fprintf('Num of Subject Not Match\n'); error('Check the input directory'); end
LadyFirst.nSubj=nSubj(Idx); LadyFirst.nSubj=LadyFirst.nSubj(1);

for ii=1:length(Idx)
    for jj=1:ii-1
        tmp=cellfun(@strcmpi,SubjID(ii),SubjID(jj),'UniformOutput',false);
        if(~all(tmp{:}))
            fprintf('Subject ID not match between %s and %s\n',Dirs{ii},Dirs{jj});
            fprintf('\t%s: ',Dirs{ii});
            for kk=1:length(LadyFirst.nSubj)
                fprintf('\t\t%s',SubjID{ii}{kk});
            end
            fprintf('\n');
            fprintf('\t%s: ',Dirs{jj});
            for kk=1:length(LadyFirst.nSubj)
                fprintf('\t\t%s',SubjID{jj}{kk});
            end
            fprintf('\n');
            error('Check the input directory');
        end
    end
end
LadyFirst.SubjID=SubjID{ii};

end


function [SubjID,nSubj,inDirs]=subfun1(SingleFolderFlag,DirList,options)

SubjID={}; inDirs={};
if SingleFolderFlag,
    subDir=dir(DirList{1}); subDir=subDir(3:end); subDir=subDir([subDir(:).isdir]);
    for ii=1:length(subDir)
        tmpDir=[DirList{1},filesep,subDir(ii).name];
        tmpFileList=dzRecurDir(tmpDir,options);
        if isempty(tmpFileList), continue; end
        subjID=subDir(ii).name;
        inDirs=[inDirs;{[DirList{1},filesep,subDir(ii).name]}];
        SubjID=[SubjID,{subjID}];
    end
else
    inDirs=DirList;
    for ii=1:length(DirList)
        tmpDir=DirList{ii};
        tmpFileList=dzRecurDir(tmpDir,options);
        if isempty(tmpFileList), error('No %s files in %s',Options.fileFilter,tmpDir); end
        [useless1,subjID,useless2]=fileparts(tmpDir);
        SubjID=[SubjID,{subjID}];
    end
end
nSubj=length(SubjID);

end
