

function [FileList,DirTree]=dzRecurDir2(Dir,Options,DirTree,TreeDepth)

% default in spm way, i.e. each file in a cell
% if specified in consice way, each directory is a cell

if nargin<2, Options=struct(); end
if ~isfield(Options,'fileFilter')||isempty(Options.fileFilter),Options.fileFilter='*';   end
if ~isfield(Options,'dirFilter') ||isempty(Options.dirFilter), Options.dirFilter ='*';   end
if ~isfield(Options,'method')    ||isempty(Options.method),    Options.method    ='spm'; end

if nargin<3, DirTree=struct('name','','num',0); end
if nargin<4, TreeDepth=1; end
if isdir(Dir)
    files=dir([Dir,filesep,Options.fileFilter]); files=files(~logical([files(:).isdir])); nFiles=length(files);
    dirs =dir([Dir,filesep,Options.dirFilter]);  dirs=dirs(logical([dirs(:).isdir]));     nDirs =length(dirs);
    if length(dirs)==1&&strcmpi(dirs(1).name,'.')||strcmpi(dirs(1).name,'..'), nDirs=0;  end
    if length(dirs)==2&&strcmpi(dirs(1).name,'.')||strcmpi(dirs(2).name,'..'), nDirs=0;  end
    DirTree(TreeDepth).num=nDirs;
    newDir=[files;dirs];
    fileFlag=0;
    idxCell=1; % file
    FileList={}; % in case of empty output
    for ii=1:length(newDir)
        idx=ii;
        if ii>nFiles&&strcmpi(newDir(ii).name,'.')||strcmpi(newDir(ii).name,'..'), continue; end
        if newDir(ii).isdir
            DirTree(TreeDepth).name=[DirTree.name,{newDir(ii).name}];
            if idxCell>1&&ischar(FileList(1)), FileList={FileList}; end % in case FileList is char and FileList{2}=* expression is not valid
            [FileList{idxCell,1},DirTree]=dzRecurDir2([Dir,filesep,newDir(ii).name],Options,DirTree,TreeDepth+1);
            if isempty(DirTree(TreeDepth+1).name), DirTree(TreeDepth+1)=[]; end
            if ~isempty(FileList{idxCell,1}),
                idxCell=idxCell+1;
            else
                FileList(idxCell)=[]; % delete empty ones
            end
        elseif strcmpi(Options.method,'concise')   % in a concise way
            if ~fileFlag % first file flag
                iiStart=idx;
                idxCell=iiStart+1;
                if ~nDirs % no dorectory
                    FileList=[Dir,filesep,newDir(ii).name];
                else
                    FileList{iiStart,1}=[Dir,filesep,newDir(ii).name];
                end
            else
                if ~nDirs % no dorectory
                    FileList=strvcat(FileList,[Dir,filesep,newDir(ii).name]);
                else
                    FileList{iiStart,1}=strvcat(FileList{iiStart},[Dir,filesep,newDir(ii).name]);
                end
            end
            fileFlag=fileFlag+1;
        elseif strcmpi(Options.method,'spm')     % in spm way, each file in a cell
            FileList{idxCell,1}=[Dir,filesep,newDir(ii).name];
            idxCell=idxCell+1;
        end
    end
else
    FileList={Dir};
end
end