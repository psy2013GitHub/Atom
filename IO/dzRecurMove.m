


function fileout=dzRecurMove(file,varargin)

% MOVE FILE IN a RECURSIVE WAY
DestDir=varargin{1}; % Destination Directory
RawDir =varargin{2};
if iscell(file),
    fileout=cell(size(file)); 
    for n1=1:numel(file),
        fileout{n1}=dzRecurMove(file{n1},DestDir,RawDir);
    end
elseif ischar(file),
    [fpath,ffile,fext]=fileparts(file); 
    sepPos=strfind(fpath,RawDir); fpathlen=length(fpath);
    newDir=[DestDir,filesep,fpath(sepPos+length(RawDir)+1:fpathlen)];
    if ~exist(newDir,'dir'), mkdir(newDir); end
    movefile(file,newDir,'f');
    fileout=[newDir,filesep,ffile,fext];
else
    fileout=file;
end