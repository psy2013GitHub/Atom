


function Move_con(root, Subj_lst, DestDir, filter, flag, varargin)

if strcmpi(flag,'dirname'), keyname=varargin{1}; end
if strcmpi(flag,'sublist'), sublist=varargin{1}; IDfilter=varargin{2}; end
if strcmpi(flag,'sublist')&&ischar(sublist), sublist=load(sublist); end


if exist(DestDir,'dir')~=7, mkdir(DestDir); end
tmp_DestDir=[DestDir, filesep, keyname];

if strcmpi(flag,'dirname')&&exist(tmp_DestDir,'dir')~=7, 
    mkdir(tmp_DestDir); 
end

empty_flag=1;
for ss=1:length(Subj_lst)
    tmp_subdir=Subj_lst(ss);
    if ~tmp_subdir.isdir, continue; end
    if strcmpi(flag,'dirname')&&isempty(strfind(lower(tmp_subdir.name), keyname)), continue; end
    if strcmpi(flag,'sublist')
        s=regexp(tmp_subdir.name,IDfilter);
        if ~any(sublist==s.ID), continue; end
    end
    
    tmp_dir=[root,filesep,tmp_subdir.name];
    display(tmp_dir)
    cons=spm_select('List',tmp_dir,filter);
    if isempty(cons), continue; end
    empty_flag=0;
    
    if strcmpi(flag,'sublist'), tmp_DestDir=DestDir; end
    fprintf('...cp %s\n',tmp_subdir.name);
    for cc=1:size(cons,1)
        src = [tmp_dir,filesep,deblank(cons(cc,:))];
        display(src)
        dest = [tmp_DestDir,filesep,tmp_subdir.name,'_',cons(cc,:)];
        display(dest)
        copyfile(src,dest,'f');
%         copyfile(frm,tmp_DestDir,'f');
%         % pause;
%         movefile([tmp_DestDir,filesep,cons(cc,:)],[tmp_DestDir,filesep,tmp_subdir.name,'_',cons(cc,:)],'f');
    end
    clear tmp_subdir;
end

return
end