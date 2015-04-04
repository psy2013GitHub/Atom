
function fileList=subList2Path(subList,PDir,prefix,suffix,method,maxNum,varargin)

% method==1
% */rest_control_012_....nii
% method==2
% */rest_control_012/....nii

if exist(PDir,'dir')~=7
    disp(['Error: ',PDir,' doesnt exist, Check']);
end

if nargin>=5
    txtfile=varargin{1};
    fid_mode='w';
    txt_open=0;
elseif nargin>=6
    fid_mode=varargin{2};
    txt_open=0;
elseif nargin>=7
    txt_open=varargin{3};
end

if exist('txtfile','var')
    if (isempty(strfind(fid_mode,'a'))||isempty(strfind(fid_mod,'A')))&&exist(txtfile,'file')==2, delete(txtfile); end
    fid=fopen(txtfile,fid_mode);
end

fileList={};
if method==1
    maxLen=length([PDir,filesep,prefix,suffix])+maxNum;
elseif method==2
    maxLen=length([PDir,filesep,prefix,filesep,suffix])+maxNum;
end
for ss=1:length(subList)
    tmp_ss=subList(ss);
    if maxNum==0
        if method==1
            fname=[PDir,filesep,prefix,num2str(tmp_ss)];
        elseif method==2
            fname=[PDir,filesep,prefix,num2str(tmp_ss),filesep];
        end
    else
        lpadzero=repmat('0',[1,maxNum-1-floor(log10(tmp_ss))]);
        if method==1
            fname=[PDir,filesep,prefix,lpadzero,num2str(tmp_ss)];
        elseif method==2
            fname=[PDir,filesep,prefix,lpadzero,num2str(tmp_ss),filesep];
        end
    end
    if strfind(suffix,'*') % wilcard
        filelist=dir([fname,suffix]);
        if isempty(filelist)
            disp(['Error: ','no files satisfy ',fname,suffix]);
            continue;
        else
            tmp_suffix=filelist(1).name;
        end
    else
        tmp_suffix=suffix;
    end
    fname=[fname,tmp_suffix];
    if exist(fname,'file')~=2
        disp([fname,' doesnt exist, please Check']);
        continue;
    end
    if exist('fid','var')
        fprintf(fid,sprintf('%s\n',fname));
    end
    tmpLen=length(fname);
    rpadblank=repmat(' ',[1,maxLen-tmpLen]);
    fname=[fname,rpadblank];
    if exist(fname,'file')==2
        fileList=[fileList;fname];
    end
end

fclose(fid);
if exist('txt_open','var')&&txt_open==1
    open(txtfile);
end

return
end