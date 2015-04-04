
function dzLoopDirectory(varargin)
%-------------------------------------------------------------------------%
% Input:
%       SubjectDir, E-Prime RESULT direcitory
%       OutputDir,
%       Attribute,  defined in E-Prime List
%     *Optional: *
%       minRT,      mininum RT
%       maxRT,      maxinum RT
%       nstd,       std threshold

% Usage: dzLoopSubjects(SubjectDir,OutputDir,Attribute,...)


%-------------------------------------------------------------------------%
% Zhou Deng, Southwest University, Chongqing, China
% superdengzhou@gmail.com
%-------------------------------------------------------------------------%

if nargin<4
    error('Usage: dzLoopDirectory(SubjectDir,OutputDir,Attribute,Prefix,...)');
end
if nargin>=4,
    SubjectDir=varargin{1};
    OutputDir =varargin{2};
    Attribute =varargin{3};
    Prefix =varargin{4};
end
if nargin>=5, minRT    =varargin{5}; end
if nargin>=6, maxRT    =varargin{6}; end
if nargin>=7, nstd     =varargin{7}; end
if nargin>=8, AccExlude=varargin{8}; end
if nargin>=9, Options  =varargin{9}; end

if ~exist('minRT','var'), minRT=[];  end
if ~exist('maxRT','var'), maxRT=[];  end
if ~exist('nstd','var'),  nstd =[];  end
if ~exist('AccExlude','var'), AccExlude=0; end

if exist(SubjectDir,'dir')~=7, error('%s doesnt exist',SubjectDir); end
SubjectList=dir(SubjectDir); SubjectList=SubjectList(3:end);

diary([OutputDir,filesep,datestr(now,'dd-mmm-yyyy'),'.log']);
Session=[]; nSubject=0; CondAllSub=[];
for s=1:length(SubjectList)
    fname=SubjectList(s).name;
    [p,f,e]=fileparts(fname);
    if ~strcmpi(e,'.txt'), continue; end
    %/* optional
    LR=strfind(fname,'-');
    subID=str2num(f(LR(1)+1:LR(2)-1));
    ses=str2num(f(LR(2)+1:end));
    fprintf('\n===========================================================================\n');
    fprintf('%d——%d\n\t',subID,ses);
    %*/
    [Cond1,nTrial_raw1,meanRT1,ACC1]=dzSingleTxtRead([SubjectDir,filesep,fname],Attribute,Prefix,minRT,maxRT,nstd,AccExlude);
    nSubject=nSubject+1;
    CondAllSub=[CondAllSub;{Cond1}];
    if nSubject==1,
        Condition=Cond1(:);
    else
        Condition=unique([Condition(:);Cond1(:)]);
    end
    nCond1=length(Cond1);
    if ~any(Session==ses)
        Session=[Session,ses];
        eval(sprintf('Session_%d=[];',ses));
        for c=1:nCond1
            cond=Cond1{c};
            eval(sprintf('Session_%d_%s_raw=[];',ses,cond));
        end
    end
    eval(sprintf('Session_%d=[Session_%d;%d];',ses,ses,subID));
    for c=1:nCond1
        cond=Cond1{c};
        eval(sprintf('Session_%d_%s_raw=[Session_%d_%s_raw;%d,%f,%f,%d];',ses,cond,ses,cond,subID,meanRT1(c),ACC1(c),nTrial_raw1(c)));
    end
end

nCond=length(Condition);
Session=sort(Session);
for s=1:length(Session)
    
    ses=Session(s);
    eval(sprintf('SubjID=Session_%d;',ses));
    nSubject=length(SubjID);
    for c=1:nCond
        cond=Condition{c};
        eval(sprintf('Session_%d_%s=zeros(nSubject,size(Session_%d_%s_raw,2)-1);',ses,cond,ses,cond));
        eval(sprintf('[Uselsess,subIdx]=intersect(SubjID,Session_%d_%s_raw(:,1));',ses,cond));
        eval(sprintf('Session_%d_%s(subIdx,:)=Session_%d_%s_raw(:,2:end);',ses,cond,ses,cond));
    end
    [useless,rawOrder]=sort(SubjID);
    clear sesDat useless;
    for c=1:nCond
        cond=Condition{c};
        eval(sprintf('Session_%d=[Session_%d,Session_%d_%s];',ses,ses,ses,cond));
    end
    eval(sprintf('Session_%d=Session_%d(rawOrder,:);',ses,ses));
    RTOrder=[2:3:3*nCond-1]; ACCOrder=[3:3:3*nCond]; nTrialOrder=[4:3:3*nCond+1];
    eval(sprintf('Session_%d=Session_%d(:,[1,RTOrder,ACCOrder,nTrialOrder]);',ses,ses));
    % ! delimiter='\t';
    fid=fopen([OutputDir,filesep,'Session_',num2str(ses),'.txt'],'w');
    fprintf(fid,'SubjectID\t');
    for c=1:nCond
        cond=Condition{c};
        fprintf(fid,'%s.RT\t',cond);
    end
    for c=1:nCond
        cond=Condition{c};
        fprintf(fid,'%s.ACC\t',cond);
    end
    for c=1:nCond
        cond=Condition{c};
        fprintf(fid,'%s.nTrial\t',cond);
    end
    
    % RT
    NewOptions=dzOperations(fid,Condition,Options,'RT','(:,1)');
    eval(sprintf('Session_%d_OperRT =[];',ses));
    if NewOptions.Operation_flag.RT,
        for o=1:length(NewOptions.Operations.RT)
            NewOptions.Operations.RT{o}=strrep(NewOptions.Operations.RT{o},'%d',num2str(ses));
            resultOper=eval(NewOptions.Operations.RT{o});
            eval(sprintf('Session_%d_OperRT =[Session_%d_OperRT,resultOper];',ses,ses));
        end
    end
    
    % ACC
    NewOptions=dzOperations(fid,Condition,Options,'ACC','(:,2)');
    eval(sprintf('Session_%d_OperACC =[];',ses));
    if NewOptions.Operation_flag.ACC,
        for o=1:length(NewOptions.Operations.ACC)
            NewOptions.Operations.ACC{o}=strrep(NewOptions.Operations.ACC{o},'%d',num2str(ses));
            resultOper=eval(NewOptions.Operations.ACC{o});
            eval(sprintf('Session_%d_OperACC =[Session_%d_OperACC,resultOper];',ses,ses));
        end
    end
    
    % Write TxT
    fprintf(fid,'\n'); fclose(fid);
    dlmwrite([OutputDir,filesep,'Session_',num2str(ses),'.txt'],...
        eval(sprintf('[Session_%d,Session_%d_OperRT,Session_%d_OperACC];',ses,ses,ses)),...
        'delimiter','\t','-append');
    
end
diary('off');

return
end


function NewOptions=dzOperations(fid,Condition,RawOptions,field,suffix)

nCond=length(Condition);
prefix='Session_%d_'; len_prefix=length(prefix);
len_suffix=length(suffix);

eval(sprintf('NewOptions.Operation_flag.%s=0;',field));

fullfield=['RawOptions.Operations.',field];
if exist('RawOptions','var')&&isfield(RawOptions,'Operations')&&isfield(RawOptions.Operations,field)&&~isempty(eval(fullfield))
    if ischar(eval(fullfield))
        eval(sprintf('%s=cellstr(%s);',fullfield,fullfield));
    end
    if iscell(eval(fullfield))
        eval(sprintf('NewOptions.Operation_flag.%s=1;',field));
        Operations=eval(fullfield);
        nOper=length(Operations);
    else
        eval(sprintf('NewOptions.Operation_flag.%s=0;',field));
        fprintf('Unkonwn');
    end
end
if eval(sprintf('NewOptions.Operation_flag.%s',field))
    eval(sprintf('NewOptions.Operations.%s=cell(nOper,1);',field));
    for o=1:nOper
        oper=Operations{o};
        fprintf(fid,'(%s).%s\t',oper,field);
        try
            for c=1:nCond
                cond=Condition{c};
                start=strfind(oper,cond);
                cumlp=0;
                for ii=1:length(start)
                    p=start(ii)+cumlp;
                    oper=[oper(1:p-1),prefix,oper(p:p+length(cond)-1),suffix,oper(p+length(cond):end)];
                    cumlp=cumlp+len_prefix+len_suffix;
                end
            end
            eval(sprintf('NewOptions.Operations.%s{%d}=oper;',field,o));
        catch
            [errmag,errID]=lasterr;
            fprintf('Warning: wrong with %s operation: %s, %s',field,oper,errmag);
        end
    end
end

return
end