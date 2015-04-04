
function [Condition, nTrial_raw, meanRT, ACC] = dzSingleTxtRead(varargin)
% Input
%     fname: file path
%     minRT
%     maxRT
% Output
%     meanRT
%     ACC
%
% Usage: [Condition, meanRT, ACC]=dzSingleTxtRead(fname,attribute,minRT,maxRT,nStd)
%
%------------------------------------------------------------------------%
% Zhou Deng, Southwest University, Chongqing, China
% superdengzhou@gmail.com
%------------------------------------------------------------------------%


if nargin<3
    error('Usage: [Condition, meanRT, ACC]=dzSingleTxtRead(fname,prefix,attribute,...)');
end
if nargin>=3,
    fname=varargin{1};
    attribute=varargin{2};
    prefix=varargin{3};
end
if nargin>=4, minRT    =varargin{4}; end
if nargin>=5, maxRT    =varargin{5}; end
if nargin>=6, nStd     =varargin{6}; end
if nargin>=7, AccExlude=varargin{7}; end

if ~exist('minRT','var'),     minRT=[];    end
if ~exist('maxRT','var'),     maxRT=[];    end
if ~exist('nStd','var'),      nStd =[];    end
if ~exist('AccExlude','var'), AccExlude=0; end


% Read Txt
nTrial=0; Condition={};
fid=fopen(fname);
while(~feof(fid))
    l=fgetl(fid);
    if isempty(l), continue; end
    [start,remain]=strtok(l);
    if strcmpi(start,'***')
        [start,remain]=strtok(remain);
        % Trial Block Start
        if strcmpi(start,'LogFrame')&&~isempty(strfind(remain,'Start'))
            fprintf('.');
            nTrial=nTrial+1;
            while ~feof(fid)
                l=fgetl(fid);
                [start,remain]=strtok(l);
                switch start
                    case [attribute,':']
                        cond=strtrim(remain);
                        cond=strrep(cond,' ','_'); % replace backspace
                        cond_idx=strmatch(cond,Condition,'exact');
                        if isempty(cond_idx)
                            Condition=[Condition,{cond}];
                            % init
                            eval(sprintf('Cond_%s_RT=[];',cond));
                            eval(sprintf('Cond_%s_ACC=[];',cond));
                        end
                    case [prefix,'.RT:']
                        rt=str2num(remain);
                        eval(sprintf('Cond_%s_RT=[Cond_%s_RT;%d];',cond,cond,rt)); % int as default in E-Prime
                    case [prefix,'.ACC:']
                        acc=str2num(remain);
                        eval(sprintf('Cond_%s_ACC=[Cond_%s_ACC;%d];',cond,cond,acc)); %
                end
                if strcmpi(start,'***')
                    [start,remain]=strtok(remain);
                    if strcmpi(start,'LogFrame')&&~isempty(strfind(remain,'End'))
                        break;
                    end
                end
            end
        end
        % Trial Block End
    end
end
fprintf('\n');
fclose(fid);

%------------------------------------------------------------------------%
% Post-Statistical
Condition=sort(Condition);
nCondition=length(Condition);
meanRT=zeros(1,nCondition);
ACC=zeros(1,nCondition);

allrt=[]; nTrial_raw=[]; nTrial_origACC=[]; nTrial_minRT=[]; nTrial_maxRT=[]; nTrial_nStd=[]; nTrial_rt=[];
for c=1:nCondition
    
    cond=Condition{c};
    eval(sprintf('nTrial_raw=[nTrial_raw,length(Cond_%s_ACC)];',cond));
    eval(sprintf('rt=Cond_%s_RT;',cond));
    
    eval(sprintf('rt(Cond_%s_ACC==0)=[];',cond));
    nTrial_origACC=[nTrial_origACC;length(rt)];
    
    if ~isempty(minRT),
        eval('nTrial_minRT=[nTrial_minRT;sum(rt<minRT)];');
        rt=rt(rt>=minRT);
    end
    
    if ~isempty(maxRT),
        eval('nTrial_maxRT=[nTrial_maxRT;sum(rt>maxRT)];');
        rt=rt(rt<=maxRT);
    end
    
    nTrial_rt=[nTrial_rt;length(rt)];
    allrt=[allrt;rt(:)];
    
end
if ~isempty(nStd)
    meanrt=mean(allrt);
    stdrt=std(allrt);
end
up=1;
for c=1:nCondition
    
    cond=Condition{c};
    rt=allrt(up:up+nTrial_rt(c)-1);
    
    if ~isempty(nStd),
        if length(nStd)==1
            if nStd(1)>0
                eval('nTrial_nStd=[nTrial_nStd;sum(rt>meanrt+nStd*stdrt)];');
                rt=rt(rt<=meanrt+nStd*stdrt);
            elseif nStd(1)<0
                eval('nTrial_nStd=[nTrial_nStd;sum(rt<meanrt-nStd*stdrt)];');
                rt=rt(rt>=meanrt-nStd*stdrt);
            end
        elseif length(nStd)==2&&prod(nStd)<0
            nStd=sort(nStd);
            eval('nTrial_nStd=[nTrial_nStd;sum((rt<meanrt+nStd(1)*stdrt)|(rt>meanrt+nStd(2)*stdrt))];');
            rt=rt(rt<=meanrt+nStd(2)*stdrt&rt>=meanrt+nStd(1)*stdrt);
        end
    end
    
    meanRT(c)=mean(rt(:));
    if  AccExlude, eval(sprintf('ACC(%d)=length(rt)/nTrial_raw(c);',c,cond));        end
    if ~AccExlude, eval(sprintf('ACC(%d)=nTrial_origACC(c)/nTrial_raw(c);',c,cond)); end
    
    up=up+nTrial_rt(c);
    
end

%------------------------------------------------------------------------%
% Screen Disp
rowName=[{'Total Trials'},{'False Trials'}]; rowDat=[nTrial_raw(:)';(nTrial_raw(:)-nTrial_origACC(:))'];
if ~isempty(minRT), rowName=[rowName,{sprintf('Trials < %d',minRT)}];    rowDat=[rowDat;nTrial_minRT(:)'];   end
if ~isempty(maxRT), rowName=[rowName,{sprintf('Trials > %d',maxRT)}];    rowDat=[rowDat;nTrial_maxRT(:)'];   end
if ~isempty(nStd),
    if length(nStd)==1
        if nStd(1)>0
            rowName=[rowName,{sprintf('Trials > %d*std',nStd)}]; rowDat=[rowDat;nTrial_nStd(:)'];
        elseif nStd(1)<0
            rowName=[rowName,{sprintf('Trials < %d*std',nStd)}]; rowDat=[rowDat;nTrial_nStd(:)'];
        end
    elseif length(nStd)==2&&prod(nStd)<0
        if abs(nStd(1))==nStd(2)
            rowName=[rowName,{sprintf('Trials <|> %d*std',nStd(2))}]; rowDat=[rowDat;nTrial_nStd(:)'];
        else
            rowName=[rowName,{sprintf('Trials < %d*std | Trials < %d*std',abs(nStd(1)),nStd(2))}]; rowDat=[rowDat;nTrial_nStd(:)'];
        end
    end
end

dzTableDisp(Condition,rowName,rowDat);
fprintf('\n');

return
end

function dzTableDisp(colName,rowName,rowDat)

clen=length(colName)+1;
rlen=length(rowName)+1;

maxBox=max(max(cellfun(@length,colName),max(cellfun(@length,rowName))));

line=repmat('-',[1,maxBox*(clen)+clen]);
linlen=length(line);
for r=1:rlen
    fprintf('%s\n',line);
    for c=1:clen
        % first box
        if r==1&&c==1,    fprintf([repmat(' ',[1,maxBox]),'|']); end
        % first row
        if r==1&&c~=1
            [formerspa,latterspa]=fillspace(maxBox,length(colName{c-1}));
            fprintf([formerspa,colName{c-1},latterspa,'|']);
        end
        % first column
        if r~=1&&c==1,
            fprintf([rowName{r-1},repmat(' ',[1,maxBox-length(rowName{r-1})]),'|']);
        end
        % fill data
        if r~=1&&c~=1
            num=num2str(rowDat(r-1,c-1));
            [formerspa,latterspa]=fillspace(maxBox,length(num));
            fprintf([formerspa,num,latterspa,'|']);
        end
    end
    fprintf('\n');
end
fprintf('%s\n',line);

return
end
function [formerspa,latterspa]=fillspace(maxBox,existlen)

formerlen=floor((maxBox-existlen)/2); latterlen=maxBox-existlen-formerlen;
formerspa=repmat(' ',[1,formerlen]);  latterspa=repmat(' ',[1,latterlen]);

return
end