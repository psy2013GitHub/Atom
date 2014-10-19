


function dzFCseed(resultDir, EPIfiles, seed, Options)



% init
if nargin<2
    fprintf('Error: Usage, dzFC(resultDir, EPIfiles, seed...)\n');
end


RESTp=fileparts(which('rest')); mask=[RESTp,filesep,'mask',filesep,'BrainMask_05_61x73x61.img'];
CovTxt=[];
nDetrend=[];
bufFilter=[];
banFilter=[];
nBlock=10;
R2Z=1;
TR=2;
if exist('Options','var')
    if isfield(Options,'mask')&&~isempty(Options.mask);                mask       =Options.mask;      end
    if isfield(Options,'CovTxt')&&~isempty(Options.CovTxt);            CovTxt     =Options.CovTxt;    end
    if isfield(Options,'TR')&&~isempty(Options.TR);                    TR         =Options.TR;        end
    if isfield(Options,'nDetrend')&&~isempty(Options.nDetrend);        nDetrend   =Options.nDetrend;  end
    if isfield(Options,'bufFilter')&&~isempty(Options.bufFilter);      bufFilter  =Options.bufFilter; end
    if isfield(Options,'banFilter')&&~isempty(Options.banFilter);      banFilter  =Options.banFilter; end
    if isfield(Options,'nBlock')&&~isempty(Options.nBlock);            nBlock     =Options.nBlock;    end
    if isfield(Options,'R2Z')&&~isempty(Options.R2Z);                  R2Z        =Options.R2Z;       end
end
    
% mask & seed File check
if exist(mask,'file')~=2, fprintf('Error: mask dosnt exist\n'); return; end
if exist(seed,'file')~=2, fprintf('Error: seed dosnt exist\n'); return; end

% Load Covariables
nCov=length(CovTxt);
for cc=1:nCov
    cov=load(CovTxt{cc});
    if cc==1
        txtSize=size(cov);
        CovNum=zeros(txtSize(1),nCov,txtSize(2)); % [nTp, nCov, nSubj];
    else
        if ~all(size(cov)==txtSize), error('Covariable dimension not match\n'); end
    end
    CovNum(:,cc,:)=cov;
end
try, CovSize=size(CovNum); end

% maskName, maskHdr
[p,maskName]=fileparts(mask); maskHdr=spm_vol(mask);
% seedName, seedHdr
[p,seedName]=fileparts(seed);seedHdr=spm_vol(seed);


D4flag=0;
for ff=1:size(EPIfiles,1)
    
    tic;
    file=EPIfiles(ff,:);
    [p,f,e]=fileparts(file); [p,subj]=fileparts(p);
    fprintf('============================ Sub %d %s\n',ff,subj);
    fprintf('Read Data\n');
    hdr=spm_vol(file);
    if ff==1
        
        % data dimension, 3D or 4D?
        if length(hdr)>1 %4D
            D4flag=1;
            nTp=size(data,4);
        elseif length(hdr)==1 %3D
            D3flag=1;
        end
        
        dataHdr=hdr(1); % hdr(1) in case of 4D
        M=dataHdr.dim(1); N=dataHdr.dim(2); O=dataHdr.dim(3);
        
        % maskDat
        if ~all(maskHdr.mat==dataHdr.mat),
            fprintf('Warning: mask dimension not match, resize\n');
            Outfname=dzImgResize(maskHdr,dataHdr);
        else
            Outfname=mask;
        end
        maskDat=single(spm_read_vols(spm_vol(Outfname)));
        maskDat(isnan(maskDat))=0; maskDat=logical(maskDat); maskDat=maskDat(:);
        if ~all(maskHdr.mat==dataHdr.mat), delete(Outfname); end
        nVox=sum(maskDat);
        
        % seedDat
        if ~all(seedHdr.mat==dataHdr.mat),
            fprintf('Warning: seed dimension not match, resize\n');
            Outfname=dzImgResize(seedHdr,dataHdr);
        else
            Outfname=seed;
        end
        seedDat=single(spm_read_vols(spm_vol(Outfname)));
        seedDat(isnan(seedDat))=0; seedDat=logical(seedDat); seedDat=seedDat(:);
        if ~all(seedHdr.mat==dataHdr.mat), delete(Outfname); end
        
    else % hdr(1) in case of 4D
        if ~all(dataHdr.mat==hdr(1).mat), fprintf('Error: image dimension not match\n'); return; end
    end
    
    % read data
    data=spm_read_vols(hdr);
    
    if D4flag
        if length(size(data))~=4 
            fprintf('Error: image dimension not match\n'); return;
        end
        data=reshape(data,[M*N*O,nTp]); data=data';
        
        % Remove Covariable
        if exist('NumCov','var'),
            fprintf('\tRemove Covariable\n');
            [useless,data]=dzRegress(data, NumCov);           
        end
        
        % Detrend
        if ~isempty(nDetrend),
            fprintf('\tDetrend %d\n',nDetrend);
            data=spm_detrend(data,nDetrend);                   
        end
        
        % Despike
        
        
        % Filter
        if ~isempty(bufFilter),   
            [bfilter, afilter] = butter(5, bufFilter); data = filtfilt(bfilter, afilter, data); 
        end
        if ~isempty(banFilter),
            fprintf('\tBandpass Filter %s-%s\n',num2str(banFilter(1)),num2str(banFilter(2)));
            data = rest_IdealFilter(data, TR, banFilter);
        end
        
        % seed & maskDat
        seed=mean(data(:,seedDat),2);
        if std(seed)<eps, fprintf('Error, standard devation of seed is zero\n'); return; end
        data=data(:,maskDat);

        % FC Calc
        fprintf('\tFC Calc\t');
        blocksize=floor(nVox/nBlock);
        R=zeros(M*N*O,1);
        r=zeros(nVox,1);
        for bb=1:nBlock
            fprintf('.');
            left=(bb-1)*blocksize+1;
            right=bb*blocksize;
            if bb==nBlock, right=nVox; end
            tmp_data=data(:,left:right);            
            tmp_data_std=std(tmp_data); tmp_data_std(tmp_data_std==0)=Inf; 
            r(left:right)=(seed - mean(seed))' * (tmp_data-ones(nTp,1)*mean(tmp_data)) ./ std(seed) ./ tmp_data_std ./ (nTp-1);
        end
        if R2Z
            fprintf('\n\tFisher R to Z\n');
            r=0.5*log((1+r)./(1-r)); 
        end
        R(maskDat)=r;
        R=reshape(R,[M,N,O]);
        clear data tmp_data r;
        
        fprintf('\tSave File\n');
        fcHdr=dataHdr;
        fcHdr.fname=[resultDir,filesep,subj,e];
        fcHdr.descrip=sprintf('Functional Connectivity: seed, %s; mask %s;R2Z=%d',R2Z, seedName, maskName);
        spm_write_vol(fcHdr,R);
        fcHdr.dt(1)=32;
        clear R;
        
    elseif D3flag
        fprintf('Unsupported yet\n'); return;
    end
    
    toc;
    fprintf('============================ Done\n');
end

return
end


function data = dzDespike(data, TR)

c1 = 2.5;
c2 = 3;


tc = tc(:);

[lestimates] = icatb_regress(tc,[ones(length(tc),1) (-1:2/(length(tc)-1):1)']);
[qestimates,  modelq] = icatb_myquadfun(tc,TR);
[splestimates,  models] = icatb_mysplinefun(tc,TR);


ylfit =  lestimates(1) + lestimates(2)*(-1:2/(length(tc)-1):1)';
yqfit = icatb_getQuadFit(qestimates,length(tc),TR);
ysfit = icatb_getSplineFit(splestimates,length(tc),TR);

err = [icatb_gfit2(tc,ylfit,'1') icatb_gfit2(tc,yqfit,'1') icatb_gfit2(tc,ysfit,'1')];

[mnerr mnID] = min(err);

if mnID == 1
    yfit =  ylfit;
elseif mnID == 2
    yfit = yqfit;
else
    yfit = ysfit;
end

res = tc - yfit;
mad_res = median(abs(res - median(res))); % median absolute deviation of residuals
sigma = mad_res* sqrt(pi/2);
s = res/sigma;
s_out = s;

ind = find(abs(s) > c1);
for uu = 1:length(ind)
    s_out(ind(uu)) = sign(s(ind(uu)))*(c1+((c2-c1)*tanh((abs(s(ind(uu)))-c1)/(c2-c1))));
end

tc_out = yfit + s_out*sigma;
end