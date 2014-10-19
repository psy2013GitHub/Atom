

function AllVolume=dzCovRegress(AllVolume,Covariables,Options)

%- Covariables
%  .Names/.PCA/.Derivative/.Path

%- Options.mask


% Init
tic;
try, nBlock=Options.nBlock;       catch, nBlock=Options.nBlock; end
try, PolyTrend=Options.PolyTrend; catch, PolyTrend=1;           end

[nTp,nVox]=size(AllVolume);

CovNames=lower(Covariables.Names);
CovDim  =Covariables.Dim;   % 0 for mean, >1 for PCA dimensions
CovThresh=Covariables.Thresh;
CovErode=Covariables.Erode;
CovDer=lower(Covariables.Derivative);
CovPath=Covariables.Path;
CovData=[];

%- dzRegress
for ii=1:length(CovNames)
    fprintf('.');
    tmpcovname=lower(CovNames{ii});
    % motion
    if ~isempty(strfind(tmpcovname,'motion'))
        motion=load(CovPath{ii});
        if CovDim(ii)==6
            break;
        elseif CovDim(ii)==24
            motion=[motion, [zeros(1,size(motion,2));motion(1:end-1,:)], motion.^2, [zeros(1,size(motion,2));motion(1:end-1,:)].^2];
        else
            fprintf('Unknow Motion Parameter');
        end
        if CovDer(ii), motion=[motion,[zeros(1,size(motion,2)),diff(motion(:,1))]]; end
        CovData=[CovData,motion];
    elseif ~isempty(strfind(tmpcovname,'wm'))||~isempty(strfind(tmpcovname,'gm'))||~isempty(strfind(tmpcovname,'csf'))
        % tmpMask
        if ischar(CovPath), tmpMask=CovPath(ii,:); elseif(iscell(CovPath)), tmpMask=CovPath{ii}; end
        if ~CovDim(ii)
            tmpOptions.Mask=tmpMask; tmpOptions.dataHdr=Options.dataHdr;
            dataMean=dzReadData_file(AllVolume,tmpOptions); dataMean=mean(dataMean,2);
            Covariables=[Covariables,dataMean];
        else % aCompCor
            tmpmaskHdr=spm_vol(tmpMask);
            if any(tmpmaskHdr.mat~=Options.dataHdr.mat)
                fprintf(' !!!warning: %s mask dimension not match, resize\n',upper(tmpcovname));
                tmpOptions.Mask=dzImgResize(tmpmaskHdr,Options.dataHdr,['resliced_',num2str(abs(det(Options.dataHdr.mat(1:3,1:3)))^(1/3)),'mm_']);
            end
            % Mask
            tmpOptions.Thresh=CovThresh(ii); tmpOptions.Erode=CovErode(ii);
            if ~isempty(strfind(tmpcovname,'wm')), tmpOptions.withinALVIN=1; end % by default, apply within ALVIN_csf mask
            tmpOptions.Mask=dzCompCor_mask(tmpOptions.Mask,tmpOptions);
            % PCA
            tmpOptions.dataHdr=Options.dataHdr;
            tmpOptions.Flag.Detrend=1; % detrend level
            tmpOptions.Flag.SigNorm=1; % bool
            tmpOptions.PCADim=CovDim(ii);
            try, tmpOptions.VoxIdxInRawMatrix=Options.VoxIdxInRawMatrix; end % in case if 'AllVolume' has already been masked!!!
            dataPCA=CPAC_PCA(AllVolume,tmpOptions);
            CovData=[CovData,dataPCA];
        end
    end
end
fprintf('\n');

% dzRegress
BlockSize=floor(nVox/nBlock);
for bb=1:nBlock
    left=1+(bb-1)*BlockSize; right=bb*BlockSize; if nVox<right, right=nVox; end
    [coef,AllVolume(:,left:right)]=dzRegress(AllVolume(:,left:right),CovData);
end

toc;
fprintf('============================ Done\n');

return
end