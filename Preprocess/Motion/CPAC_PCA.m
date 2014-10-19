
function dataPCA=CPAC_PCA(Data,Options)

if ~isfield(Options,'PCADim')||isempty(Options.PCADim), error('please specify PCA dimension\n'); end

[AllVolume]=dzReadData_file(Data,Options);

nBlock=10; if isfield(Options,'nBlock')&&~isempty(Options.nBlock)&&Options.nBlock>0, nBlock=Options.nBlock; end
try, Detrend=Options.Flag.Detrend;  catch, Detrend=1; end
try, SigNorm=Options.Flag.SigNorm;  catch, SigNorm=1; end
    
[nTp,nVox]=size(AllVolume); 
fprintf('\t\t%d voxels remained for PCA\n',nVox);
BlockSize=floor(nVox/nBlock);
fprintf('PCA computing\n');
dataCov  = zeros(nTp,nTp); % Covariance of tmpData
fprintf(' step1 covariance\n\t');
for bb=1:nBlock
    % Blockmize
    fprintf('.');
    left=(bb-1)*BlockSize+1; right=bb*BlockSize; if bb==nBlock, right=nVox; end
    tmpData=AllVolume(:,left:right);
%     % rexPCA - step1, compute covariance
%     dataMean = zeros(1,size(tmpData,2));
%     % nan replaced by mean
%     if(any(isnan(tmpData)))
%         nanMat = find(isnan(tmpData));
%         tmpData(nanMat) = 0;
%         dataMean = sum(tmpData,1)./max(eps,sum(~nanMat,1)); % if all nan ,than very large
%         tmpData(nanMat) = dataMean();
%     else
%         datamean = mean(tmpData,1);
%     end
    % detrend
    if Detrend
        tmpData = spm_detrend(tmpData,1);
    end
    % signal normalization
    if SigNorm
        tmpData = (tmpData-repmat(mean(tmpData),size(tmpData,1),1))./repmat(std(tmpData),size(tmpData,1),1);
        tmpData(isnan(tmpData))=0;
    end
    % covariance
    dataCov = dataCov + tmpData * tmpData';
end
fprintf('\n');
% rexPCA - step2, svd & dim
fprintf(' step2 svd\n');
[U,S]=svd(dataCov); % caution!, seemingly strang of svd usage in matlab, if U=svd(dataCov), then U is sigular vector!
if(1<Options.PCADim) % aCompCor 
    if(size(U,2)<Options.PCADim), fprintf('warning: U column %d < Specified PCA dimension %d, select mininum\n',size(U,2),Options.PCADim); end
    dataPCA=U(:,1:min(size(U,2),Options.PCADim));
elseif (Options.PCADim<1)
    dzdummy();
end

fprintf(' done\n');

return
end

