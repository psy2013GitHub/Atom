function [ output_args ] = rexPCA( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ischar(Options.BrainMask), [BrainMaskHdr,BrainMaskDat]=dzReadMask(Options.BrainMask); else, fprintf('Error: Brain Mask Need To Be File Name\n'); return; end % default BrainMask numeric
mask=CovPath(ii); [maskHdr,maskDat]=dzReadMask(mask);
if ~all(maskHdr.mat==BrainMaskHdr.mat),
    fprintf('Warning: %s not match Brain Mask, Resizing\n',mask);
    Outfname=dzImgResize(maskHdr,BrainMaskHdr);
    [maskHdr,maskDat]=dzReadMask(Outfname);
    delete(Outfname); % delete immediate file
end
AllVolume=AllVolume(:,maskDat&Options.BrainMask); % intersect
[nTp,nVox]=size(AllVolume); BlockSize=floor(nVox/nBlock);
for bb=1:nBlock
    % Blockmize
    left=(bb-1)*BlockSize+1; right=bb*nBlock; if bb==nBlock, right=nVox; end
    tmpData=AllVolume(:,left:right);
    
    % rexPCA - step1, compute covariance
    dataMean = zeros(1,size(tmpData,2));
    dataCov  = zeros(nTp,nTp); % Covariance of tmpData
    % nan replaced by mean
    if(any(isnan(tmpData)))
        nanMat = find(isnan(tmpData));
        tmpData(nanMat) = 0;
        dataMean = sum(tmpData,1)./max(eps,sum(~nanMat,1)); % if all nan ,than very large
        tmpData(nanMat) = dataMean();
    else
        datamean = mean(tmpData,1);
    end
    % covariance
    dataCov = dataCov + tmpData * tmpData';
end
% rexPCA - step2, svd & dim
[q1,q2,nill]=svd(dataCov); clear dataCov;
Dim=min(size(q1,2),CovDim(min(ii,length(CovDim))));
dataPCA=q1(:,1:Dim)*diag(sqrt(diag(q2(1:Dim,1:Dim))));
% rexPCA - step3,
basis=rexPCA_step3(dataPCA,AllVolume); % C version, see /src/rexPCA_step3.c
% rexPCA - step4,
basis=basis*diag(1./max(eps,sqrt(sum(basis.^2,1))));
tmp=sign(sum(basis,1))./max(eps,sum(abs(basis),1));
dataPCA=dataPCA*diag(tmp);

end

