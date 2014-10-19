function [AllVolume,dataHdrAll,maskHdr,maskDat]=dzReadAllVolume(funImgs,mask)


% delete ... delete ... delete ... delete ...
if nargin<2, mask=''; end
if iscell(funImgs), nscans=length(funImgs); elseif ischar(funImgs), nscans=size(funImgs,1); end 
fprintf('\nRead Data\n\t');
dataHdrAll=struct();
for ii=1:nscans
    if nscans<20, fprintf('.'); elseif mod(ii,5), fprintf('.'); end
    try, Img=funImgs{ii}; catch, Img=funImgs(ii,:); end
    hdr=spm_vol(Img); dat=spm_read_vols(hdr); dat=dat(:)';
    if ii==1, 
        dataHdr=hdr; dataHdrAll=hdr;
        if isempty(mask), 
            BrainMask=dzDefaultBrainMask(dataHdr); % will create a default mask
            [maskHdr,maskDat]=dzReadMask(BrainMask);
        end
        AllVolume=zeros(nscans,sum(maskDat));
    else
        dataHdrAll=[dataHdrAll;hdr];
    end
    AllVolume(ii,:)=dat(maskDat);
end
AllVolume(isnan(AllVolume)|isinf(AllVolume))=0; % zero nan/inf
fprintf('\n');

end