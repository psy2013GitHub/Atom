


function DVARS=dzPowerDVARS(AllVolume,Options)

if nargin<2, Options=struct(); end
% AllVolume. 4D
nTp=size(AllVolume,4);
AllVolume=reshape(AllVolume,[],nTp);
AllVolume=AllVolume';
if isfield(Options,'mask')&&~isempty(Options.mask),
    if ischar(Options.mask)
       maskHdr=spm_vol(Options.mask); mask=spm_read_vols(maskHdr);
       mask=mask(:);
    elseif isnumeric(Options.mask)
       mask=Options.mask; 
    end
    mask(isnan(mask))=0; mask(isinf(mask))=0; mask=logical(mask);
    AllVolume=AllVolume(:,mask);
end
nVox=size(AllVolume,2);
% convert to signal change
AllVolume=bramila_bold2perc(AllVolume);
% DVARS
Di=[zeros(1,nVox);diff(AllVolume)];
DVARS=sqrt(mean(Di.^2,2));

end

function [y m]=bramila_bold2perc(ts)
% BRAMILA_BOLD2PERC - Converts a time series with mean into a time series of percentage changes.
%   - Usage:
%   ts_perc = bramila_bold2perc(ts) ts is a matrix NxM where N is the
%   number of time points. Values returned are in percentages
%   - Notes:
%   If the mean is zero, then the absolute maximum is used.
%
%   The formula used follows the SPM convention, i.e. we first normalize the
%   time series so that they have 100 as mean value.

%	EG 2014-10-01
m=mean(ts,1);
T=size(ts,1);

% if we have a signal with zero mean, we need to treat it in a special
% way. Zero in our case is 1e5*eps i.e. roughly 10^-11
if(any(m<1e5*eps))
    ids=find(m<1e5*eps);
    m(ids)=max(abs(ts(:,ids)));
    ts(:,ids)=ts(:,ids)+repmat(m(ids),T,1);
end
y=100*(ts./repmat(m,T,1))-100;
y(isnan(y))=0;
end