

function varargout=dzDetrend(AllVolume,Options)

if nargin<2, Options=struct(); end

fprintf('detrending....\n');
% 4D to 2D
if length(size(AllVolume))==4
    nTp=size(AllVolume,4);
    AllVolume=reshape(AllVolume,[],nTp);
    AllVolume=AllVolume';
end

if length(size(AllVolume))~=2, fprintf('Error: Data Dimension Not 2D or 4D\n'); end

nTp =size(AllVolume,1);
nVox=size(AllVolume,2);

% PolyTrend
p=1; if isfield(Options,'PolyTrend')&&~isempty(Options.PolyTrend)&&isnumeric(Options.PolyTrend), p=Options.PolyTrend; end

% Detrend
AllVolume=spm_detrend(AllVolume,p);
AllVolume=AllVolume+ones(nTp,1)*mean(AllVolume); % remember to add mean

% Varargout
% if 1, return AllVolume
if nargout==1, varargout{1}=AllVolume; fprintf('Done....\n'); return; end
% if 0 ,write to img
if ~nargout
    if ~isfield(Options,'Hdr')||isempty(Options.Hdr)||isstruct(Options.Hdr), fprintf('Error: Specify the Header\n'); end
    dzWrite4DNIfTI(AllVolume,Options.Hdr);
end
fprintf('Done....\n');

end