function [AllVolume] =dzReadData_file(Data,Options)

% read Images || mask Data

% identify type
if ischar(Data)||iscell(Data)
    type='Imgs';
elseif isnumeric(Data)
    type='Dat';
end

%- Load Data
if ~strcmpi(type,'Dat')
    nImgs=size(Data,1); fprintf('Read Data\n');
    for ff=1:nImgs
        if(~mod(ff,5)), fprintf('.'); end
        D3flag=0; D4flag=0;
        try, file=Data(ff,:); catch, file=Data{ff}; end
        hdr=spm_vol(file);
        % data dimension, 3D or 4D?
        if length(hdr)>1 %4D
            D4flag=1;
            nTp=size(data,4);
        elseif length(hdr)==1 %3D
            D3flag=1;
        end
        if ff==1
            % init AllVolume
            dataHdr=hdr(1); % hdr(1) in case of 4D
            M=dataHdr.dim(1); N=dataHdr.dim(2); O=dataHdr.dim(3);
            if D3flag
                AllVolume=zeros(nImgs,M*N*O);
            elseif D4flag
                AllVolume=cell(nImgs,1);
            end
        else % hdr(1) in case of 4D
            if ~all(dataHdr.mat==hdr(1).mat), error('Error: image dimension not match\n'); end
        end
        % read data
        data=spm_read_vols(hdr);
        if D3flag
            AllVolume(ff,:)=data(:)';
        elseif D4flag
            AllVolume{ff}=data; % each 4D as a cell
        end
    end
else % Dat Formalt
    if(~isfield(Options,'dataHdr')&&isempty(Options.dataHdr)), error('Please supply data header'); end
    if size(Data,2)~=prod(Options.dataHdr.dim)
       try
           VoxIdxInRawMatrix=Options.VoxIdxInRawMatrix;
       catch
           error('Please supply voxel idx of raw data in %d*%d*%d',Options.dataHdr.dim(1),Options.dataHdr.dim(2),Options.dataHdr.dim(3));
       end
    end
    dataHdr=Options.dataHdr;
end

%- maskDat
if isfield(Options,'Mask')&&~isempty(Options.Mask)
    fprintf('\nApply mask');
    for mm=1:size(Options.Mask,1)
        try, mask=Options.Mask(mm,:); catch, mask=Options.Mask{mm}; end
        resizeFlag=0; maskHdr=spm_vol(mask);
        if any(maskHdr.mat~=dataHdr.mat)
            fprintf(' !!!warning: mask %d dimension not match, resize\n',mm);
            resizeFlag=1; mask=dzImgResize(maskHdr,dataHdr,['resliced_',num2str(det(maskHdr.mat(1:3,1:3))),'mm_']);
        end
        [maskHdr,tmpmaskDat]=dzReadMask(mask);
%       if resizeFlag, delete(mask); end
        if mm==1, maskDat=tmpmaskDat; else, maskDat=maskDat&tmpmaskDat; end
    end
    nVox=sum(maskDat);
end
if strcmpi(type,'Dat')&&exist('VoxIdxInRawMatrix','var')
    if all(size(VoxIdxInRawMatrix)~=1)&&islogical(VoxIdxInRawMatrix)
        VoxIdxInRawMatrix=find(VoxIdxInRawMatrix);
    else
        VoxIdxInRawMatrix=sort(VoxIdxInRawMatrix);
    end
    maskVoxIdx=find(maskDat);
    [useless1,useless2,maskVoxIdx]=intersect(maskVoxIdx,VoxIdxInRawMatrix);
    clear useless*;
else
    maskVoxIdx=find(maskDat);
end
clear maskDat;

%- Mask AllVolume
if strcmpi(type,'Dat')
    AllVolume=Data(:,maskVoxIdx);
    fprintf('\nDone\n');
    return;
end
if iscell(AllVolume)
    for cc=1:length(AllVolume)
        AllVolume{cc}=AllVolume{cc}(:,maskVoxIdx);
    end
else
    AllVolume=AllVolume(:,maskVoxIdx);
end
fprintf('\nDone\n');
return;


end

