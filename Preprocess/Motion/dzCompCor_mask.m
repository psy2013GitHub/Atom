


function outMasks=dzCompCor_mask(Masks,Options)

% Masks

if nargin<2, Options=struct(); end
try, Erode  =Options.Erode;  catch, Erode =1;  end
try, Thresh =Options.Thresh; catch, Thresh=1;  end
try, withinALVIN=Options.withinALVIN; catch, withinALVIN=0; end
if iscell(Masks), Masks=strvcat(Masks{:}); end

fprintf('aComCor making mask\n');

mlen=size(Masks,1);
for m=1:mlen
    tmpmask=Masks(m,:);
    [p,f,ext]=fileparts(tmpmask);
    ext=deblank(ext); % strange enough, but still a bug, so 'deblank'
    if ~isempty(strmatch(ext,{'.img','.nii','.hdr'}))
        Hdr=spm_vol(tmpmask); % mask
        Dat=spm_read_vols(Hdr);
        idx1=find(Dat(:)>Thresh);
        [idxx,idxy,idxz]=ind2sub(size(Dat),idx1);
        idxt=find(idxx>Erode&idxx<size(Dat,1)+1-Erode&idxy>Erode&idxy<size(Dat,2)+1-Erode&idxz>Erode&idxz<size(Dat,3)+1-Erode);
        for n1=1:length(idxt),
            if (sum(sum(sum(Dat(idxx(idxt(n1))+(-Erode:Erode),idxy(idxt(n1))+(-Erode:Erode),idxz(idxt(n1))+(-Erode:Erode))<Thresh,3),2),1))>1,
                idxt(n1)=0;
            end
        end
        idxt=idxt(idxt>0);
        idx1=idx1(idxt);
        X1=zeros(size(Dat));
        if withinALVIN
             progp=fileparts(which('Atom.m'));
            [maskHdr,maskDat]=dzReadMask([progp,filesep,'Templates',filesep,'CSF',filesep,'ALVIN_',num2str(abs(det(Hdr.mat))^(1/3)),'mm_v1.nii']);
            idx1=intersect(idx1,find(maskDat));
        end
        X1(idx1)=1;
        Hdr.fname=[p,filesep,'e',f,ext];
        outMasks{m}=Hdr.fname;
        spm_write_vol(Hdr,X1);
    else
        dzdummy(); % fix me
    end
end
if ischar(Masks), outMasks=strvcat(outMasks); end