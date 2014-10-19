function [maskHdr,maskDat]=dzReadMask(mask)

maskHdr=spm_vol(mask); maskDat=spm_read_vols(maskHdr);
maskDat(isnan(maskDat)|isinf(maskDat))=0; maskDat=logical(maskDat(:));

end