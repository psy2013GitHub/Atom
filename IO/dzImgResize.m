function Outfname=dzImgResize(fromHdr,toHdr,prefix)

if nargin<3, prefix='tmp_'; end

tmp_arg=struct('mean',false,'interp' ,0,'which',1, 'prefix',prefix);
spm_reslice([toHdr fromHdr],tmp_arg);
[fromP,fromF,fromE]=fileparts(fromHdr.fname);
Outfname=[fromP,filesep,tmp_arg.prefix,fromF,fromE];

end