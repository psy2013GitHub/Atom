

function mask=dzDefaultBrainMask(Hdr)

voxelsize=abs(det(Hdr.mat(1:3,1:3)))^(1/3);
ProgramPath=fileparts(which('Atom.m'));

mask=[ProgramPath,filesep,'Templates',filesep,'Brain',filesep,'BrainMask_',num2str(voxelsize),'mm.nii'];

end