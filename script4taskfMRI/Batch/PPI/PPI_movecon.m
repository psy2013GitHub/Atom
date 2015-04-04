

cwd=pwd;


ConDir ='/home/dengzhou/WM_fMRI/PPI_1stlevel';
DestDir='/home/dengzhou/WM_fMRI/PPI_2ndlevel';
% ConDir ='\\QIU&ZHANG\dengzhou\WM_fMRI\PPI_1stlevel';
% DestDir='\\QIU&ZHANG\dengzhou\WM_fMRI\PPI_2ndlevel';


filter='^con.*(001|002).*\.(img|hdr|nii)$';
flag='dirname';

keyname=[{'fail'},{'unr'},{'suc'}];

%--------------------------------------------------------------------------
cd ..; 
cd Utility/Move_con;
% cd Utility\Move_con; 
for kk=1:length(keyname), Move_con(ConDir,DestDir,filter,flag,keyname{kk}); end; delete *asv; cd(cwd);