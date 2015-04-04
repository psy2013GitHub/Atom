
PDir=['/home/dengzhou/Depression_35_35',filesep,'3DRaw'];
prefix='control_';
suffix='smwc1*';
method=2;
maxNum=3;
txtfile=[pwd,filesep,prefix,'smwc1.txt'];
fid_mode='w';
txt_open=1;

%--------------------------------------------------------------------------
fileList=subList2Path(subList,PDir,prefix,suffix,method,maxNum,txtfile,fid_mode,txt_open);