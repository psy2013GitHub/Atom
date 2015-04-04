

clear; clc;
%-------------------------------------------------------------------------%

SubjectDir='J:\temp\EmoConf_txt'; % E-prime??????
OutputDir=pwd;       % ??????
Attribute='sti';
Prefix='stimulus';
minRT=[];            % ?????
maxRT=[];            % ?????
nStd=[-2,2];             % ?????
ACCExlude=1;
Options.Operations.RT =[{}];
Options.Operations.ACC=[{}];


%-------------------------------------------------------------------------%
dzLoopDirectory(SubjectDir,OutputDir,Attribute,Prefix,minRT,maxRT,nStd,ACCExlude,Options);