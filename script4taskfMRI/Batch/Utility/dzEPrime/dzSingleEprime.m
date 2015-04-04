

clear; clc;
%-------------------------------------------------------------------------%

TxtFname='E:\lhj'; % E-prime??????
minRT=[];            % ?????
maxRT=[];            % ?????
nstd=[];             % ?????
attribute='consistency';


%-------------------------------------------------------------------------%
[Condition, meanRT, ACC] = dzSingleTxtRead(TxtFname,attribute,minRT,maxRT,nstd);