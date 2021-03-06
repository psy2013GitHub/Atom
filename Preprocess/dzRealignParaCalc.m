function dzRealignParaCalc(rpTxt,meanImage,OutputDir,ExcludeDir,ExcludeInfo)

%---------------%
% - maxHeadMotion & meanHeadMotion & FD_Jenkinson & Power

%- FD_Jenkinson
[p,f,e]=fileparts(rpTxt);
[relRMS, absRMS]=dzFD_Jenkinson(rpTxt,meanImage);
save([OutputDir,filesep,'FD_Jenkinson.mat'],'relRMS','absRMS');

%- Power 2012 Scrubbing
% convert degrees into motion in mm;
radius=50; % 50mm as default
HeadMotion=load(rpTxt);
HeadMotion(:,4:6)=HeadMotion(:,4:6)*(2*radius*pi/360); % Caution!, spm/AFNI/FSL might be different in estimating rotation parameter!, see YAN, 2013, NeuroImage
FD=diff(HeadMotion);
FD=[zeros(1,size(FD,2)); FD];  % first element is a zero, as per Power et al 2014
FD=sum(abs(FD),2);
save([OutputDir,filesep,'FD_Power.mat'],'FD');


%- maxHeadMotion
maxHeadMotion=max(HeadMotion);
maxHeadMotion(4:6)=maxHeadMotion(4:6)*180/pi;
save([OutputDir,filesep,'maxHeadMotion.mat'],'maxHeadMotion');
%---------------%
%- ExcludeSub_TXT
for ExcludingCriteria=3:-0.5:0.5
    fid = fopen([ExcludeDir,filesep,'ExcludeSubjects_',num2str(ExcludingCriteria),'.txt'],'at+');
    if any(maxHeadMotion>ExcludingCriteria)
       fprintf(fid,'\n%s',ExcludeInfo);
    end
end
fclose(fid);
end