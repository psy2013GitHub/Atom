


function FD=dzPowerFD(SixRegidBodyPara,Options)

% SixRegidBodyPara, tx,ty,tx,rx,ry,rz

if nargin<2, Options=struct(); end

radius=50; % 50mm by default
if isfield(Options,'radius')&&~isempty(Options.radius)&&isnumeric(Options.radius), radius=Options.radius; end

% convert from degree to mm
SixRegidBodyPara(4:6)=SixRegidBodyPara(4:6)*(2*pi*radius/360);

% FD
Dts=[zeros(1,6);diff(SixRigidBodyPara)];
FD=sum(abs(Dts),2);

end