

function dzSPM_orthoview()
% (underlay,overlay,options)
% screenResolution=get(0,'ScreenSize');

underlay='ch2.nii';
overlay='GMD_10_STI_0.001.img';
options.actions='render';
options.position=[];
options.brt=1;

%-------------------------------------------------------------------------%
if strcmpi(options.actions,'addcolouredblobs')||strcmpi(options.actions,'addtruecolourimage')
    % Prepare
    spm_orthviews('Reset');
    spm_figure('Clear','Graphics');
    
    % Underlay
    P=spm_vol(underlay);
    fg = spm_figure('GetWin','Graphics');
    spm_image('Reset');
    spm_orthviews('Image', P, [0.0 0.45 1 0.55]);
    WS = spm('WinScale');
end

% Overlay
% read data
overHdr=spm_vol(overlay); overDat=spm_read_vols(overHdr);
% Usage: spm_orthviews('AddColouredBlobs',handle,XYZ,Z,mat,colour,name)
% overXYZvox=[overXYZmni(1,:)',overXYZmni(2,:)',overXYZmni(3,:)',ones(size(
% overXYZmni,2),1)]*(inv(overHdr.mat))';
switch(lower(options.actions))
    case 'addcolouredblobs'
        % mesh
        RCP=meshRCP(overHdr);
        overXYZmni=overHdr.mat(1:3,:)*RCP;
        % Overlay
        % cnames  = 'Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs';
        % colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
        spm_orthviews('AddColouredBlobs',1,RCP,overDat(:)',overHdr.mat,[1,0,0],'overlay'); % spm_orthviews('AddColouredImg',...)
        % spm_orthviews('AddColouredBlobs',handle,XYZ,Z,mat,colour,name)
    case 'addtruecolourimage'
        mycmap = colormap([gray(64); jet(64)]);
        spm_orthviews('addtruecolourimage',1,overlay,mycmap,0.4);
    case 'render'
        dat.XYZ=meshRCP(overHdr); dat.t=overDat; dat.mat=overHdr.mat; dat.dim=overHdr.dim;
        h = spm_render(dat,options.brt,fullfile(spm('dir'), 'rend', 'render_single_subj.mat'));
end
% hold on; colorbar('location','south'); hold off;

% Reposition
if isempty(options.position), position=overXYZmni(:,abs(overDat(:))==max(abs(overDat(:)))); else, position=options.position; end
if any(size(position)==1), position=position(:); else, position=position(:,1); end
spm_orthviews('Reposition', position); % jump to a certain coordinate


% Save jpg
if isfield(options,'save_flag')&&options.save_flag
    fs = get(0,'Children');
    res = getframe(fs(1));
    imwrite(res.cdata, 'niftifile.jpg'); % save the window as JPG
end

end


function [RCP]=meshRCP(hdr)
[R,C,P]  = ndgrid(1:hdr.dim(1),1:hdr.dim(2),1:hdr.dim(3));
RCP      = [R(:)';C(:)';P(:)'];
clear R C P
RCP(4,:) = 1;
end
