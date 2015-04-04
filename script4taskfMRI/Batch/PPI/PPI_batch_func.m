
function PPI_batch_func(PreDir, GLMDir, ResultDir, nii_filter, Batch_PPI_Dir, Batch_1stlevel_Dir, Batch_Contrast_Dir, voi, MultiRun_flag, MergeRuns_flag, ppi_specify)
%+++++++++++++++++
% Caution: multiple run should be fixed in future release
%+++++++++++++++++
cwd=pwd;
subj_dirs=dir(GLMDir); subj_dirs=subj_dirs(3:end);
tmp_log_job=[Batch_PPI_Dir,filesep,'ppi_log.txt']; if exist(tmp_log_job,'file'), delete(tmp_log_job); end;
log_fid=fopen(tmp_log_job,'w'); dz_Write_batch_header(log_fid); subj_without_voi=0;
for ii=1:length(subj_dirs)
    tmp_subj_dir=subj_dirs(ii);
    if ~tmp_subj_dir.isdir
        continue;
    end
    fprintf('>> Batch for %s...\n',tmp_subj_dir.name);
    tmp_voi_job=[Batch_PPI_Dir,filesep,'tmp_voi_job.m'];
    tmp_ppi_job=[Batch_PPI_Dir,filesep,'tmp_ppi_job.m'];
    if exist(tmp_ppi_job,'file'), delete(tmp_ppi_job); end
    if exist(tmp_voi_job,'file'), delete(tmp_voi_job); end
    ppi_fid=fopen(tmp_ppi_job,'w'); voi_fid=fopen(tmp_voi_job,'w');
    dz_Write_batch_header(ppi_fid); dz_Write_batch_header(voi_fid);
    tmp_subj_dir_path=[GLMDir,filesep,tmp_subj_dir.name];
    tmp_dest_subj_dir_path=[ResultDir,filesep,tmp_subj_dir.name];
    if exist(tmp_dest_subj_dir_path,'dir')==7, rmdir(tmp_dest_subj_dir_path,'s'); end
    if exist(tmp_dest_subj_dir_path,'dir')~=7, mkdir(tmp_dest_subj_dir_path); end
    nii_files=spm_select('List',tmp_subj_dir_path,'^SPM\.mat$');
    if size(nii_files,1)~=1
        fprintf('%d SPM.mat found in %s\n',size(nii_files,1),tmp_subj_dir.name);
        continue;
    end
    % voi timeseries extraction
    no_voi_flag=dz_Write_batch_voi(tmp_voi_job, voi_fid, tmp_subj_dir_path, tmp_dest_subj_dir_path, voi, cwd);
    if no_voi_flag
        if subj_without_voi==0, fprintf(log_fid,'Sub without voi:\n'); end
        fprintf(log_fid,'                %s\n',tmp_subj_dir.name);
        subj_without_voi=subj_without_voi+1;
        continue;
    end
    % create ppi interaction term
    dz_Write_batch_ppi(tmp_ppi_job, ppi_fid, tmp_subj_dir_path, tmp_dest_subj_dir_path, voi, ppi_specify);
    % 1st level model specification, estimation and contrast,
    % in Batch_1stlevel_Dir
    [append_idx,ppi_specify]=dz_Writw_batch_ppi_specify(tmp_dest_subj_dir_path,ppi_specify);
    dz_Write_batch_ppi_model([PreDir,filesep,tmp_subj_dir.name],ResultDir,nii_filter,Batch_1stlevel_Dir,Batch_Contrast_Dir,MultiRun_flag,MergeRuns_flag,ppi_specify,append_idx);
    ppi_specify.sess.reg(end-2:end)=[];
    fclose(voi_fid); fclose(ppi_fid);
    fprintf('.. ... .\n');
end
fclose(log_fid);
return
end


function dz_Write_batch_header(fid)

fprintf(fid,'\n\n%%----------------------------------------------------\n');
fprintf(fid,'%% Created at %s\n',datestr(now,'yyyymmmdd HH:MM:SS'));
fprintf(fid,'%%----------------------------------------------------\n\n\n');

return
end

% voi timeserie extraction
function no_voi_flag = dz_Write_batch_voi(voi_job, fid, subj_dir_path, dest_subj_dir_path, voi, cwd)

load(fullfile(subj_dir_path,'SPM.mat'));
for vv=1:length(voi)
    fprintf(fid,'\n\nclear matlabbatch;\n');
    fprintf(fid,'matlabbatch{1}.spm.util.voi.spmmat = {''%s''};\n',fullfile(subj_dir_path,'SPM.mat'));
    fprintf(fid,'matlabbatch{1}.spm.util.voi.adjust = %d;\n',voi{vv}.ajust_con);
    fprintf(fid,'matlabbatch{1}.spm.util.voi.session = %d;\n',voi{vv}.sesseion);
    fprintf(fid,'matlabbatch{1}.spm.util.voi.name = ''%s'';\n',voi{vv}.name);
    for rr=1:length(voi{vv}.roi)
        if isfield(eval(sprintf('voi{%d}.roi{%d}',vv,rr)),'spm')
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.spm.spmmat = {''%s''};\n',rr,fullfile(subj_dir_path,'SPM.mat'));
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.spm.contrast = %d;\n',rr,voi{vv}.roi{rr}.spm.con);
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.spm.threshdesc = ''%s'';\n',rr,voi{vv}.roi{rr}.spm.thresdesc);
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.spm.thresh = %d;\n',rr,voi{vv}.roi{rr}.spm.thresh);
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.spm.extent = %d;\n',rr,voi{vv}.roi{rr}.spm.extend);
            % specify inclusive or exlusive contrast herer
            if isfield(eval(sprintf('voi{%d}.roi{%d}.spm',vv,rr)),'mask')
                % ----
                fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.spm.mask.contrast = [',rr);
                for cc=1:length(voi{vv}.roi{rr}.spm.mask.contrast)
                    fprintf(fid,'%d',voi{vv}.roi{rr}.mask.contrast(cc));
                end
                fprintf(fid,' ];\n');
                % ----
                fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.spm.mask.thresh = %d;\n',rr,voi{vv}.roi{rr}.mask.thresh);
                fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.spm.mask.mtype = %d;\n',rr,voi{vv}.roi{rr}.mask.mtype);
            end
        end
        %- self define 'peak_in_mask'
        if isfield(eval(sprintf('voi{%d}.roi{%d}',vv,rr)),'peak_in_mask')&&isfield(eval(sprintf('voi{%d}.roi{%d}.peak_in_mask',vv,rr)),'mask_file')&&~isempty(voi{vv}.roi{rr}.peak_in_mask.mask_file)
            diy_mask=voi{vv}.roi{rr}.peak_in_mask.mask_file;
            if size(diy_mask,1)>1, disp('two or more masks is not implemented yet'); return; end
            mask_hdr=spm_vol(diy_mask); data_hdr=SPM.xY.VY(1); %fix me
            if ~all(isequal(mask_hdr.mat,data_hdr.mat)), mask_dat=dzDCM_mask_resize(data_hdr,mask_hdr);
            else mask_dat=spm_read_vols(mask_hdr); mask_dat(isnan(mask_dat))=0; mask_dat=logical(mask_dat); end
            if ~isfield(eval(sprintf('voi{%d}.roi{%d}.peak_in_mask',vv,rr)),'spm_con')||isempty(voi{vv}.roi{rr}.peak_in_mask.spm_con), fprintf('Error: by ''peak_in_mask'' option, U need tp specify statistcal map\n'); return; end
            stadat=spm_read_vols(spm_vol([SPM.swd,filesep,SPM.xCon(voi{vv}.roi{rr}.peak_in_mask.spm_con).Vspm.fname])); stadat(isnan(stadat))=0; tmp_stadat=stadat(mask_dat);
            max_idx=find(stadat(:)==max(tmp_stadat(:))); [max_X,max_Y,max_Z]=ind2sub(data_hdr.dim,max_idx);
            if length(max_X)~=1, fprintf('%d peaks found in mask, confused\n',length(max_X)); return; end
            maxXYZ=data_hdr.mat*[max_X;max_Y;max_Z;1]; maxXYZ=maxXYZ(1:3,:);
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.centre = [%0.1f %0.1f %.1f];\n',rr,maxXYZ(1), maxXYZ(2),maxXYZ(3));
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.radius = %d;\n',rr,radius);
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.move.fixed = 1;\n',rr);
        end
        %- self define 'peak_in_sphere'
        if isfield(eval(sprintf('voi{%d}.roi{%d}',vv,rr)),'peak_in_sphere') % &&isfield(eval(sprintf('voi{%d}.roi{%d}.peak_in_sphere',vv,rr)),'mask_file')&&~isempty(voi{vv}.roi{rr}.peak_in_mask.mask_file)
            if ~isfield(eval(sprintf('voi{%d}.roi{%d}.peak_in_sphere',vv,rr)),'spm_con')||isempty(voi{vv}.roi{rr}.peak_in_sphere.spm_con), fprintf('Error: by ''peak_in_sphere'' option,you need tp specify statistcal map\n'); return; end
            data_hdr=SPM.xY.VY(1);
            [R,C,P]  = ndgrid(1:data_hdr.dim(1),1:data_hdr.dim(2),1:data_hdr.dim(3));
            RCP      = [R(:)';C(:)';P(:)'];
            clear R C P
            RCP(4,:) = 1;
            XYZmm    = data_hdr.mat(1:3,:)*RCP;
            if ~isfield(eval(sprintf('voi{%d}.roi{%d}.peak_in_sphere',vv,rr)),'centre')||isempty(voi{vv}.roi{rr}.peak_in_sphere.centre), fprintf('Error: by ''peak_in_sphere'' option,you need tp specify centre\n'); return; end
            centre=voi{vv}.roi{rr}.peak_in_sphere.centre; centre=centre(:);
            if ~isfield(eval(sprintf('voi{%d}.roi{%d}.peak_in_sphere',vv,rr)),'radius')||isempty(voi{vv}.roi{rr}.peak_in_sphere.radius), fprintf('Error: by ''peak_in_sphere'' option,you need tp specify radius\n'); return; end
            radius=voi{vv}.roi{rr}.peak_in_sphere.radius;
            mask_idx=find(XYZmm-centre*ones(1,size(XYZmm,2))<=radius);
            stadat=spm_read_vols(spm_vol([SPM.swd,filesep,SPM.xCon(voi{vv}.roi{rr}.peak_in_mask.spm_con).Vspm.fname])); stadat(isnan(stadat))=0; tmp_stadat=stadat(mask_idx);
            max_idx=find(stadat(:)==max(tmp_stadat(:))); [max_X,max_Y,max_Z]=ind2sub(data_hdr.dim,max_idx);
            if length(max_X)~=1, fprintf('%d peaks found in mask, confused\n',length(max_X)); return; end
            maxXYZ=data_hdr.mat*[max_X;max_Y;max_Z;1]; maxXYZ=maxXYZ(1:3,:);
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.centre = [%0.1f %0.1f %.1f];\n',rr,maxXYZ(1), maxXYZ(2),maxXYZ(3));
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.radius = %d;\n',rr,voi{vv}.roi{rr}.peak_in_mask.sphere.radius);
            %             if ~isfield(voi{vv}.roi{rr}.peak_in_mask.sphere,'move')||~isfield(voi{vv}.roi{rr}.peak_in_mask.sphere.move,'global')||~isfield(voi{vv}.roi{rr}.peak_in_mask.sphere.move.global,'spm')...
            %                     ||isempty(voi{vv}.roi{rr}.peak_in_mask.sphere.move.global.spm)
            %                 voi{vv}.roi{rr}.peak_in_mask.sphere.move.global.spm=1;
            %             end
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.move.fixed = 1;\n',rr);
        end
        %- sphere
        if isfield(eval(sprintf('voi{%d}.roi{%d}',vv,rr)),'sphere')
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.centre = [%0.1f %0.1f %.1f];\n',rr,voi{vv}.roi{rr}.sphere.center(1), voi{vv}.roi{rr}.sphere.center(2),voi{vv}.roi{rr}.sphere.center(3));
            fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.radius = %d;\n',rr,voi{vv}.roi{rr}.sphere.radius);
            if isfield(eval(sprintf('voi{%d}.roi{%d}.sphere.move',vv,rr)),'fixed')
                fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.move.fixed = 1;\n',rr);
            elseif isfield(eval(sprintf('voi{%d}.roi{%d}.sphere.move',vv,rr)),'global')
                fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.move.global.spm = %d;\n',rr,voi{vv}.roi{rr}.sphere.move.global.spm);
                fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.move.global.mask = '''';\n',rr);
            elseif isfield(eval(sprintf('voi{%d}.roi{%d}.sphere.move',vv,rr)),'local')
                fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.move.local.spm = %d;\n',rr,voi{vv}.roi{rr}.sphere.move.local.spm);
                fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.move.local.mask = '''';\n',rr);
            elseif isfield(eval(sprintf('voi{%d}.roi{%d}.sphere.move',vv,rr)),'supra')
                fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.move.supra.spm = %d;\n',rr,voi{vv}.roi{rr}.sphere.move.supra.spm);
                fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.sphere.move.supra.mask = '''';\n',rr);
            end
        end
        % fprintf(fid,'matlabbatch{1}.spm.util.voi.expression = ''%s'';\n',voi.expression);
    end
    % to avoid voxels outside of the brain, mask is added in default
    tmp_mask=spm_select('List',subj_dir_path,'^mask\.(img|nii)'); if size(tmp_mask,1)>1, frintf('%d masks found in %s, Please Check\n',size(tmp_mask,1),subj_dir_path); end
    fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.mask.image = {''%s''};\n',rr+1,[subj_dir_path,filesep,tmp_mask]);
    fprintf(fid,'matlabbatch{1}.spm.util.voi.roi{%d}.mask.threshold = 0.5;\n',rr+1); % binary, any vaule in [0:1] is ok
    fprintf(fid,'matlabbatch{1}.spm.util.voi.expression = ''%s&i%d'';\n',voi{1}.expression,rr+1);
    fprintf('\n\n');
end

% move voi_files to PPI directory
no_voi_flag=0;
dzWrite_batch_run_job(voi_job);
voi_filter=[]; for vv=1:length(voi), voi_filter=[voi_filter,'|',voi{vv}.name]; end; voi_filter=['^VOI_.*',voi_filter(2:end),'?.*\.(mat|img|hdr|nii)$']; 
voi_files=spm_select('List',subj_dir_path,voi_filter);
voi_mat=spm_select('List',subj_dir_path,['^VOI.*',voi{1}.name,'?.*\.mat$']);
if size(voi_mat,1)~=1, fprintf(' \nWarning: %d VOI.*mat found in %s, Please Check\n',size(voi_mat,1),subj_dir_path); no_voi_flag=1; return; end
for ff=1:size(voi_files,1), movefile([subj_dir_path,filesep,deblank(voi_files(ff,:))],[dest_subj_dir_path,filesep,deblank(voi_files(ff,:))]); end
cd(cwd);
return
end

% ppi variable construction
function dz_Write_batch_ppi(ppi_job, fid, subj_dir_path, dest_subj_dir_path, voi, ppi_specify)
load(fullfile(subj_dir_path,'SPM.mat'));
%-------------------------
% Suppot only one voi 
%

fprintf(fid,'\n\n');
fprintf(fid,'matlabbatch{1}.spm.stats.ppi.spmmat = {''%s''};\n',[subj_dir_path,filesep,'SPM.mat']);
voi_file=spm_select('List',dest_subj_dir_path,['^VOI.*',voi{1}.name,'.*\.mat$']);
fprintf(fid,'matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {''%s''};\n',[dest_subj_dir_path,filesep,deblank(voi_file)]);
fprintf(fid,'matlabbatch{1}.spm.stats.ppi.type.ppi.u = [');
for rr=1:size(ppi_specify.u,1)
    for cc=1:size(ppi_specify.u,2)
        if rr==1&&cc==1, fprintf(fid,'%d',ppi_specify.u(1)); continue; end
        if cc==1, fprintf(fid,'                                          %d',ppi_specify.u(rr,cc)); continue; end
        if rr~=size(ppi_specify.u,1)&&cc~=size(ppi_specify.u,2), fprintf(fid,' %d',ppi_specify.u(rr,cc)); continue; end
        if rr~=size(ppi_specify.u,1)&&cc==size(ppi_specify.u,2), fprintf(fid,' %d\n',ppi_specify.u(rr,cc)); continue; end
        if rr==size(ppi_specify.u,1)&&cc~=size(ppi_specify.u,2), fprintf(fid,' %d',ppi_specify.u(rr,cc)); continue; end
        if rr==size(ppi_specify.u,1)&&cc==size(ppi_specify.u,2), fprintf(fid,' %d];\n',ppi_specify.u(rr,cc)); continue; end
    end
end
fprintf(fid,'\n');
fprintf(fid,'matlabbatch{1}.spm.stats.ppi.name = ''%s'';\n',ppi_specify.name);
fprintf(fid,'matlabbatch{1}.spm.stats.ppi.disp = %d;\n',ppi_specify.disp);
fprintf(fid,'\n\n');
dzWrite_batch_run_job(ppi_job);
ppi_files=spm_select('List',subj_dir_path,['^PPI_.*',ppi_specify.name,'?.*\.mat$']);
for ff=1:size(ppi_files,1), movefile([subj_dir_path,filesep,ppi_files(ff,:)],[dest_subj_dir_path,filesep,ppi_files(ff,:)]); end
return
end

% ppi regressors specification
function [append_idx,ppi_specify]=dz_Writw_batch_ppi_specify(dest_subj_dir_path,ppi_specify)
ppi_files=spm_select('List',dest_subj_dir_path,['^PPI_.*',ppi_specify.name,'?.*\.mat$']);
if size(ppi_files,1)>1, fprintf('%d PPI.*mat found in %s, Please Check\n',size(ppi_files,1),tmp_subj_dir.name); end
load([dest_subj_dir_path,filesep,ppi_files]);
if isfield(ppi_specify,'sess')&&isfield(ppi_specify.sess,'reg')&&~isempty([ppi_specify.sess.reg(:)])
    append_idx=length(ppi_specify.sess.reg);
else
    append_idx=0;
end
ppi_specify.sess.reg(append_idx+1).name='PPI';  ppi_specify.sess.reg(append_idx+1).val=PPI.ppi;
ppi_specify.sess.reg(append_idx+2).name='Psy';  ppi_specify.sess.reg(append_idx+2).val=PPI.Y;
ppi_specify.sess.reg(append_idx+3).name='Phy';  ppi_specify.sess.reg(append_idx+3).val=PPI.P;

return
end

% ppi modeal specification and estimation
function dz_Write_batch_ppi_model(PreSubjDir,ResultDir,nii_filter,Batch_1stlevel_Dir,Batch_Contrast_Dir,MultiRun_flag,MergeRuns_flag,ppi_specify,append_idx)
cwd=pwd;
[useless1,Subj,useless2]=fileparts(PreSubjDir);
% Model Specification and Estimation
if ~exist(nii_filter,'var'), nii_filter='^swr.*\.(img|nii)$'; end
cd(Batch_1stlevel_Dir); SingleSubj_1stlevel_func(PreSubjDir,ResultDir,nii_filter,Batch_1stlevel_Dir,MultiRun_flag, MergeRuns_flag, ppi_specify);
% Contrast Manager
contrast.consess{1}.tcon.name='posPPI_interaction';
contrast.consess{1}.tcon.convec=[zeros(1,append_idx), 1,0]; % right padd with 0.
contrast.consess{2}.tcon.name='negPPI_interaction';
contrast.consess{2}.tcon.convec=[zeros(1,append_idx),-1,0]; % right padd with 0.
contrast.delete=1;
cd(Batch_Contrast_Dir); SingleSubj_Contrast_func([ResultDir,filesep,Subj],Batch_Contrast_Dir,MultiRun_flag, contrast);
% run template;

cd(cwd);
return
end


function mask_dat=dzDCM_mask_resize(data_hdr,mask_hdr)
disp('resizing mask');
tmp_arg=struct('mean',false,'interp' ,0,'which',1, 'prefix','tmp_' );
spm_reslice([data_hdr mask_hdr],tmp_arg);
[mask_path,mask_fname,mask_ext]=fileparts(mask_hdr.fname);
tmp_mask_fname=[mask_path,filesep,tmp_arg.prefix,mask_fname,mask_ext];
mask_dat=spm_read_vols(spm_vol(tmp_mask_fname));
mask_dat=logical(mask_dat); mask_dat(isnan(mask_dat))=0;
delete(tmp_mask_fname);
return
end


function dzWrite_batch_run_job(job_name)

nrun=1;
jobfile = {job_name};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('serial', jobs, '', inputs{:});

return
end


