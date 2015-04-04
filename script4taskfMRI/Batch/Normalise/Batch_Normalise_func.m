function Batch_Normalise_func(root, Subj_lst, nii_filter, MultiRun_flag, Batch_Normalise_Dir, eoptions, roptions)
% nii_filter.img
% nii_filter.src
% 
for ss=1:length(Subj_lst)
    tmp_subj=Subj_lst(ss);
    if ~tmp_subj.isdir, continue; end

    if MultiRun_flag
        tmp_Dir = [root, filesep, Subj_lst(ss).name]; 
        run_lst = dir(tmp_Dir); run_lst = run_lst(3:end); run_lst = run_lst([run_lst(:).isdir] == 1);
        root1 = [root, filesep, Subj_lst(ss).name];
        Batch_Normalise_func(root1, run_lst,nii_filter,0,Batch_Normalise_Dir, eoptions, roptions) % set MultiRun_flag = 0
    else
        s = strrep(root,filesep,'_');
        tmp_file=[Batch_Normalise_Dir, filesep, s, '_',  tmp_subj.name, '_', 'tmp_normalise_job.m'];
        display(tmp_file)
        fprintf('Writing Batch for %s...\n',tmp_subj.name);
        if exist(tmp_file,'file')==2, delete(tmp_file); end; 
        tmp_fid=fopen(tmp_file,'w');
        fprintf(tmp_fid,'\n\n%%----------------------------------------------------\n');
        fprintf(tmp_fid,'%% Created at %s %s\n',now, date);
        fprintf(tmp_fid,'%%----------------------------------------------------\n\n\n');
        
        tmp_dir=[root, filesep, tmp_subj.name];
        tmp_nii_files=spm_select('List',tmp_dir,nii_filter.img);
        tmp_nii_files=[repmat([tmp_dir,filesep],[size(tmp_nii_files,1),1]),tmp_nii_files];
        nii_files={tmp_nii_files};
        clear tmp_nii_files;
        source_img=spm_select('List',tmp_dir,nii_filter.src); if isempty(source_img), continue; end
        if size(source_img,1)>1, fprintf('%d mean imgs found in %s\n',size(source_img,1),tmp_subj.name); end
        fprintf(tmp_fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source = {''%s''};\n',[tmp_dir,filesep,source_img]);
        % fix me
        wtsrc=''; fprintf(tmp_fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.subj.wtsrc = ''%s'';\n',wtsrc);
        if isempty(nii_files), continue; end
        
        dzWrite_data(tmp_fid, nii_files);
        dzWrite_options(tmp_fid, eoptions, roptions);
        fclose(tmp_fid);
        display(tmp_file)
        
        template_run_jobs(tmp_file);
        delete(tmp_file);

    end
end
return
end


function dzWrite_data(fid, nii_files)
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {\n');
for rr=1:length(nii_files)
    tmp_nii_files=nii_files{rr};
    for ff=1:size(tmp_nii_files,1)
        fprintf(fid,'                                                               ''%s''\n',tmp_nii_files(ff,:));
    end
end
fprintf(fid,'                                                               };\n');
return
end

function dzWrite_options(fid, eoptions, roptions)

% eoptions
if ~isfield(eoptions,'template')||isempty(eoptions.template), p=fileparts(which('spm')); eoptions.template=[p,filesep,'templates',filesep,'EPI.nii']; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = {''%s''};\n',eoptions.template);
if ~isfield(eoptions,'weight')||isempty(eoptions.weight), eoptions.weight=''; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.weight = ''%s'';\n',eoptions.weight);
if ~isfield(eoptions,'smosrc')||isempty(eoptions.smosrc), eoptions.smosrc=8; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smosrc = %d;\n',eoptions.smosrc);
if ~isfield(eoptions,'smoref')||isempty(eoptions.smoref), eoptions.smoref=0; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smoref =%d;\n',eoptions.smoref);
if ~isfield(eoptions,'regtype')||isempty(eoptions.regtype), eoptions.regtype='mni'; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.regtype = ''%s'';\n',eoptions.regtype);
if ~isfield(eoptions,'cutoff')||isempty(eoptions.cutoff), eoptions.cutoff=25; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.cutoff = %d;\n',eoptions.cutoff);
if ~isfield(eoptions,'nits')||isempty(eoptions.nits), eoptions.nits=16; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.nits = %d;\n',eoptions.nits);
if ~isfield(eoptions,'reg')||isempty(eoptions.reg), eoptions.reg=1; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = %d;\n',eoptions.reg);

% roptions
if ~isfield(roptions,'preserve')||isempty(roptions.preserve), roptions.preserve=0; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.preserve = %d;\n',roptions.preserve);
if ~isfield(roptions,'bb')||isempty(roptions.bb), roptions.bb=[-90,-126,-72;90,90,108]; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = [%d %d %d\n',roptions.bb(1,1),roptions.bb(1,2),roptions.bb(1,3));
fprintf(fid,'                                                             %d %d %d];\n',roptions.bb(2,1), roptions.bb(2,2), roptions.bb(2,3));
if ~isfield(roptions,'vox')||isempty(roptions.vox), roptions.vox=[3,3,3]; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [%d %d %d];\n',roptions.vox(1),roptions.vox(2),roptions.vox(3));
if ~isfield(roptions,'interp')||isempty(roptions.interp), roptions.interp=1; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = %d;\n',roptions.interp);
if ~isfield(roptions,'wrap')||isempty(roptions.wrap), roptions.wrap=[0,0,0]; end
fprintf(fid,'matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [%d %d %d];\n',roptions.wrap(1),roptions.wrap(2),roptions.wrap(3));
if ~isfield(roptions,'prefix')||isempty(roptions.prefix), roptions.prefix='w'; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = ''%s'';\n',roptions.prefix);

return
end