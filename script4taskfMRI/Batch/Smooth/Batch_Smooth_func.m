function Batch_Smooth_func(root, Subj_lst, nii_filter, Batch_Smooth_Dir, MultiRun_flag, smooth)

for ss=1:length(Subj_lst)
    tmp_subj=Subj_lst(ss);
    if ~tmp_subj.isdir, continue; end

    if MultiRun_flag
        tmp_Dir = [root, filesep, Subj_lst(ss).name]; 
        run_lst = dir(tmp_Dir); run_lst = run_lst(3:end); run_lst = run_lst([run_lst(:).isdir] == 1);
        root1 = [root, filesep, Subj_lst(ss).name];
        Batch_Smooth_func(root1, run_lst,nii_filter,Batch_Smooth_Dir,0, smooth) % recursively
    else
        s = strrep(root,filesep,'_');
        tmp_file=[Batch_Smooth_Dir, filesep, s,'_',  tmp_subj.name, '_', 'tmp_smooth_job.m'];
        fprintf('Writing Batch for %s...\n',tmp_subj.name);
        if exist(tmp_file,'file')==2, delete(tmp_file); end; tmp_fid=fopen(tmp_file,'w');
        fprintf(tmp_fid,'\n\n%%----------------------------------------------------\n');
        fprintf(tmp_fid,'%% Created at %s %s\n',now, date);
        fprintf(tmp_fid,'%%----------------------------------------------------\n\n\n');

        tmp_dir=[root, filesep, tmp_subj.name];
        tmp_nii_files=spm_select('List',tmp_dir,nii_filter);
        nii_files=[repmat([tmp_dir,filesep],[size(tmp_nii_files,1),1]),tmp_nii_files];
        if isempty(tmp_nii_files), continue; end
        nii_files={nii_files};
        
        dzWrite_data(tmp_fid, nii_files);
        dzWrite_options(tmp_fid, smooth);
        fclose(tmp_fid);
        
        template_run_jobs(tmp_file);
        delete(tmp_file);
        
    end
end
return
end


function dzWrite_data(fid, nii_files)
fprintf(fid,'matlabbatch{1}.spm.spatial.smooth.data = {\n');
for rr=1:length(nii_files)
    tmp_nii_files=nii_files{rr};
    for ff=1:size(tmp_nii_files,1)
        fprintf(fid,'                                          ''%s''\n',tmp_nii_files(ff,:));
    end
end
fprintf(fid,'                                          };\n');
return
end

function dzWrite_options(fid, smooth)

if ~isfield(smooth,'fwhm')||isempty(smooth.fwhm), smooth.fwhm=[8,8,8]; end
fprintf(fid,'matlabbatch{1}.spm.spatial.smooth.fwhm = [%d %d %d];\n',smooth.fwhm(1),smooth.fwhm(2),smooth.fwhm(3));
if ~isfield(smooth,'dtype')||isempty(smooth.dtype), smooth.dtype=0; end
fprintf(fid,'matlabbatch{1}.spm.spatial.smooth.dtype = %d;\n',smooth.dtype);
if ~isfield(smooth,'im')||isempty(smooth.im), smooth.im=0; end
fprintf(fid,'matlabbatch{1}.spm.spatial.smooth.im = %d;\n',smooth.im);
if ~isfield(smooth,'prefix')||isempty(smooth.prefix), smooth.prefix='s'; end
fprintf(fid,'matlabbatch{1}.spm.spatial.smooth.prefix = ''%s'';\n',smooth.prefix);

return
end
