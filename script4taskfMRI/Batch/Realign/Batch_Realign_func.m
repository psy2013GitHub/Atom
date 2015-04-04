function Batch_Realign_func(root, Subj_lst,nii_filter,Batch_Realign_Dir,MultiRun_flag, nDiscard_pre, nDiscard_end, eoptions, roptions)
if nargin<5, eoptions=struct(); end
if nargin<6, roptions=struct(); end

% subj_list=dir(DataDir); subj_list=subj_list(3:end);
for ss=1:length(Subj_lst)
    tmp_subj=Subj_lst(ss);
    
    if ~tmp_subj.isdir, continue; end
    
    if MultiRun_flag
        tmp_Dir = [root, filesep, Subj_lst(ss).name]; 
        run_lst = dir(tmp_Dir); run_lst = run_lst(3:end); run_lst = run_lst([run_lst(:).isdir] == 1);
        root1 = [root, filesep, Subj_lst(ss).name];
        Batch_Realign_func(root1, run_lst, nii_filter, Batch_Realign_Dir, 0, nDiscard_pre, nDiscard_end, eoptions, roptions) % trucate files at the begining of every run fix me
    else
        s = strrep(root,filesep,'_');
        tmp_file=[Batch_Realign_Dir,filesep, s, '_', tmp_subj.name, '_', '_tmp_realign_job.m'];
        fprintf('Writing Batch for %s...',tmp_subj.name);
        if exist(tmp_file,'file')==2, delete(tmp_file); end; tmp_fid=fopen(tmp_file,'w');
        fprintf(tmp_fid,'\n\n%%----------------------------------------------------\n');
        fprintf(tmp_fid,'%% Created at %s %s\n',now, date);
        fprintf(tmp_fid,'%%----------------------------------------------------\n\n\n');

        tmp_dir=[root, filesep, tmp_subj.name];
        tmp_nii_files=spm_select('List',tmp_dir,nii_filter);
        tmp_nii_files=tmp_nii_files(nDiscard_pre+1:end-nDiscard_end,:);
        if isempty(tmp_nii_files), continue; end
        
        nii_files=[repmat([tmp_dir,filesep],[size(tmp_nii_files,1),1]),tmp_nii_files];
        nii_files={nii_files};
        
        dzWrite_data(tmp_fid, nii_files);
        dzWrite_options(tmp_fid, eoptions, roptions);
        fclose(tmp_fid);
        
        template_run_jobs(tmp_file);
        delete(tmp_file);
        
    end 
end
return
end


function dzWrite_data(fid, nii_files)
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.data = {\n');
for rr=1:length(nii_files)
    tmp_nii_files=nii_files{rr};
    fprintf(fid,'                                                    {\n');
    for ff=1:size(tmp_nii_files,1)
        fprintf(fid,'                                                    ''%s''\n',tmp_nii_files(ff,:));
    end
    fprintf(fid,'                                                    }\n');
end
fprintf(fid,'                                                    }'';\n');
return
end

function dzWrite_options(fid, eoptions, roptions)
% eoptions
if ~isfield(eoptions,'quality')||isempty(eoptions.quality), eoptions.quality=0.9; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = %d;\n',eoptions.quality);
if ~isfield(eoptions,'sep')||isempty(eoptions.sep), eoptions.sep=4; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = %d;\n',eoptions.sep);
if ~isfield(eoptions,'fwhm')||isempty(eoptions.fwhm), eoptions.fwhm=5; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = %d;\n',eoptions.fwhm);
if ~isfield(eoptions,'rtm')||isempty(eoptions.rtm), eoptions.rtm=1; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = %d;\n',eoptions.rtm);
if ~isfield(eoptions,'interp')||isempty(eoptions.interp), eoptions.interp=2; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = %d;\n',eoptions.interp);
if ~isfield(eoptions,'wrap')||isempty(eoptions.wrap), eoptions.wrap=[0 0 0]; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [%d %d %d];\n',eoptions.wrap(1),eoptions.wrap(2),eoptions.wrap(3));
if ~isfield(eoptions,'weight')||isempty(eoptions.weight), eoptions.weight=''; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = ''%s'';\n',eoptions.weight);

% roptions
if ~isfield(roptions,'which')||isempty(roptions.which), roptions.which=[2 1]; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [%d %d];\n',roptions.which(1),roptions.which(2));
if ~isfield(roptions,'interp')||isempty(roptions.interp), roptions.interp=4; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = %d;\n',roptions.interp);
if ~isfield(roptions,'wrap')||isempty(roptions.wrap), roptions.wrap=[0 0 0]; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [%d %d %d];\n',roptions.wrap(1),roptions.wrap(2),roptions.wrap(3));
if ~isfield(roptions,'mask')||isempty(roptions.mask), roptions.mask=1; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = %d;\n',roptions.mask);
if ~isfield(roptions,'prefix')||isempty(roptions.prefix), roptions.prefix='r'; end
fprintf(fid,'matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = ''%s'';\n',roptions.prefix);

return
end