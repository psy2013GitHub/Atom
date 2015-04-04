
function SingleSubj_Contrast_func(SubjDir, Batch_Contrast_Dir, MultiRun_flag, contrast)
%
%contrast: like:
%               contrast     
%


if exist(SubjDir,'dir')~=7, return; end

s = strrep(SubjDir,filesep,'_');

[DataDir,Subj,u2]=fileparts(SubjDir);

tmp_file=[Batch_Contrast_Dir, filesep, '_', s, '_', 'tmp_contrast_job.m'];
fprintf('Writing Contrast Batch for %s...\n',Subj);

if exist(tmp_file,'file')==2, delete(tmp_file); end; 
tmp_fid=fopen(tmp_file,'w');
fprintf(tmp_fid,'\n\n%%----------------------------------------------------\n');
fprintf(tmp_fid,'%% Created at %s %s\n',now, date);
fprintf(tmp_fid,'%%----------------------------------------------------\n\n\n');

if MultiRun_flag
    run_list=dir([ResultDir,filesep,Subj]); run_list=run_list(3:end);
    for rr=1:length(run_list)
        tmp_run=run_list(rr);
        if ~tmp_run.isdir||strcmpi(tmp_run.name,'run'), continue; end
        tmp_dir=[ResultDir,filesep,Subj,filesep,tmp_run.name];
        SPMmat_path=[tmp_dir,filesep,'SPM.mat'];
    end
else
    SPMmat_path=[SubjDir,filesep,'SPM.mat'];
end

fprintf(tmp_fid,'matlabbatch{1}.spm.stats.con.spmmat = {''%s''};\n',SPMmat_path);
dzWrite_options(tmp_fid, contrast, SPMmat_path);
fclose(tmp_fid);

display(tmp_file)
template_run_jobs(tmp_file);
delete(tmp_file);

return
end

function dzWrite_options(fid, contrast, SPMmat_path)
for cc=1:length(contrast.consess)
    tmp_contrast=contrast.consess{cc};
    if isfield(tmp_contrast,'tcon')
        % auto-get 'convec' field
        if isfield(tmp_contrast.tcon,'cond_list')
            if ~isfield(tmp_contrast.tcon,'cond_weight')
                error('you must provide ''cond_list'' along with ''weight'' field')
            end
            display('loading SPMmat')
            SPMmat = load(SPMmat_path);
            nCond_in_SPMmat = length(SPMmat.SPM.Sess.U);
            if isfield(tmp_contrast.tcon,'convec')
                error('you cant provide ''cond_list'' , ''convec'' together')
            end
            contrast.consess{cc}.tcon.convec = zeros(1,nCond_in_SPMmat);
            for c = 1:length(tmp_contrast.tcon.cond_list)
                found_flag = 0;
                for i = 1:nCond_in_SPMmat
                    if strcmpi(tmp_contrast.tcon.cond_list{c}, SPMmat.SPM.Sess.U(i).name{1})
                         found_flag = 1;
                        break
                    end
                end
                if ~found_flag % not found
                    error('no %s found in SPMmat.SPM.Sess.U', tmp_contrast.tcon.cond_list{c})
                end
                contrast.consess{cc}.tcon.convec(i) = tmp_contrast.tcon.cond_weight(c);
            end
        end
        tmp_contrast=contrast.consess{cc};
        
        fprintf(fid,'matlabbatch{1}.spm.stats.con.consess{%d}.tcon.name = ''%s'';\n',cc,tmp_contrast.tcon.name);
        fprintf(fid,'matlabbatch{1}.spm.stats.con.consess{%d}.tcon.convec = [%d',cc,tmp_contrast.tcon.convec(1));
        for vv=2:length(tmp_contrast.tcon.convec)
            fprintf(fid,' %d',tmp_contrast.tcon.convec(vv));
        end
        fprintf(fid,'];\n');
        if ~isfield(tmp_contrast.tcon,'sessrep')||isempty(tmp_contrast.tcon.sessrep), tmp_contrast.tcon.sessrep='none'; end
        fprintf(fid,'matlabbatch{1}.spm.stats.con.consess{%d}.tcon.sessrep = ''%s'';\n',cc,tmp_contrast.tcon.sessrep);
    elseif isfield(tmp_contrast,'fcon')
        fprintf(fid,'matlabbatch{1}.spm.stats.con.consess{%d}.fcon.name = ''%s'';\n',cc,tmp_contrast.fcon.name);
        fprintf(fid,'matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = {\n');
        fprintf(fid,'                                                       [%d',tmp_contrast.fcon.convec(1));
        for ii=1:size(tmp_contrast.fcon.convec,1)
            for jj=1:size(tmp_contrast.fcon.convec,2)
                if ii==1&&jj==1, continue; end
                fprintf(fid,' %d',tmp_contrast.fcon.convec(ii,jj));
            end
            if ii~=size(tmp_contrast.fcon.convec,1), fprintf(fid,'\n                                                      '); end
            if ii==size(tmp_contrast.fcon.convec,1), fprintf(fid,'];\n'); end
        end
        fprintf(fid,'                                                       }'';\n');
        if ~isfield(tmp_contrast.fcon,'sessrep')||isempty(tmp_contrast.fcon.sessrep), tmp_contrast.fcon.sessrep='none'; end
        fprintf(fid,'matlabbatch{1}.spm.stats.con.consess{%d}.fcon.sessrep = ''%s'';\n',cc,tmp_contrast.fcon.sessrep);
    end
end
if ~isfield(contrast,'delete')||isempty(contrast.delete), contrast.delete=1; end
fprintf(fid,'matlabbatch{1}.spm.stats.con.delete = %d;\n',contrast.delete);
return
end
