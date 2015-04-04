




function exist_ps_flag = SingleSubj_ResutDisp_func(SubjDir, DestDir, Batch_ResultDisp_Dir,MultiRun_flag, conspec)

if exist(SubjDir,'dir')~=7, return; end
[DataDir,Subj,u2]=fileparts(SubjDir);
tmp_file=[Batch_ResultDisp_Dir,filesep,'tmp_resultdisp_job.m'];
fprintf('Writing ResultDisp Batch for %s...\n',Subj);
if exist(tmp_file,'file')==2, delete(tmp_file); end; tmp_fid=fopen(tmp_file,'w');
fprintf(tmp_fid,'\n\n%%----------------------------------------------------\n');
fprintf(tmp_fid,'%% Created at %s %s\n',now, date);
fprintf(tmp_fid,'%%----------------------------------------------------\n\n\n');
if MultiRun_flag
    run_list=dir([SubjDir,filesep,Subj]); run_list=run_list(3:end);
    for rr=1:length(run_list)
        tmp_run=run_list(rr);
        if ~tmp_run.isdir||strcmpi(tmp_run.name,'run'), continue; end
        tmp_dir=[SubjDir,filesep,Subj,filesep,tmp_run.name];
        SPMmat_path=[tmp_dir,filesep,'SPM.mat'];
    end
else
    SPMmat_path=[SubjDir,filesep,'SPM.mat'];
end
dzWrite_options(tmp_fid, SPMmat_path, conspec); try, run_job(tmp_file); end; 
exist_ps_flag=dzWrite_mvps([SubjDir,filesep,'spm_',datestr(now,'yyyymmmdd'),'.ps'],DestDir,Subj);
fclose(tmp_fid);
if ~exist_ps_flag, return; end
return
end

function dzWrite_options(fid, SPMmat_path, conspec)

fprintf(fid,'matlabbatch{1}.spm.stats.results.spmmat = {''%s''};\n',SPMmat_path);
fprintf(fid,'matlabbatch{1}.spm.stats.results.conspec.titlestr = ''%s'';\n',conspec.title);
fprintf(fid,'matlabbatch{1}.spm.stats.results.conspec.contrasts = %d;\n',conspec.contrast);
fprintf(fid,'matlabbatch{1}.spm.stats.results.conspec.threshdesc = ''%s'';\n',conspec.thresdesc);
fprintf(fid,'matlabbatch{1}.spm.stats.results.conspec.thresh = %d;\n',conspec.thresh);
fprintf(fid,'matlabbatch{1}.spm.stats.results.conspec.extent = %d;\n',conspec.extend);
fprintf(fid,'matlabbatch{1}.spm.stats.results.conspec.mask = struct(''contrasts'', {}, ''thresh'', {}, ''mtype'', {});\n');
fprintf(fid,'matlabbatch{1}.spm.stats.results.units = 1;\n');
fprintf(fid,'matlabbatch{1}.spm.stats.results.print = %d;\n',conspec.print);

return
end

function exist_ps_flag=dzWrite_mvps(psfile,DestDir,Subj)

exist_ps_flag=1;
if exist(psfile,'file')~=2, fprintf('Error: %s not found\n',psfile); exist_ps_flag=0; return; end
if exist(DestDir,'dir')~=7, mkdir(DestDir); end
[p,f,e]=fileparts(psfile);
movefile(psfile,DestDir); movefile([DestDir,filesep,f,e],[DestDir,filesep,[Subj,'_',f,e]]);

return
end


function run_job(job_file)

nrun = 1; % enter the number of runs here
jobfile = {sprintf('%s',job_file)};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('serial', jobs, '', inputs{:});

return
end


