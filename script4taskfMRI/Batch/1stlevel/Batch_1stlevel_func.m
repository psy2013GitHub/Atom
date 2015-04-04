


function Batch_1stlevel_func(root, Subj_lst, ResultDir, nii_filter, Batch_1stlevel_Dir, MultiRun_flag, MergeRun_flag, specify, subj_specific_covariates)

RawCondLen = 0;
if isfield(specify, 'sess') && isfield(specify.sess, 'cond')
    RawCondLen = length(specify.sess.cond);
end
display(['RawCondLen ', num2str(RawCondLen)]);

RawRegLen = 0;
if isfield(specify, 'sess') && isfield(specify.sess, 'regress')
    RawRegLen = length(specify.sess.regress);
end
display(['RawRegLen ', num2str(RawRegLen)]);

if MultiRun_flag
    for ss = 1:length(Subj_lst)
        subjDir = [root, filesep, Subj_lst(ss).name];
        display(subjDir)
        run_lst = dir(subjDir); run_lst = run_lst(3:end); run_lst = run_lst([run_lst(:).isdir]==1);
        
        %##--## Merge Runs ##--##
        if ~MergeRun_flag
            % fix ResultDir
            Batch_1stlevel_func(root, run_lst, ResultDir, nii_filter, Batch_1stlevel_Dir, 0, MergeRun_flag, specify, subj_specific_covariates)
            return
        end
        
        %##--## No Merge Runs ##--##
        if subj_specific_covariates
            covariates_dir = [subjDir, filesep, 'covariates']; % fix with subject-specific covariates
            subj_cond = dzReadCovariates(covariates_dir);
            specify.sess.cond(RawCondLen+1:RawCondLen+length(subj_cond)) = subj_cond; 
        end
        % rp parameter
        if specify.regress_motion
            motion_regressor = merge_rp_para(subjDir, run_lst);
            specify.sess.regress(RawRegLen+1 : RawRegLen+length(motion_regressor)) = motion_regressor;
        end
        % data files
        nii_files = {};
        for run = 1:length(run_lst)
            run_dir = [subjDir, filesep, run_lst(run).name];
            tmp_nii_files = spm_select('List', run_dir, nii_filter);
            if isempty(tmp_nii_files)
                display(['no files in ', run_dir]) 
                continue
            end
            tmp_nii_files = [repmat([run_dir, filesep],[size(tmp_nii_files,1),1]), tmp_nii_files];
            nii_files = [nii_files, {tmp_nii_files}];
        end
        if isempty(nii_files)
           display(['no files in ', subjDir])
           return % skip this subject
        end
        % batch file name
        s = strrep(subjDir, filesep, '_');
        tmp_file = [Batch_1stlevel_Dir, filesep, s, '_tmp_1stlevel_job.m'];
        single_1stlevel_batch(tmp_file, [ResultDir, filesep, Subj_lst(ss).name], nii_files, specify);
        
    end
else % single run
    if subj_specific_covariates
        covariates_dir = [subjDir, filesep, 'covariates']; % fix with subject-specific covariates
        subj_cond = dzReadCovariates(covariates_dir);
        specify.sess.cond(RawCondLen+1:RawCondLen+length(subj_cond)) = subj_cond;
    end
    % rp parameter
    if specify.regress_motion
        motion_regressor = merge_rp_para(root, {subjDir});
        specify.sess.regress(RawRegLen+1 : RawRegLen+length(motion_regressor)) = motion_regressor;
    end
    % data
    nii_files = spm_select('List', subjDir, nii_filter);
    if isempty(nii_files)
        return
    end
    s = strrep(subjDir, filesep, '_');
    tmp_file = [Batch_1stlevel_Dir, filesep, s, '_tmp_1stlevel_job.m'];
    single_1stlevel_batch(tmp_file, [ResultDir, filesep, Subj_lst(ss).name], nii_files, specify);
    
end
return
end


function subj_cond = dzReadCovariates(covariates_dir)

subj_cond = struct('name',[],'onset',[],'duration',[]);
cond_lst = {};
cond_name_file = [covariates_dir, filesep, '_cond_name'];
fid = fopen(cond_name_file);
if fid<0
    display([cond_name_file, ' unaccessble']);
    error('hha');
end
while 1
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    cond_lst = [cond_lst, {deblank(line)}];
end

for c = 1:length(cond_lst)
    
   c_name = cond_lst{c}; 
   c_onset = load([covariates_dir, filesep, c_name, '_onset']);
   c_duration = load([covariates_dir, filesep, c_name, '_duration']);
   
   if length(c_onset) ~= length(c_duration)
       error(['wrong with onset & duration file of ',c_name])
   end
   
   if isempty(c_onset)
       continue
   end
   
   cond.name = c_name;
   cond.onset = c_onset;
   cond.duration = c_duration;
   
   subj_cond(end+1) = cond;

end

try
    subj_cond = subj_cond(2:end);
catch
    subj_cond = [];
end

end


function single_1stlevel_batch(tmp_file, ResultDir, nii_files, specify)

if exist(tmp_file,'file')==2, delete(tmp_file); end;
tmp_fid=fopen(tmp_file,'w');


fprintf(tmp_fid,'\n\n%%----------------------------------------------------\n');
fprintf(tmp_fid,'%% Created at %s %s\n',now, date);
fprintf(tmp_fid,'%%----------------------------------------------------\n\n\n');

if exist(ResultDir,'dir')==7, rmdir(ResultDir,'s'); end
mkdir(ResultDir);

fprintf(tmp_fid,'matlabbatch{1}.spm.stats.fmri_spec.dir = {''%s''};\n',ResultDir);

dzWrite_timing(tmp_fid, specify);
dzWrite_data(tmp_fid, nii_files);
dzWrite_specify_design(tmp_fid,specify);
dzWrite_specify_options(tmp_fid,specify);
dzWrite_estimate_options(tmp_fid);
fclose(tmp_fid);

template_run_jobs(tmp_file);
delete(tmp_file);


return
end

function dzWrite_timing(fid, specify)
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.timing.units = ''%s'';\n',specify.timing.units);
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.timing.RT = %d;\n',specify.timing.RT);
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = %d;\n',specify.timing.fmri_t);
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = %d;\n',specify.timing.fmri_t0);
return
end

function dzWrite_data(fid, nii_files)
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {\n');
if iscell(nii_files), rrl=length(nii_files); else, rrl=1; nii_files={nii_files}; end % merge runs or not merge respectivly
for rr=1:rrl
    tmp_nii_files=nii_files{rr};
    for ff=1:size(tmp_nii_files,1)
        fprintf(fid,'                                          ''%s''\n',tmp_nii_files(ff,:));
    end
end
fprintf(fid,'                                                 };\n');
return
end


function dzWrite_specify_design(fid,specify)

% Condition Section

if ~isfield(specify,'sess')||~isfield(specify.sess,'cond')||isempty([specify.sess.cond(:)]),
    fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct(''name'', {}, ''onset'', {}, ''duration'', {}, ''tmod'', {}, ''pmod'', {});\n');
end
if isfield(specify,'sess')&&isfield(specify.sess,'cond')&&~isempty([specify.sess.cond(:)])
    dzWrite_specify_conditions(fid,specify);
end
if ~isfield(specify,'sess')||~isfield(specify.sess,'multi')||isempty([specify.sess.multi(:)]), specify.sess.multi=''; end
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''%s''};\n',specify.sess.multi);

% Regressor Section
if ~isfield(specify,'sess')||~isfield(specify.sess,'regress')||isempty([specify.sess.regress(:)]),
    fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct(''name'', {}, ''val'', {});\n');
end
if isfield(specify,'sess')&&isfield(specify.sess,'regress')&&~isempty([specify.sess.regress(:)])
    dzWrite_specify_regressors(fid,specify);
end
if ~isfield(specify,'sess')||~isfield(specify.sess,'multi_reg')||isempty([specify.sess.multi_reg(:)]), specify.sess.multi_reg=''; end
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''%s''};\n',specify.sess.multi_reg);

% Factor Section fix me
if ~isfield(specify,'fact')||isempty([specify.fact(:)]),
    fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.fact = struct(''name'', {}, ''levels'', {});\n');
end

% hpf
if ~isfield(specify,'sess')||~isfield(specify.sess,'hpf')||isempty(specify.sess.hpf), specify.sess.hpf=128; end
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = %d;\n',specify.sess.hpf);

return
end


function dzWrite_specify_regressors(fid,specify)

% regress_section
for rr=1:length(specify.sess.regress)
    fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.regress(%d).name = ''%s'';\n',rr,specify.sess.regress(rr).name);
    % onset
    if length(specify.sess.regress(rr).val)==1
        fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.regress(%d).val = %d;\n',rr,specify.sess.regress(rr).val(1));
    else
        fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.regress(%d).val = [%d\n',rr,specify.sess.regress(rr).val(1));
        for vv=2:length(specify.sess.regress(rr).val)-1
            fprintf(fid,'                                                         %d\n',specify.sess.regress(rr).val(vv));
        end
        if isempty(vv), vv = 1; end %!!! if length(...)==2
        fprintf(fid,'                                                         %d];\n',specify.sess.regress(rr).val(vv));
    end
end

return
end

function dzWrite_specify_conditions(fid,specify)
% Condition Section
for cc=1:length(specify.sess.cond)
    fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.cond(%d).name = ''%s'';\n',cc,specify.sess.cond(cc).name);
    % onset
    specify.sess.cond(cc).name
    specify.sess.cond(cc).onset
    if length(specify.sess.cond(cc).onset)==1
        fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.cond(%d).onset = %d;\n',cc,specify.sess.cond(cc).onset(1));
    else
        fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.cond(%d).onset = [%d\n',cc,specify.sess.cond(cc).onset(1));
        for oo=2:length(specify.sess.cond(cc).onset)-1
            fprintf(fid,'                                                         %d\n',specify.sess.cond(cc).onset(oo));
        end
        if isempty(oo), oo = 1; end % !!! if length(specify.sess.cond(cc).onset)==2
        fprintf(fid,'                                                         %d];\n',specify.sess.cond(cc).onset(oo+1));
    end
    % duration
    if length(specify.sess.cond(cc).duration)==1
        fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.cond(%d).duration = %d;\n',cc,specify.sess.cond(cc).duration(1));
    else
        fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.cond(%d).duration = [%d\n',cc,specify.sess.cond(cc).duration(1));
        for dd=2:length(specify.sess.cond(cc).duration)-1
            fprintf(fid,'                                                         %d\n',specify.sess.cond(cc).duration(dd));
        end
        if isempty(dd), dd = 1; end
        fprintf(fid,'                                                         %d];\n',specify.sess.cond(cc).duration(dd+1));
    end
    if ~isfield(specify.sess.cond(cc),'tmod')||isempty(specify.sess.cond(cc).tmod), specify.sess.cond(cc).tmod=0; end
    fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.cond(%d).tmod = %d;\n',cc,specify.sess.cond(cc).tmod);
    % fix me
    fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.sess.cond(%d).pmod = struct(''name'', {}, ''param'', {}, ''poly'', {});\n',cc);
end

return
end

function dzWrite_specify_options(fid,specify)

% fix me
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];\n');
if ~isfield(specify,'volt')||isempty(specify.volt), specify.volt=1; end
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.volt = %d;\n',specify.volt);
if ~isfield(specify,'global')||isempty(specify.global), specify.global='None'; end
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.global = ''%s'';\n',specify.global);
if ~isfield(specify,'mask')||isempty(specify.mask), specify.mask=''; end
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.mask = {''%s''};\n',specify.mask);
if ~isfield(specify,'cvi')||isempty(specify.cvi), specify.cvi='AR(1)'; end
fprintf(fid,'matlabbatch{1}.spm.stats.fmri_spec.cvi = ''%s'';\n',specify.cvi);

return
end

function dzWrite_estimate_options(fid)

fprintf(fid,'matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;\n');
fprintf(fid,'matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = ''Select SPM.mat'';\n');
fprintf(fid,'matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = ''filter'';\n');
fprintf(fid,'matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = ''mat'';\n');
fprintf(fid,'matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = ''strtype'';');
fprintf(fid,'matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = ''e'';\n');
fprintf(fid,'matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = ''fMRI model specification: SPM.mat File'';\n');
fprintf(fid,'matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct(''.'',''val'', ''{}'',{1}, ''.'',''val'', ''{}'',{1}, ''.'',''val'', ''{}'',{1});\n');
fprintf(fid,'matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct(''.'',''spmmat'');\n');
fprintf(fid,'matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;\n');

return
end