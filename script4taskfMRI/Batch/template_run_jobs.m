
function template_run_jobs(jobfile)

% List of open inputs
nrun = 1; % enter the number of runs here
%jobfile = {'/home/houyuling/emotion_conflict/script/Batch/Contrast/tmp_contrast_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('serial', jobs, '', inputs{:});

end