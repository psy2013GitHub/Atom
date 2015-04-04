
function Batch_Contrast_func(root, Subj_lst, Batch_Contrast_Dir, MultiRun_flag, contrast)

for ss=1:length(Subj_lst)
    tmp_subj=Subj_lst(ss);
    if ~tmp_subj.isdir, continue; end
    SubjDir=[root, filesep, Subj_lst(ss).name];
    SingleSubj_Contrast_func(SubjDir, Batch_Contrast_Dir, MultiRun_flag, contrast);
end
return
end
