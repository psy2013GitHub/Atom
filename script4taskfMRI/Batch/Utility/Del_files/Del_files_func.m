function Del_files_func(FileDir,FileFilter,flag)
cwd=pwd;
subdir_list=dir(FileDir); subdir_list=subdir_list(3:end);
for ss=1:length(subdir_list)
    tmp_subdir=subdir_list(ss); if ~tmp_subdir.isdir, continue; end
    cd([FileDir,filesep,tmp_subdir.name]);
    cnt=[];
    Files=spm_select('List','.',FileFilter);
    if strcmpi(flag,'delete'), 
        if isempty(Files), fprintf('%d files deleted in %s\n',size(Files,1), tmp_subdir.name); continue; end
        for ff=1:size(Files,1),
            try
                delete(Files(ff,:));
                fprintf('%d files deleted in %s\n',size(Files,1), tmp_subdir.name);
            catch
                [msgstr, msgid] = lastwarn;
                warning('off', msgid);
            end
        end
        % check
        Files=spm_select('List','.',FileFilter);
        cnt=[cnt,size(Files,1)];
    end
    if strcmpi(flag,'check'), fprintf('%d files found in %s\n',size(Files,1),tmp_subdir.name); end
end
while strcmpi(flag,'delete')&&any(cnt~=0), cd(cwd); Del_files_func(FileDir,FileFilter,flag); end
cd(cwd);
fprintf('\n\n');
return
end