
function Batch_ResultDisp_func(DataDir,DestDir,Batch_ResultDisp_Dir,MultiRuns_flag, conspec)

cwd=pwd; missps_cnt=0;
if exist(DestDir,'dir')~=7, mkdir(DestDir); end
subj_list=dir(DataDir); subj_list=subj_list(3:end);
for ss=1:length(subj_list)
    tmp_subj=subj_list(ss);
    if ~tmp_subj.isdir, continue; end
    log_file=[Batch_ResultDisp_Dir,filesep,'resultdisp_log.txt'];
    log_fid=fopen(log_file,'w');dz_Write_batch_header(log_fid)
    SubjDir=[DataDir,filesep,tmp_subj.name];
    exist_ps_flag=SingleSubj_ResutDisp_func(SubjDir, DestDir, Batch_ResultDisp_Dir,MultiRuns_flag, conspec);
    cd(cwd);
    if ~exist_ps_flag, dz_Write_batch_missps(log_fid,missps_cnt,tmp_subj.name);  missps_cnt=missps_cnt+1; continue; end
end
delete .*asv;
return
end

function dz_Write_batch_header(fid)

fprintf(fid,'\n\n%%----------------------------------------------------\n');
fprintf(fid,'%% Created at %s\n',datestr(now,'yyyymmmdd HH:MM:SS'));
fprintf(fid,'%%----------------------------------------------------\n\n\n');

return
end

function dz_Write_batch_missps(fid,missps_cnt,Subj)

if missps_cnt==0, fprintf(fid,'Subj without SPM*.ps file:'); end;
fprintf(fid,'                          %s\n',Subj);
return
end

