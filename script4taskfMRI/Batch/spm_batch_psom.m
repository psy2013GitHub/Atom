function spm_batch_psom(nTask,DataDir,Steps,path_logs)
% 
%   e,g, spm_batch_psom(8,'/home/...',{'',''},'/home/...')

% further_DataDir, used to note whether change DataDir

if exist(path_logs, 'file')
    rmdir(path_logs, 's');
end

run spm_pipeline_init;

% config setting
pipeline_opt.path_logs=path_logs;
pipeline_opt.mode='batch';
pipeline_opt.mode_pipeline_manager='session';
pipeline_opt.max_queued=nTask;
pipeline_opt.nb_resub=0;
pipeline_opt.flag_update=1;
pipeline_opt.flag_verbose=1;
pipeline_opt.flag_debug=1;
pipeline_opt.flag_pause=0;


for step_idx = 1:length(Steps)
    display(step_idx)
    step = Steps{step_idx};
    
    
    clear further_DataDir; % each step each set
    Subj_lst = getSubj_lst(DataDir);
    nSubj = length(Subj_lst);
    nsubj_perTask = floor(nSubj/nTask);
    
    clear opt;
    opt.spm_pipeline  = spm_pipeline;
    
    opt.root = DataDir;
    %display(DataDir)
    for ii=1:nTask                                      % search in mask space
        
        prefix = step;
        
        %- Job name
        tmp_job_name=sprintf('psom_%sTask_%d',prefix,ii);
        
        %- Command
        subj_left_idx=(ii-1)*nsubj_perTask+1;
        if ii~=nTask
            subj_right_idx=ii*nsubj_perTask;
        else
            subj_right_idx=nSubj;
        end
 
        %- Specific steps 
        opt.Subj_lst = Subj_lst(subj_left_idx : subj_right_idx);

        if strcmpi(step,'Realign')
            command = ...
                'Batch_Realign_func(opt.root, opt.Subj_lst, opt.nii_filter, opt.spm_pipeline.Batch_Realign_Dir, opt.MultiRun_flag, opt.nDiscard_pre, opt.nDiscard_end, opt.eoptions, opt.roptions);';
            
            opt.MultiRun_flag = 1;
            opt.nii_filter = '^a.*\.(img|nii)$';
            opt.nDiscard_pre=5; opt.nDiscard_end=0;
            opt.eoptions.rtm=1; opt.roptions=struct();
        
        elseif strcmpi(step, 'Normalize')
            command = ...
                'Batch_Normalise_func(opt.root, opt.Subj_lst, opt.nii_filter, opt.MultiRun_flag, opt.spm_pipeline.Batch_Normalise_Dir, opt.eoptions, opt.roptions);';
            
            opt.MultiRun_flag = 1;
            opt.nii_filter.img='^ra.*\.(img|nii)$';
            opt.nii_filter.src='^mean.*\.(img|nii)$';
            opt.roptions.vox=[3,3,3];
            opt.eoptions=struct();
            
        elseif strcmpi(step, 'Smooth')
            command = ...
                'Batch_Smooth_func(opt.root, opt.Subj_lst, opt.nii_filter, opt.spm_pipeline.Batch_Smooth_Dir, opt.MultiRun_flag, opt.smooth);';
            
            opt.MultiRun_flag = 1;
            opt.nii_filter='^wr.*\.(img|nii)$';
            opt.smooth.fwhm=[8,8,8];
            
        elseif strcmpi(step, '1stlevel')
            command = ...
                'Batch_1stlevel_func(opt.root, opt.Subj_lst, opt.ResultDir, opt.nii_filter, opt.spm_pipeline.Batch_1stlevel_Dir, opt.MultiRun_flag, opt.MergeRun_flag, opt.specify, opt.subj_specific_covariates);';
            
            opt.ResultDir = '/home/houyuling/emotion_conflict/GLM_1stlevel';
            opt.nii_filter='^swr.*\.(img|nii)$';
            opt.MultiRun_flag=1;
            opt.MergeRun_flag=1;
            opt.subj_specific_covariates=1;
            
            %- model parameter specific
            opt.specify.timing.units='secs'; % 'scans' or 'secs'
            opt.specify.timing.RT=2;
            opt.specify.timing.fmri_t=32;
            opt.specify.timing.fmri_t0=16;
            opt.specify.sess.hpf = inf;      % default 128

            %- design specific
            %condition 
            opt.specify.sess.cond(1).name='mid_Rest';
            opt.specify.sess.cond(1).onset=201;
            opt.specify.sess.cond(1).duration=4;
            
            opt.specify.sess.cond(2).name='Fixation';
            opt.specify.sess.cond(2).onset=[
                0,3,8,12,16,22,26,31,36,42,47,51,55,61,67,72,76,82,87,91,97,101,105,110,116,121,127,132,136,141,147,152,156,162,168,173,177,183,188,192,196,205,210,214,220,226,231,235,240,246,250,256,260,264,269,274,278,283,289,295,299,304,309,315,319,325,330,335,341,345,349,354,360,366,370,375,381,386,390,394,400,405
                ];
            
            opt.specify.sess.cond(2).duration=[
                2,4,3,3,5,3,4,4,5,4,3,3,5,5,4,3,5,4,3,5,3,3,4,5,4,5,4,3,4,5,4,3,5,5,4,3,5,4,3,3,5,4,3,5,5,4,3,4,5,3,5,3,3,4,4,3,4,5,5,3,4,4,5,3,5,4,4,5,3,3,4,5,5,3,4,5,4,3,3,5,4,5
                ];
            %regressor
            opt.specify.regress_motion=1;
            %this multi-run regressor added as you wish
%             total_nImg_lst = [210,210]; 
%             nImgDiscard_pre = 5; 
%             nImgDiscard_pro = 0;
%             regressor_cell = {'linear_trend' , 'nonlinear_fourier_trend' , 'cofound_mean'}; % 'linear_trend' | 'nonlinear_fourier_trend' | 'cofound_mean'
%             
%             run_specific_regressor = design_specific_regressor(total_nImg_lst, nImgDiscard_pre, nImgDiscard_pro, opt.specify.timing.RT, regressor_cell);
%             if ~isempty(run_specific_regressor)
%                 opt.specify.sess.regress = run_specific_regressor;
%             end
            
            % !!!! chang root, change by yourself as needed
            further_DataDir = opt.ResultDir;
            
            
        elseif strcmpi(step, 'Contrast')
            
            command = ...
                'Batch_Contrast_func(opt.root, opt.Subj_lst, opt.spm_pipeline.Batch_Contrast_Dir, opt.MultiRun_flag, opt.contrast);';
            
            % !!! if merged, set to 0
            opt.MultiRun_flag = 0;
            
            opt.contrast.consess{1}.tcon.name='incon1 - con1';
            opt.contrast.consess{1}.tcon.cond_list={'incon_1','con_1'};
            opt.contrast.consess{1}.tcon.cond_weight=[1,-1];
            
            opt.contrast.consess{2}.tcon.name='incon1 - Fixation';
            opt.contrast.consess{2}.tcon.cond_list={'incon_1','Fixation'};
            opt.contrast.consess{2}.tcon.cond_weight=[1,-1];
            
            opt.contrast.consess{3}.tcon.name='con1 - Fixation';
            opt.contrast.consess{3}.tcon.cond_list={'con_1','Fixation'};
            opt.contrast.consess{3}.tcon.cond_weight=[1,-1];
            
            
        elseif strcmpi(step, 'mvContrast')
            
            command = ...
                'for kk=1:length(opt.keyname), Move_con(opt.root, opt.Subj_lst, opt.ResultDir, opt.filter, opt.flag, opt.keyname{kk}); end;';
            
            opt.ResultDir='/home/houyuling/emotion_conflict/GLM_2ndlevel';
            opt.filter='^con_(0001|0002|0003).*\.(img|hdr|nii)$'; % '^con_(0002|0003).*\.(img|nii)$'; !!! 'img/hdr' paire
            opt.flag='dirname';
            opt.keyname=[{'elder'}]; % [{'unr'},{'suc'},{'fail'}];
            
            % !!!! chang root, change by yourself as needed
            further_DataDir = opt.ResultDir;
            
        else
            error('Unknown %s', step);
            
        end
        
        
        %- General Settings
        fileout = sprintf('%s%s%sTask_%d',DataDir,filesep,step,ii);
        command = [command,create_step_done_file(fileout)];
        pipeline.(tmp_job_name).command=command;
        
        %files_in & files_out
        pipeline.(tmp_job_name).files_out= fileout;
        if step_idx > 1
            last_step = Steps{step_idx-1};
            pipeline.(tmp_job_name).files_in = sprintf('%s%s%sTask_%d',DataDir,filesep,last_step,ii);
        end
        
        %opt
        pipeline.(tmp_job_name).opt= opt;
        
    end
    
    % step after step
    psom_run_pipeline(pipeline, pipeline_opt);
    
    clear pipeline
    
    if exist('further_DataDir', 'var')
        DataDir = further_DataDir;
    end
    
end

% cleanup immediate file
% for ii=1:nTask
%     tmp_job_name=sprintf('%sTask_%d',prefix,ii);
%     pipeline.cleanup.files_clean{ii}=pipeline.(tmp_job_name).files_out;
% end
% command='delete(files_clean{:})';
% pipeline.cleanup.command=command;

disp('---------------------------------------');

return
end


function Subj_lst = getSubj_lst(DataDir)

Subj_lst = dir(DataDir);
Subj_lst = Subj_lst(3:end);
Subj_lst = Subj_lst([Subj_lst(:).isdir] == 1);

return
end

function cmd_str = create_step_done_file(fileout)

str1 = sprintf('fid=fopen(''%s'',''w'');',fileout);

cmd_str = ... 
    [str1,...
     'fprintf(fid,''%% Created at %s %s\n'',now, date);',...
     'fclose(fid);'];

return
end