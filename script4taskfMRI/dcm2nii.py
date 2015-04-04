
import sys
import os
import shutil

# config
merge = 0

dcm_Dir = '/home/houyuling/emotion_conflict/FunRaw'
nii_Dir = '/home/houyuling/emotion_conflict/FunImg'
dcm2nii_path = '/home/houyuling/emotion_conflict/script/dcm2nii'

# main
Dir = dcm_Dir 
if not os.path.exists(nii_Dir):
   os.mkdir(nii_Dir)
for subj in os.listdir(Dir):
     subj_path = Dir + os.sep + subj
     print subj_path
     subj_nii_dir = nii_Dir + os.sep + subj 
     if not os.path.isdir(subj_path):
        continue
     if not os.path.isdir(subj_nii_dir):
        os.mkdir(subj_nii_dir)
     for run in os.listdir(subj_path):
         run_dcm_dir = subj_path + os.sep + run
         dcm_file = run_dcm_dir + os.sep + os.listdir(run_dcm_dir)[0]
 
         # merge or not -------                 
         if merge==0:
            run_nii_dir = subj_nii_dir + os.sep + run
            if not os.path.exists(run_nii_dir):
               os.mkdir(run_nii_dir)
         else:
            run_nii_dir = subj_nii_dir
         # merge -------                 
 
         # unix ---
         dcm2nii_lnx_path = dcm2nii_path + os.sep + 'dcm2nii_linux'
         dcm2nii_lnx_ini_path = dcm2nii_path + os.sep + 'dcm2nii_linux_3DImg.ini'
         Option = ' -b '+ dcm2nii_lnx_ini_path
         cmd = 'chmod +x ' + dcm2nii_lnx_path
         os.popen(cmd)
         cmd = dcm2nii_lnx_path + Option + ' -o '+ run_nii_dir + ' ' +  dcm_file
         print os.popen(cmd).readlines()
         # unix ---
