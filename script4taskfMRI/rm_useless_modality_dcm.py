
import sys
import os
import shutil

# config
elder_Dir = "/home/houyuling/emotion_conflict/Raw/elder"
control_Dir = "/home/houyuling/emotion_conflict/Raw/elder_control_201402"

# main
Dir = control_Dir
for subj in os.listdir(Dir):
  task_run_lst = []
  if os.path.isdir(Dir+os.sep+subj):
     for modality in os.listdir(Dir+os.sep+subj):
         modality_path = Dir+os.sep+subj+os.sep+modality
         if modality.rfind("task") != -1: # task
            os.rename(modality_path, Dir+os.sep+subj+os.sep+"task")
            continue  
         if os.path.isdir(modality_path):
            n_dcm = len(os.listdir(modality_path)) 
         if n_dcm == 242:
            os.rename(modality_path, Dir+os.sep+subj+os.sep+"rest")
         elif n_dcm == 176: 
            os.rename(modality_path, Dir+os.sep+subj+os.sep+"sMRI")
         elif n_dcm == 62:
            os.rename(modality_path, Dir+os.sep+subj+os.sep+"dMRI")
         elif n_dcm % 210 == 0: # task
            task_run_lst.append(modality_path)
         else:
            shutil.rmtree(modality_path)
  if task_run_lst:
     task_dir = Dir+os.sep+subj+os.sep+"task"
     if not os.path.exists(task_dir):
        os.mkdir(task_dir)
     for run in task_run_lst:
        run_p, run_n = os.path.split(run)
        shutil.move(run,task_dir+os.sep+run_n)
       
