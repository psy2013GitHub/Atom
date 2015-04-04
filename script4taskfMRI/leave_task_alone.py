
import sys
import os
import shutil

# config
test_Dir = "/home/houyuling/emotion_conflict/Raw/test"
elder_Dir = "/home/houyuling/emotion_conflict/Raw/elder"
control_Dir = "/home/houyuling/emotion_conflict/Raw/elder_control_201402"

# main
Dir = control_Dir
for subj in os.listdir(Dir):
     subj_path = Dir + os.sep + subj
     if not os.path.isdir(subj_path):
        continue
     for modality in os.listdir(subj_path):
         modality_path = subj_path+os.sep+modality
         if modality.rfind("task") != -1: # task
            for run in os.listdir(modality_path):
                run_path = modality_path + os.sep + run
                shutil.move(run_path,subj_path+os.sep+run)
            shutil.rmtree(modality_path)
            continue  
         else:
            shutil.rmtree(modality_path)
