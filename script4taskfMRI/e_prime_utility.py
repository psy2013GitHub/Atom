

# ! remmber E-prime use utf-18 by default in Chinese Windows
#

import sys
import os
import shutil
import codecs
open=codecs.open  # for usage of special encoding style

class parser(object):
     '''
      'ACC/RT = -1' means ACC/RT unexist 
     '''

     def __init__(self, input):
         self.input = input
     
     def parse(self):
         if os.path.isdir(self.input):
            Subj_lst = self.loop_direc(self.input)
            return Subj_lst 
         else:
            self.extract_txt(self.input) # fix me
     
     def loop_direc(self, Dir):
         Subj_lst={} # dict
         for txt in os.listdir(Dir):
             txt_path = Dir + os.sep + txt

             r_idx = len(txt) - txt.rfind('.txt')
             if r_idx == 4:
                Subj, Exp, subj_trial_lst, subj_ACC_lst, subj_RT_lst, subj_Cond_lst, subj_Congruency_lst = self.extract_txt(txt_path)
      
                # cat all subjs
                Subj_lst.setdefault(Subj,{}) 
                Subj_lst[Subj].setdefault(Exp,{})

                Subj_lst[Subj][Exp].setdefault('Trial',subj_trial_lst);
                Subj_lst[Subj][Exp].setdefault('ACC',subj_ACC_lst);
                Subj_lst[Subj][Exp].setdefault('RT',subj_RT_lst);
                Subj_lst[Subj][Exp].setdefault('Cond',subj_Cond_lst);
                Subj_lst[Subj][Exp].setdefault('Congruency',subj_Congruency_lst);

             else:
                continue

         return Subj_lst

     def extract_txt(self, txt):
         fid = open(txt, 'r', 'utf-16') 

         Task_lst = []; ACC_lst = []; RT_lst = []; Cond_lst = []; Congruency_lst = []
         out_of_trial = 1
         for line in fid:
             line = line.lstrip() # delete left blank

             # Subject
             k_str = 'Subject:'
             subj_idx = line.find(k_str)
             if subj_idx != -1:
                Subj = int(line[subj_idx+len(k_str):])
                continue

             # Exp
             k_str = 'Experiment:'
             exp_idx = line.find(k_str)
             if exp_idx != -1:
                Exp = str(line[exp_idx+len(k_str):].lstrip().rstrip())
                continue

             if line.find('*** LogFrame Start ***') != -1:
                while line.find('*** LogFrame End ***') == -1:

                      line = fid.next()
                      line = line.lstrip() # delete left blank
                      # trial idx
                      k_str = 'task:'
                      trial_idx = line.find(k_str)
                      if trial_idx != -1:
                         trial_idx = int(line[trial_idx+len(k_str):])
                         Task_lst.append(trial_idx)
                         continue
  
                      # cond
                      k_str = 'Condition:'
                      cond_idx = line.find(k_str)
                      if cond_idx != -1:
                         cond = str(line[cond_idx+len(k_str):].lstrip().rstrip())
                         Cond_lst.append(cond)
                         continue
      
                      # congruency
                      k_str = 'congruency:'
                      congruency_idx = line.find(k_str)
                      if congruency_idx != -1:
                         congruency = str(line[congruency_idx+len(k_str):].lstrip().rstrip())
                         Congruency_lst.append(congruency)
                         continue

                      # stim.ACC / RT
                      stim_idx = line.find('stim.')
                      if stim_idx != -1:
                         # ACC
                         k_str = 'ACC:'
                         acc_idx = line.find(k_str)
                         if acc_idx != -1:
                            acc = int(line[acc_idx+len(k_str):].lstrip())
                            ACC_lst.append(acc)
                            continue
                         # RT
                         k_str = 'RT:'
                         rt_idx = line.find(k_str)
                         if rt_idx != -1:
                            rt = float(line[rt_idx+len(k_str):].lstrip())
                            RT_lst.append(rt)
                            continue

                # for special task
                if len(Task_lst) > len(ACC_lst):
                   acc = -1; ACC_lst.append(acc)
                if len(Task_lst) > len(RT_lst):
                   rt = -1;  RT_lst.append(rt)
                if len(Task_lst) > len(Congruency_lst):
                   congruency = 'NaN'; Congruency_lst.append(congruency)
                if len(Task_lst) > len(Cond_lst):
                   cond = 'NaN'; Cond_lst.append(cond)
                   
         return Subj, Exp, Task_lst, ACC_lst, RT_lst, Cond_lst, Congruency_lst,

if __name__ == "__main__":
  txt_Dir = "/home/houyuling/emotion_conflict/behav/test"
  Subj_lst = parser(txt_Dir).parse()
  print Subj_lst
