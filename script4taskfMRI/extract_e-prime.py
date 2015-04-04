

# ! remmber E-prime use utf-18 by default in Chinese Windows
#

import sys
import os
import shutil
import codecs
open=codecs.open    # for usage of special encoding style


def match_subj_output_dir(OutputDir, subj):
    match_dir = []
    already_match = None
    for sub_dir in os.listdir(OutputDir):
        if os.isdir(sub_dir):
           if sub_dir.find(subj) != -1:
              if not already_match:
                 already_match = 1
              else:
                 print 'Error: ', subj, ' multi-match ',  
                 print '\t', match_dir
                 print '\t', OutputDir + os.sep + sub_dir
                 return None
              match_dir = OutputDir + os.sep + sub_dir
    return match_dir

def desing_matrix(Subj_lst, OutputDir):
    for subj, Exp_dict in Subj_lst.iteritems():
        subj_match_dir = match_subj_output_dir(OutputDir, subj):
        for Exp in Exp_dict.iterkeys():
            # init
            trial_lst      = Exp_dict[Exp]['Trial']
            ACC_lst        = Exp_dict[Exp]['ACC']
            RT_lst         = Exp_dict[Exp]['RT']
            Cond_lst       = Exp_dict[Exp]['Cond']
            Congruency_lst = Exp_dict[Exp]['Congruency']
            # Cond_category by hash dict
            cond_Category = {}
            for c in Cond_lst:
                cond_Category[c] = 1
            cond_Category = cond_Category.keys()
            # loop category
            subj_exp_design_mat = {}
            for c in cond_Category:
               if not c:
                  continue
               ACC = [0,1]
               for acc in ACC:
                   new_cond = c + acc
                   cond_trial_idx = [] 
                   subj_exp_design_mat 
                   for i, trial_c in enumerate(Cond_lst):
                       if ACC_lst[i] == acc:
                           cond_trial_idx.append(trial_lst[i]) 
                   cond_trial_idx.sort() 
                             
                 with open(OutputDir + subj)

def loop_direc(Dir):
    Subj_lst={} # dict
    for txt in os.listdir(Dir):
        txt_path = Dir + os.sep + txt

        idx = txt.rfind('.txt')
        if idx != -1:
           subj_dir = output_Dir +os.sep + txt[:idx]
           Subj, Exp, subj_trial_lst, subj_ACC_lst, subj_RT_lst, subj_Cond_lst, subj_Congruency_lst = extract_txt(txt_path)
           
           # cat all subjs
           Subj_lst.setdefault(Subj,{}) 
           Subj_lst[Subj].setdefault(Exp,{})

           Subj_lst[Subj][Exp].setdefault('Trial',[]);        Subj_lst[Subj][Exp]['Trial'].append(subj_trial_lst)
           Subj_lst[Subj][Exp].setdefault('ACC',[]);          Subj_lst[Subj][Exp]['ACC'].append(subj_ACC_lst)
           Subj_lst[Subj][Exp].setdefault('RT',[]);           Subj_lst[Subj][Exp]['RT'].append(subj_RT_lst)
           Subj_lst[Subj][Exp].setdefault('Cond',[]);         Subj_lst[Subj][Exp]['Cond'].append(subj_Cond_lst)
           Subj_lst[Subj][Exp].setdefault('Congruency',[]);   Subj_lst[Subj][Exp]['Congruency'].append(subj_Congruency_lst)

        else:
           continue

    return Subj_lst

def extract_txt(txt):
    fid = open(txt, 'r', 'utf-16') 

    Task_lst = []; ACC_lst = []; RT_lst = []; Cond_lst = []; Congruency_lst = []
    out_of_trial = 1
    iter = 0
    for line in fid:
        iter += 1
        line = line.lstrip() # delete left blank
        #print line.find('*** LogFrame Start ***')
        if line.find('*** LogFrame Start ***') != -1:
           out_of_trial = 0 # now in    
        if line.find('*** LogFrame End ***') != -1:
           out_of_trial = 1 # now out    

        if out_of_trial == 0: # in trial
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
              Exp = str(line[exp_idx+len(k_str):].rstrip())
              continue

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

    return Subj, Exp, Task_lst, ACC_lst, RT_lst, Cond_lst, Congruency_lst,

if __name__ == "__main__":
   txt_Dir = "/home/houyuling/emotion_conflict/behav/test"
   output_Dir = "/home/houyuling/emotion_conflict/behav/test_out"
   Dir = txt_Dir
   Subj_lst = loop_direc(Dir)
   print Subj_lst
