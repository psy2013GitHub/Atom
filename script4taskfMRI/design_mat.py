
#!coding=utf8
import os
import shutil

import e_prime_utility


def design_matrix(Subj_lst, OutputDir, Exp_order, Exp_trial_onset_time, Exp_cond_duration, Exp_total_time, Exp_merge=True):
    '''
     Subj_lst: src txt directory
     OutpuDir: des img directory
     Exp_order: experiment(whole name!) order, like ['run1','run2']
     Exp_trial_onset_time: onset of each recorded trial of each Exp in the order of 'Exp_order'
     Exp_cond_duration: dict in the form of {'run1':{'cond1':,'cond2':},'run2':{...}}, duration of each recorded trial of each Exp 
     cond_duration_time: duration of each stimulus, all the same if scaler
     Exp_merge: merge experiment or not
     Exp_total_time: list type with each element refers to the total time of each Exp 
      
    '''
    nExp = len(Exp_order)
    if nExp != len(Exp_trial_onset_time) and Exp != len(Exp_total_time):
       raise typeerror('unequal length')
   
    # if merge Exps:
    if nExp > 1 and Exp_merge:
       cum_time = 0
       for e in xrange(1, len(Exp_total_time)):
           cum_time += Exp_total_time[e-1]
           Exp_trial_onset_time[e] = [i+cum_time for i in Exp_trial_onset_time[e]]
       
    # iterate over subjects
    for subj in Subj_lst.keys():
        Exp_dict = Subj_lst[subj]
        subj_output_dir = Output_subj_dir_match(OutputDir, str(subj))
        if not subj_output_dir:
           print subj, ' not found'
           continue
        else:
           print subj, '\t', subj_output_dir
        subj_covariates_output_dir =subj_output_dir +  os.sep + 'covariates'

        if os.path.exists(subj_covariates_output_dir):
           shutil.rmtree(subj_covariates_output_dir)
        os.mkdir(subj_covariates_output_dir)
 
        cond_onset_dict = {}; cond_duration_dict = {}; 
        file_cond_lst = [] # cond name which will be written into file
        for Exp_i, Exp in enumerate(Exp_order):
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
            for c in Exp_cond_duration[Exp].keys(): # raw conditon
                exp_cond_duration = Exp_cond_duration[Exp][c]
                if not c:
                   continue
                ACC = [0,1] # interact condition with ACC
                for acc in ACC:
                    new_cond = c + '_' + str(acc)
                    cond_onset_dict.setdefault(new_cond, [])  # default as []
                    cond_duration_dict.setdefault(new_cond, [])  # default as []
                    cond_trial_idx = []
                    for i, trial_c in enumerate(Cond_lst):
                        if Cond_lst[i] == c and ACC_lst[i] == acc:
                           cond_trial_idx.append(trial_lst[i])

                    cond_trial_idx.sort()
                    if Exp_merge: # fix me in future add term 'Exp_cond_offset_from_trial_onset'
                       cond_onset_dict[new_cond] += [Exp_trial_onset_time[Exp_i][i-1] for i in cond_trial_idx] # !!!notion!, E-prime starts from 1, while python starts from 0
                       if len(exp_cond_duration) == 1:
                          cond_duration_dict[new_cond] += [exp_cond_duration[0] for i in xrange(len(cond_trial_idx))]# repeation if merge
                    else:
                       raise typeerror('unsupported yet')

            # write to file
            # _onset &  _duration
            fid_cond = open(subj_covariates_output_dir + os.sep + '_cond_name', 'w')
            for c in cond_onset_dict.keys():
               # cond
               fid_cond.writelines('%s \n' % c)
               # onset
               fid = open(subj_covariates_output_dir + os.sep + c + '_onset', 'w')
               fid.writelines('%s \n' % i for i in cond_onset_dict[c])
               fid.close()
               # duration
               fid = open(subj_covariates_output_dir + os.sep + c + '_duration', 'w')
               fid.writelines('%s \n' % i for i in cond_duration_dict[c])
               fid.close()
            fid_cond.close()

           
def Output_subj_dir_match(OutputDir, subj):
    '''
      search the directory of subj in OutputDir
    '''
    already_match = 0
    for sub_dir in os.listdir(OutputDir):
        if os.path.isdir(OutputDir + os.sep + sub_dir) and sub_dir.find(subj) != -1:
           if already_match == 0:
              match_dir = OutputDir + os.sep + sub_dir
              already_match = 1
           else:
              print 'Error: ', subj, ' multi-match'
              print '\t ', match_dir
              print '\t ', OutputDir + os.sep + sub_dir 
              raise Typeeror
    if already_match == 0:
       return None
    else:
       return match_dir 


if __name__ == "__main__":
   # *** parse e-prime txt *** #
   txt_Dir = "/home/houyuling/emotion_conflict/behav2"
   Subj_lst = e_prime_utility.parser(txt_Dir).parse()
   #print Subj_lst
   # *** design_matrix     *** #
   # ImgDir
   OutputDir =  "/home/houyuling/emotion_conflict/FunImg2"
   # Exp_order
   Exp_order = ['Emotional_conflict_run_1','Emotional_conflict_run_2']
   # Exp_trial_onset_time, in the order of Exp
   Exp_trial_onset_time = [[2,7,11,15,21,25,30,35,41,46,50,54,60,66,71,75,81,86,90,96,100,104,109,115,120,126,131,135,140,146,151,155,161,167,172,176,182,187,191,195,201,209,213,219,225,230,234,239,245,249,255,259,263,268,273,277,282,288,294,298,303,308,314,318,324,329,334,340,344,348,353,359,365,369,374,380,385,389,393,399,404]]
   # Exp_cond_duration, !!!caution, condition must be list
   Exp_cond_duration = {'Emotional_conflict_run_1':{'con':[1],'incon':[1]},'Emotional_conflict_run_2':{'con':[1],'incon':[1]}}
   # Exp_total_time
   Exp_total_time = [410, 410] 
   design_matrix(Subj_lst, OutputDir, Exp_order, Exp_trial_onset_time, Exp_cond_duration, Exp_total_time, Exp_merge=True)
