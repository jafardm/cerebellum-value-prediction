%% Bulid behave dataset
function JDM_bulid_behave_dataset_session

JDM_params_funcs
user = 'JDM';
loc  = 'ctx';
tags = 1:10;

path_data_monkey_sorted = params.path_data_monkey_sorted;
animal_list = params.animal_list;
session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);

num_animal = numel(session_list);

data_behav.sac_amp       = cell(num_animal,1);
data_behav.sac_vel       = cell(num_animal,1);
data_behav.sac_tag       = cell(num_animal,1);
data_behav.sac_ang       = cell(num_animal,1);
data_behav.sac_acc_dur   = cell(num_animal,1);
data_behav.sac_dec_dur   = cell(num_animal,1);
data_behav.sac_acc       = cell(num_animal,1);
data_behav.sac_dec       = cell(num_animal,1);
data_behav.sess_cr_ms    = cell(num_animal,1);
data_behav.sess_cr_ma    = cell(num_animal,1);
data_behav.sess_id       = cell(num_animal,1);
data_behav.sac_vec       = cell(num_animal,1);
data_behav.gl_vec        = cell(num_animal,1);
data_behav.sac_rt        = cell(num_animal,1);
data_behav.tgt_cond      = cell(num_animal,1);
data_behav.tgt_num       = cell(num_animal,1);        
data_behav.cue_x_high    = cell(num_animal,1); 
data_behav.cue_x_low     = cell(num_animal,1);
data_behav.cue_y_high    = cell(num_animal,1);
data_behav.cue_y_low     = cell(num_animal,1);  
data_behav.high_tgt_num  = cell(num_animal,1);
data_behav.low_tgt_num   = cell(num_animal,1);
data_behav.task_cond     = cell(num_animal,1);       
data_behav.choice        = cell(num_animal,1);  
data_behav.rew_cond      = cell(num_animal,1);     
data_behav.jump_cond     = cell(num_animal,1);
data_behav.eye_vm_vis    = cell(num_animal,1);
data_behav.eye_vm_sac    = cell(num_animal,1);

for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess = numel(current_sess_list);

    sac_amp_sess       = cell(num_sess,1);
    sac_vel_sess       = cell(num_sess,1);
    sac_tag_sess       = cell(num_sess,1);
    sac_ang_sess       = cell(num_sess,1);
    sac_acc_dur_sess   = cell(num_sess,1);
    sac_dec_dur_sess   = cell(num_sess,1);
    sac_acc_sess       = cell(num_sess,1);
    sac_dec_sess       = cell(num_sess,1);
    sess_criteria_ms   = cell(num_sess,1);
    sess_criteria_ma   = cell(num_sess,1);
    sess_id_sess       = cell(num_sess,1);
    sac_vec_sess       = cell(num_sess,1);
    gl_vec_sess        = cell(num_sess,1);
    sac_rt_sess        = cell(num_sess,1);
    tgt_num_sess       = cell(num_sess,1);
    tgt_cond_sess      = cell(num_sess,1);
    cue_x_high_rew_sess   = cell(num_sess,1);
    cue_x_low_rew_sess    = cell(num_sess,1);
    cue_y_high_rew_sess   = cell(num_sess,1);
    cue_y_low_rew_sess    = cell(num_sess,1);
    high_rew_tgt_num_sess = cell(num_sess,1);
    low_rew_tgt_num_sess  = cell(num_sess,1);
    task_cond_sess        = cell(num_sess,1);
    choice_sess           = cell(num_sess,1);
    rew_cond_sess         = cell(num_sess,1);
    jump_cond_sess        = cell(num_sess,1);
    eye_vm_sess_vis       = cell(num_sess,1);
    eye_vm_sess_sac       = cell(num_sess,1);

   parfor counter_sess = 1:num_sess
        current_sess = current_sess_list{counter_sess};
        disp(['     ' current_sess]);

        data_eye             = MAF_load_eye_traces(current_path,current_sess);
        units_info           = MAF_extract_cell_metadata(current_path,current_sess);
        cell_num_rec         = sum(cell2mat(cellfun(@isempty,units_info.cell_ids,'UniformOutput',false)),2);
        [~,selected_cell_id] = min(cell_num_rec);
        data                 = MAF_load_cell(current_path,current_sess,...
            units_info.cell_list{selected_cell_id});
        rec_flag = data.rec_info.rec_flag;
        eye_traces = data_eye.eye_traces(rec_flag);
        sacs_tr    = [eye_traces.sac];
        tags_tr    = [sacs_tr.tag];
        res_am = MAF_find_max_acc_dec(eye_traces,tags);

        data_rec     = data.data_recordings;
        eye          = [data_rec.eye];
        sac          = [eye.sac];
        sac_tag_     = [sac.tag];
        sac_amp_     = [sac.eye_amp];
        sac_vel_     = [sac.eye_vm_max];
        sac_ang_     = [sac.eye_ang];

        tgt_num_     = [sac.tgt_num];
        tgt_cond_    = [sac.tgt_cond];
        cue_x_high_rew_   = [sac.cue_x_high_rew];
        cue_x_low_rew_    = [sac.cue_x_low_rew];
        cue_y_high_rew_   = [sac.cue_y_high_rew];
        cue_y_low_rew_    = [sac.cue_y_low_rew];
        high_rew_tgt_num_ = [sac.high_rew_tgt_num];
        low_rew_tgt_num_  = [sac.low_rew_tgt_num];
        task_cond_        = [sac.task_cond];
        choice_           = [sac.choice];
        rew_cond_         = [sac.rew_cond];
        jump_cond_        = [sac.jump_cond];

        sac_rt_      = ([sac.time_onset] - [sac.time_visual])*1e3;
        sac_acc_dur_ = ([sac.time_vmax] - [sac.time_onset])*1e3;
        sac_dec_dur_ = ([sac.time_offset] - [sac.time_vmax])*1e3;

        assert(all(tags_tr == sac_tag_));

        sac_end_x_ = [sac.eye_px_offset];
        sac_end_y_ = [sac.eye_py_offset];

        sac_st_x_  = [sac.eye_px_onset];
        sac_st_y_  = [sac.eye_py_onset];

        vis_end_x_ = [sac.vis_px_offset];
        vis_end_y_ = [sac.vis_py_offset];

        sac_vec_   = (sac_end_x_ - sac_st_x_) + 1j*(sac_end_y_ - sac_st_y_);
        gl_vec_    = (vis_end_x_ - sac_st_x_) + 1j*(vis_end_y_ - sac_st_y_);

        res = MAF_find_sac_outliers(sac);

        sess_id_sess{counter_sess}     = current_sess;
        sess_criteria_ms{counter_sess} = res.M_mainseq;
        sess_criteria_ma{counter_sess} = res.M_accdec;

        ind_selected     = ismember(sac_tag_,tags) & ~res.rm_ind; % not outlies + tag 1-10
        ind_selected_tag = ~res.rm_ind(ismember(sac_tag_,tags));
 
        eye_vm = MAF_combine_vm_traces(eye_traces,{'visual','onset'},250,tags);
        eye_vm_vis = eye_vm.visual.eye_vm;
        eye_vm_sac = eye_vm.onset.eye_vm;


        sac_amp_sess{counter_sess}     = sac_amp_(ind_selected)';
        sac_vel_sess{counter_sess}     = sac_vel_(ind_selected)';
        sac_tag_sess{counter_sess}     = sac_tag_(ind_selected)';
        sac_ang_sess{counter_sess}     = sac_ang_(ind_selected)';
        sac_rt_sess{counter_sess}      = sac_rt_(ind_selected)';
        sac_acc_dur_sess{counter_sess} = sac_acc_dur_(ind_selected)';
        sac_dec_dur_sess{counter_sess} = sac_dec_dur_(ind_selected)';
        sac_acc_sess{counter_sess}     = res_am.am_max(ind_selected_tag);
        sac_dec_sess{counter_sess}     = res_am.dm_max(ind_selected_tag);
        sac_vec_sess{counter_sess}     = sac_vec_(ind_selected)';
        gl_vec_sess{counter_sess}      = gl_vec_(ind_selected)';

        tgt_num_sess {counter_sess}            = tgt_num_(ind_selected)';
        tgt_cond_sess {counter_sess}            = tgt_cond_(ind_selected)';
        cue_x_high_rew_sess {counter_sess}      = cue_x_high_rew_(ind_selected)';
        cue_x_low_rew_sess {counter_sess}       = cue_x_low_rew_(ind_selected)';
        cue_y_high_rew_sess {counter_sess}     = cue_y_high_rew_(ind_selected)';
        cue_y_low_rew_sess {counter_sess}       = cue_y_low_rew_(ind_selected)';
        high_rew_tgt_num_sess {counter_sess}    = high_rew_tgt_num_(ind_selected)';
        low_rew_tgt_num_sess {counter_sess}     = low_rew_tgt_num_(ind_selected)';
        task_cond_sess {counter_sess}           = task_cond_(ind_selected)';
        choice_sess {counter_sess}              = choice_(ind_selected)';
        rew_cond_sess {counter_sess}            = rew_cond_(ind_selected)';
        jump_cond_sess {counter_sess}           = jump_cond_(ind_selected)';
        eye_vm_sess_vis {counter_sess}          = eye_vm_vis;
        eye_vm_sess_sac {counter_sess}          = eye_vm_sac;
    end
    data_behav.sac_amp{counter_animal}     = sac_amp_sess;
    data_behav.sac_vel{counter_animal}     = sac_vel_sess;
    data_behav.sac_ang{counter_animal}     = sac_ang_sess;
    data_behav.sac_tag{counter_animal}     = sac_tag_sess;
    data_behav.sac_acc{counter_animal}     = sac_acc_sess;
    data_behav.sac_dec{counter_animal}     = sac_dec_sess;
    data_behav.sac_rt{counter_animal}      = sac_rt_sess;
    data_behav.sac_acc_dur{counter_animal} = sac_acc_dur_sess;
    data_behav.sac_dec_dur{counter_animal} = sac_dec_dur_sess;
    data_behav.sac_vec{counter_animal}     = sac_vec_sess;
    data_behav.gl_vec{counter_animal}      = gl_vec_sess;
    data_behav.sess_cr_ms{counter_animal}  = sess_criteria_ms;
    data_behav.sess_cr_ma{counter_animal}  = sess_criteria_ma;
    data_behav.sess_id{counter_animal}     = sess_id_sess;

    data_behav.tgt_num{counter_animal}      = tgt_num_sess;      
    data_behav.tgt_cond{counter_animal}     = tgt_cond_sess;        
    data_behav.cue_x_high{counter_animal}   = cue_x_high_rew_sess;   
    data_behav.cue_x_low{counter_animal}    = cue_x_low_rew_sess;   
    data_behav.cue_y_high{counter_animal}   = cue_y_high_rew_sess;  
    data_behav.cue_y_low{counter_animal}    = cue_y_low_rew_sess;  
    data_behav.high_tgt_num{counter_animal} = high_rew_tgt_num_sess;  
    data_behav.low_tgt_num{counter_animal}  = low_rew_tgt_num_sess;  
    data_behav.task_cond{counter_animal}    = task_cond_sess;        
    data_behav.choice{counter_animal}       = choice_sess;          
    data_behav.rew_cond{counter_animal}     = rew_cond_sess;          
    data_behav.jump_cond{counter_animal}    = jump_cond_sess; 
    data_behav.eye_vm_vis {counter_animal}  = eye_vm_sess_vis;
    data_behav.eye_vm_sac {counter_animal}  = eye_vm_sess_sac;
end

save_path = 'C:\Users\Jafar\Documents\reward\';
file_name = 'behave_data_session';
save(fullfile(save_path,'population_data',sprintf('%s.mat',file_name)), 'data_behav', '-v7.3');
end

