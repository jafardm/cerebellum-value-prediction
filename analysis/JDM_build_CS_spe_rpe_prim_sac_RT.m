function JDM_build_CS_spe_rpe_prim_sac_RT(data_path)
% =============================================================
% Build CS data binned by reaction time (RT)
% =============================================================
JDM_params_funcs
data_behav = load(fullfile(data_path,'population_data','behave_data'));
data_behav = data_behav.data_behav;

user = 'JDM'; 
loc  = 'ctx';
path_data_monkey_sorted = params.path_data_monkey_sorted;
animal_list = params.animal_list;
session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);
num_animal = numel(session_list);

cs_on_data  = cell(0);
cs_cell_ids = cell(0);
cs_rate     = [];
cs_vm       = [];

% === Fixed RT edges (applied to every session) ===
rt_edges = [60:20:200 600];     % ms


for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess = numel(current_sess_list);

    current_sess_cr_ma = data_behav.sess_cr_ma{counter_animal};
    current_sess_cr_ms = data_behav.sess_cr_ms{counter_animal};

    % --- holders per session ---
    cs_on_data_sess  = cell(num_sess,1);
    cs_cell_ids_sess = cell(num_sess,1);
    cs_rate_sess     = cell(num_sess,1);
    cs_vm_sess       = cell(num_sess,1);

    parfor counter_sess = 1:num_sess
        current_sess = current_sess_list{counter_sess};
        disp(['     ' current_sess]);

        current_cr_ma = current_sess_cr_ma{counter_sess};
        current_cr_ms = current_sess_cr_ms{counter_sess};

        % --- load data ---
        data_eye = MAF_load_eye_traces(current_path,current_sess);
        units_info = MAF_extract_cell_metadata(current_path,current_sess);
        cell_types = units_info.cell_type;
        cell_ids   = units_info.cell_list;
        cs_cell_ids_ = cell_ids(ismember(cell_types,{'CS','PC'}));

        % --- compute ---
        [cs_on_data_, cs_rate_, cs_vm_] = par_loop_cs_on_data_RT( ...
            current_path, current_sess, data_eye, cs_cell_ids_, ...
            current_cr_ma, current_cr_ms, rt_edges);

        cs_cell_ids_sess{counter_sess} = cs_cell_ids_;
        cs_on_data_sess{counter_sess}  = cs_on_data_;
        cs_rate_sess{counter_sess}     = cs_rate_;
        cs_vm_sess{counter_sess}       = cs_vm_;
    end

    % --- concatenate across sessions for this animal ---
    cs_on_data  = [cs_on_data; vertcat(cs_on_data_sess{:})];
    cs_cell_ids = [cs_cell_ids; vertcat(cs_cell_ids_sess{:})];

    if counter_animal == 1
        cs_rate = cat(1,cs_rate_sess{:});
        cs_vm   = cat(1,cs_vm_sess{:});
    else
        cs_rate = cat(1,cs_rate, cat(1,cs_rate_sess{:}));
        cs_vm   = cat(1,cs_vm, cat(1,cs_vm_sess{:}));
    end
end

% --- save ---
file_name = 'cs_on_rate_spe_rpe_prim_sac_RT';
save(fullfile(data_path,'population_data',sprintf('%s.mat',file_name)), ...
     'cs_on_data','cs_cell_ids','cs_rate','cs_vm','rt_edges','-v7.3');

end

%% =====================================================================
function [cs_on_data, cs_rate, cs_vm] = par_loop_cs_on_data_RT( ...
    current_path, current_sess, data_eye, cell_ids, ...
    current_cr_ma, current_cr_ms, rt_edges)

num_cs = numel(cell_ids);
cs_on_data = cell(num_cs,1);

sac_crit.ma = current_cr_ma;
sac_crit.ms = current_cr_ms;

n_rt_bins = numel(rt_edges) - 1;
cs_rate = nan(num_cs,500,8,2,n_rt_bins,2); % [cell × time × dir × reward × RTbin × epoch]
cs_vm   = nan(num_cs,500,8,2,n_rt_bins,2);

parfor counter_cs = 1:num_cs
    current_cs = cell_ids{counter_cs};
    data = MAF_load_cell(current_path,current_sess,current_cs);
    data_recordings = data.data_recordings;

    cs_on_data{counter_cs} = MAF_CS_on_analysis_fr( ...
        [data_recordings.eye],[data_recordings.Neural_Data], ...
        [1,4,6,7],[40,85],[-70,30]);

    res_rate = load_cs_data_tuned_RT(current_path,current_sess, ...
        current_cs,1:8,sac_crit,data_eye,rt_edges);

    % res_rate.rate_tot is [500 × 8 × ncond × n_rt_bins × 2]
    cs_rate(counter_cs,:,:,:,:,:) = res_rate.rate_tot;
    cs_vm(counter_cs,:,:,:,:,:)   = res_rate.vm_tot;
end
end

%% =====================================================================
function res = load_cs_data_tuned_RT(current_path,current_sess,current_cell,order_ang,sac_crit,data_eye,rt_edges)
JDM_params_funcs;
tags = 1;
force_conds  =  ...
 [1, 1, 1;  % High
  1, 1, 0];  % Low


n_timepoints = 500;
n_dirs       = 8;
n_forc_cond  = size(force_conds,1);
n_conds      = n_forc_cond ;
n_rt_bins = numel(rt_edges) - 1;

cs_rate_tot_sac    = nan(n_timepoints, n_dirs, n_conds,n_rt_bins);
cs_rate_tot_vis    = nan(n_timepoints, n_dirs, n_conds,n_rt_bins);

cs_vm_tot_sac    = nan(n_timepoints, n_dirs, n_conds,n_rt_bins);
cs_vm_tot_vis    = nan(n_timepoints, n_dirs, n_conds,n_rt_bins);

data       = MAF_load_cell(current_path,current_sess,current_cell);
eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);

%neural_sac_data2 = MAF_combine_dataset_sac([data.data_recordings],[data.rec_info],{'visual','onset'},250,tags);
neural_sac_data = JDM_combine_dataset_sac([data.data_recordings],[data.rec_info],{'visual','onset'},250,tags);
eye_vm = MAF_combine_vm_traces(eye_traces,{'visual','onset'},250,tags);

% vis
neur_data_vis.eye_vm    = eye_vm.visual.eye_vm;
neur_data_vis.CS        = neural_sac_data.visual.CS;
sac_data_vis.sac        = neural_sac_data.visual.sac;
sac_data_vis.fix        = neural_sac_data.visual.fix;
sac_data_vis.sac.has_ss = neural_sac_data.visual.has_ss;
sac_data_vis.sac.has_cs = neural_sac_data.visual.has_cs;
sac_data_dir_vis        = MAF_bin_by_ang(neur_data_vis,sac_data_vis,'visual');

% sac
neur_data_sac.eye_vm    = eye_vm.onset.eye_vm;
neur_data_sac.CS        = neural_sac_data.onset.CS;
sac_data_sac.sac        = neural_sac_data.onset.sac;
sac_data_sac.fix        = neural_sac_data.onset.fix;
sac_data_sac.sac.has_cs = neural_sac_data.onset.has_cs;
sac_data_dir_sac   = MAF_bin_by_ang(neur_data_sac,sac_data_sac,'onset');

%parfor counter_dir = 1:8
 for counter_dir = 1:n_dirs  
    for cond_idx = 1:n_conds
       task_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.task_cond;
       task_tag = sac_data_dir_vis(order_ang(counter_dir)).sac.tag;
       tgt_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.tgt_cond;
       has_cs_vis = sac_data_dir_vis(order_ang(counter_dir)).sac.has_cs;
       istuned_vis = ismember([task_cond' task_tag' tgt_cond'],force_conds(cond_idx,:),'rows');
       % compute RT
       time_onset = sac_data_dir_vis(order_ang(counter_dir)).sac.time_onset;
       time_visual = sac_data_dir_vis(order_ang(counter_dir)).sac.time_visual; 
       current_rt = (time_onset - time_visual)*1e3; % ms
       rt_bin_vis = discretize(current_rt, rt_edges);

       for iRT = 1:n_rt_bins
        trials_in_bin = (rt_bin_vis == iRT);
        % vis
        idx_tuned_vis = istuned_vis(:) & has_cs_vis(:) & trials_in_bin(:);

        if ~any(idx_tuned_vis)
            cs_rate_tot_vis(:, counter_dir, cond_idx, iRT) = nan;
            cs_vm_tot_vis(:, counter_dir, cond_idx, iRT)   = nan;
            continue
        end

        current_amp = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_amp(idx_tuned_vis);
        current_vel = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_vm_max(idx_tuned_vis);
        current_dec = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_offset(idx_tuned_vis) -...
            sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(idx_tuned_vis))*1e3;
        current_acc = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(idx_tuned_vis) -...
            sac_data_dir_vis(order_ang(counter_dir)).sac.time_onset(idx_tuned_vis))*1e3;
        current_raster_vis = sac_data_dir_vis(order_ang(counter_dir)).CS(:,idx_tuned_vis);
        current_vm_vis = sac_data_dir_vis(order_ang(counter_dir)).eye_vm(:,idx_tuned_vis);

        % remove trials without spike if ss also remove invalid sacs
        poly_x = sac_crit.ms(1,:);
        poly_y = sac_crit.ms(2,:);
        [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);
    
        acc_dec_ = current_dec./current_acc;
        poly_x = sac_crit.ma(1,:);
        poly_y = sac_crit.ma(2,:);
        [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);
        
        % negative acc, dec, outliers in main-seq and amp-acc/dec plot
        rm_ind = ~(in_ms | on_ms) | ~(in_acc | on_acc) | (current_acc<=0) | (current_dec<=0) ;
    
        current_raster_vis(:,rm_ind) = [];
        current_vm_vis(:,rm_ind) = [];

        cs_rate_tot_vis(:, counter_dir, cond_idx,iRT) = ESN_smooth(mean(current_raster_vis,2))*1e3;
        cs_vm_tot_vis(:, counter_dir, cond_idx,iRT) = mean(current_vm_vis,2);
    
       % sac
       task_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.task_cond;
       tgt_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.tgt_cond;
       task_tag = sac_data_dir_sac(order_ang(counter_dir)).sac.tag;
       has_cs_sac = sac_data_dir_sac(order_ang(counter_dir)).sac.has_cs;
       istuned_sac = ismember([task_cond' task_tag' tgt_cond'],force_conds(cond_idx,:),'rows');

       time_onset_sac  = sac_data_dir_sac(order_ang(counter_dir)).sac.time_onset;
       time_visual_sac = sac_data_dir_sac(order_ang(counter_dir)).sac.time_visual;
       current_rt_sac  = (time_onset_sac - time_visual_sac) * 1e3;
       rt_bin_sac      = discretize(current_rt_sac, rt_edges);
       trials_in_bin_sac = (rt_bin_sac == iRT);
       idx_tuned_sac = istuned_sac(:) & has_cs_sac(:) & trials_in_bin_sac(:);

       if ~any(idx_tuned_sac)
        cs_rate_tot_sac(:, counter_dir, cond_idx, iRT) = nan;
        cs_vm_tot_sac(:, counter_dir, cond_idx, iRT)   = nan;
        continue
       end

        current_amp = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_amp(idx_tuned_sac);
        current_vel = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_vm_max(idx_tuned_sac);
        current_dec = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_offset(idx_tuned_sac) -...
            sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(idx_tuned_sac))*1e3;
        current_acc = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(idx_tuned_sac) -...
            sac_data_dir_sac(order_ang(counter_dir)).sac.time_onset(idx_tuned_sac))*1e3;
        current_raster_sac = sac_data_dir_sac(order_ang(counter_dir)).CS(:,idx_tuned_sac);
        current_vm_sac = sac_data_dir_sac(order_ang(counter_dir)).eye_vm(:,idx_tuned_sac);
    
        % remove trials without spike if ss also remove invalid sacs
        poly_x = sac_crit.ms(1,:);
        poly_y = sac_crit.ms(2,:);
        [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);
    
        acc_dec_ = current_dec./current_acc;
        poly_x = sac_crit.ma(1,:);
        poly_y = sac_crit.ma(2,:);
        [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);
        
        % no spike, negative acc, dec, outliers in main-seq and amp-acc/dec plot
        rm_ind = ~(in_ms | on_ms) | ~(in_acc | on_acc) | (current_acc<=0) | (current_dec<=0);
    
        current_raster_sac(:,rm_ind) = [];
        current_vm_sac(:,rm_ind) = [];
  
        cs_rate_tot_sac(:, counter_dir, cond_idx,iRT) = ESN_smooth(mean(current_raster_sac,2))*1e3;
        cs_vm_tot_sac(:, counter_dir, cond_idx,iRT)   = mean(current_vm_sac,2);
      end
    end % cond end
 end % dir end
res.rate_tot = cat(5, cs_rate_tot_vis, cs_rate_tot_sac); % [500 x 8 x n_conds x nRTbins x 2]
res.vm_tot   = cat(5, cs_vm_tot_vis, cs_vm_tot_sac);

end% function end
