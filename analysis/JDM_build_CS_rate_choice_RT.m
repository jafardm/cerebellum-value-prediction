function JDM_build_CS_rate_choice_RT(data_path)
% =============================================================
% Build CS data for choice trials binned by reaction time (RT)
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
perf_H1_animal = cell(0);
perf_H2_animal = cell(0);

% === Fixed RT edges (applied to every session) ===
rt_edges = [60:20:200 600];     % ms

for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess = numel(current_sess_list);
    sub_name          = animal_list{counter_animal};

    current_sess_cr_ma = data_behav.sess_cr_ma{counter_animal};
    current_sess_cr_ms = data_behav.sess_cr_ms{counter_animal};

    % --- holders per session ---
    cs_on_data_sess  = cell(num_sess,1);
    cs_cell_ids_sess = cell(num_sess,1);
    cs_rate_sess     = cell(num_sess,1);
    cs_vm_sess       = cell(num_sess,1);
    choice_perf_RT_sess  = cell(num_sess,1); 
    choice_RT_bin_sess = cell(num_sess,1); 

    perf_H1_sess = cell(num_sess,1);
    perf_H2_sess = cell(num_sess,1);

    
    %parfor
   for counter_sess = 1:num_sess
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

        % =====================================================
        % === Compute choice performance for RT bins (per session)
        % =====================================================

        cell_num_rec = sum(cell2mat(cellfun(@isempty, ...
                                units_info.cell_ids,'UniformOutput',false)),2);
    
        [~,selected_cell_id] = min(cell_num_rec);
        data = MAF_load_cell(current_path,current_sess,...
                             units_info.cell_list{selected_cell_id});
        
        rec_flag   = data.rec_info.rec_flag;
        eye_traces = data_eye.eye_traces(rec_flag);
        
        sacs_tr    = [eye_traces.sac];
        tags_tr    = [sacs_tr.tag]';
        
        data_rec   = data.data_recordings;
        eye        = [data_rec.eye];
        sac        = [eye.sac];
        
        sac_tag_   = [sac.tag]';
        sac_amp_   = [sac.eye_amp]';
        
        assert(all(tags_tr == sac_tag_));

        sac_choice = [sac.choice]';     % 1 = High, 0 = Low
        task_cond_ = [sac.task_cond]'; 
        sac_rt_    = ([sac.time_onset] - [sac.time_visual])' * 1e3; % ms
       
        stim_high_num_ = [sac.high_rew_tgt_num]'; 
        
        % === Apply your behavioral filters ===
        idx_rxt      = sac_rt_ > 0 & sac_rt_ < 600;   % valid RT range
        tag_1        = sac_tag_ == 1;                 % only tag==1
        sac_amp_idx  = sac_amp_ > 3.0;                % amplitude threshold

        
        % Outlier removal
        res          = MAF_find_sac_outliers(sac);
        not_outlier  = ~res.rm_ind';
 
        % --- choice task index (task_cond_ == 0 means "choice trials") ---
        idx_choice = (task_cond_ == 0);
        
        % Combined good-trial index: only CHOICE trials
        idx_good = idx_rxt & not_outlier & tag_1 & sac_amp_idx & idx_choice;
        
        % Extract cleaned behavioral signals (only choice trials)
        sess_rt      = sac_rt_(idx_good);         % RT of choice trials
        sess_choice  = sac_choice(idx_good);      % 1 = High, 0 = Low
        sess_task    = task_cond_(idx_good);      % should all be 0 now

        % =====================================================
        %   PERFORMANCE FOR H1 and H2 TRIAL TYPES
        %   (two high stimuli: 0 = H1, 1 = H2)
        %   sac_choice: 1 = High, 0 = Low
        % =====================================================
 
        
        % high-value target identity for these good trials
        sess_high_id = stim_high_num_(idx_good);  
        cutoff = datetime("2024-08-20");
        d = datetime(current_sess);
        
        % initialize logical arrays
        idx_H1 = false(size(sess_high_id));
        idx_H2 = false(size(sess_high_id));
        
        % 132F
        if strcmp(sub_name,'132F')
        
            if d < cutoff
                % Before cutoff: all high trials → H1
                idx_H1 = (sess_high_id == 0) | (sess_high_id == 1);
                idx_H2 = false;
        
            else
                % After cutoff: all high trials → H2
                idx_H2 = (sess_high_id == 0) | (sess_high_id == 1);
                idx_H1 = false;
            end
        
       % 65F
        elseif strcmp(sub_name,'65F')
            % For 65F, ID=0 is H1, ID=1 is H2
            idx_H1 = sess_high_id == 0;
            idx_H2 = sess_high_id == 1;
        end
        
        % Compute performance
       
        if any(idx_H1)
            perf_H1 = sum(sess_choice(idx_H1) == 1) / sum(idx_H1);
        else
            perf_H1 = nan;
        end
        
        if any(idx_H2)
            perf_H2 = sum(sess_choice(idx_H2) == 1) / sum(idx_H2);
        else
            perf_H2 = nan;
        end

        perf_H1_sess{counter_sess} = perf_H1;
        perf_H2_sess{counter_sess} = perf_H2; 
 
      % Sanity check (optional): all(sess_task == 0)
        
        % === RT BINNING ===
        rt_bin    = discretize(sess_rt, rt_edges);
        n_rt_bins = numel(rt_edges) - 1;
        
        perf_RT         = zeros(n_rt_bins,1);  % performance per RT bin
        n_choice_RT_bin = nan(n_rt_bins,1);    % #choice trials per RT bin
        RT_bins         = cell(n_rt_bins,1);
        
        for b = 1:n_rt_bins
            idx_bin = (rt_bin == b);         % trials whose RT falls in bin b
            % these are *already* choice trials, because of idx_choice in idx_good
            n_choice_RT_bin(b) = sum(idx_bin);
            RT_bins{b} = sess_rt(idx_bin);
            if n_choice_RT_bin(b) > 0
                % performance = (# High choices in bin) / (# choice trials in bin)
                perf_RT(b) = sum(sess_choice(idx_bin) == 1) / n_choice_RT_bin(b);
            end
        end
      
        choice_perf_RT_sess{counter_sess} = perf_RT;
        choice_RT_bin_sess{counter_sess} = RT_bins;

        perf_H1_sess{counter_sess} = perf_H1;
        perf_H2_sess{counter_sess} = perf_H2;
    end
    
    perf_H1_animal{counter_animal} = perf_H1_sess;
    perf_H2_animal{counter_animal} = perf_H2_sess;


    % --- concatenate across sessions for this animal ---
    cs_on_data  = [cs_on_data; vertcat(cs_on_data_sess{:})];
    cs_cell_ids = [cs_cell_ids; vertcat(cs_cell_ids_sess{:})];
    choice_perf_RT{counter_animal} = choice_perf_RT_sess;
    choice_RT_bins{counter_animal} = choice_RT_bin_sess;
    
    if counter_animal == 1
        cs_rate = cat(1,cs_rate_sess{:});
        cs_vm   = cat(1,cs_vm_sess{:});
    else
        cs_rate = cat(1,cs_rate, cat(1,cs_rate_sess{:}));
        cs_vm   = cat(1,cs_vm, cat(1,cs_vm_sess{:}));
    end
 
end

% --- save ---
file_name = 'cs_on_rate_choice_sac_RT';
save(fullfile(data_path,'population_data',sprintf('%s.mat',file_name)), ...
     'cs_on_data','cs_cell_ids','cs_rate','cs_vm','animal_list',...
     'perf_H1_animal','perf_H2_animal','choice_perf_RT','choice_RT_bins',...
     'session_list','rt_edges','-v7.3');

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
%parfor
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
choice_conds  =  ...
 [0, 1, 1;  % task, tag, High
  0, 1, 0];  % task, tag, Low

n_timepoints = 500;
n_dirs       = 8;
n_cond  = size(choice_conds,1);
n_conds      = n_cond ;
n_rt_bins    = numel(rt_edges) - 1;

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

 for counter_dir = 1:n_dirs
    for cond_idx = 1:n_conds

        % ===================== VISUAL =====================
        task_cond   = sac_data_dir_vis(order_ang(counter_dir)).sac.task_cond;
        task_tag    = sac_data_dir_vis(order_ang(counter_dir)).sac.tag;
        ch_cond    = sac_data_dir_vis(order_ang(counter_dir)).sac.choice;
        has_cs_vis  = sac_data_dir_vis(order_ang(counter_dir)).sac.has_cs;

        istuned_vis = ismember([task_cond' task_tag' ch_cond'], ...
                               choice_conds(cond_idx,:), 'rows');

        time_onset  = sac_data_dir_vis(order_ang(counter_dir)).sac.time_onset;
        time_visual = sac_data_dir_vis(order_ang(counter_dir)).sac.time_visual;
        current_rt  = (time_onset - time_visual) * 1e3; % ms
        rt_bin_vis  = discretize(current_rt, rt_edges);

        for iRT = 1:n_rt_bins
            trials_in_bin = (rt_bin_vis == iRT);

            % ---- VIS ----
            idx_tuned_vis = istuned_vis(:) & has_cs_vis(:) & trials_in_bin(:);

            if any(idx_tuned_vis)
                current_amp      = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_amp(idx_tuned_vis);
                current_vel      = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_vm_max(idx_tuned_vis);
                current_dec      = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_offset(idx_tuned_vis) - ...
                                    sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(idx_tuned_vis)) * 1e3;
                current_acc      = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(idx_tuned_vis) - ...
                                    sac_data_dir_vis(order_ang(counter_dir)).sac.time_onset(idx_tuned_vis)) * 1e3;
                current_raster_vis = sac_data_dir_vis(order_ang(counter_dir)).CS(:,idx_tuned_vis);
                current_vm_vis     = sac_data_dir_vis(order_ang(counter_dir)).eye_vm(:,idx_tuned_vis);

                % main sequence + acc/dec filters
                poly_x = sac_crit.ms(1,:);
                poly_y = sac_crit.ms(2,:);
                [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);

                acc_dec_ = current_dec ./ current_acc;
                poly_x   = sac_crit.ma(1,:);
                poly_y   = sac_crit.ma(2,:);
                [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);

                rm_ind = ~(in_ms | on_ms) | ~(in_acc | on_acc) | (current_acc<=0) | (current_dec<=0);

                current_raster_vis(:,rm_ind) = [];
                current_vm_vis(:,rm_ind)     = [];

                cs_rate_tot_vis(:, counter_dir, cond_idx, iRT) = ESN_smooth(mean(current_raster_vis,2,'omitnan')) * 1e3;
                cs_vm_tot_vis(:,   counter_dir, cond_idx, iRT) = mean(current_vm_vis,2);
            else
                cs_rate_tot_vis(:, counter_dir, cond_idx, iRT) = nan;
                cs_vm_tot_vis(:,   counter_dir, cond_idx, iRT) = nan;
            end

            % ===================== SACCADE =====================
            task_cond   = sac_data_dir_sac(order_ang(counter_dir)).sac.task_cond;
            task_tag    = sac_data_dir_sac(order_ang(counter_dir)).sac.tag;
            ch_cond    = sac_data_dir_sac(order_ang(counter_dir)).sac.choice;
            has_cs_sac  = sac_data_dir_sac(order_ang(counter_dir)).sac.has_cs;

            istuned_sac = ismember([task_cond' task_tag' ch_cond'], ...
                                   choice_conds(cond_idx,:), 'rows');  

            time_onset_sac  = sac_data_dir_sac(order_ang(counter_dir)).sac.time_onset;
            time_visual_sac = sac_data_dir_sac(order_ang(counter_dir)).sac.time_visual;
            current_rt_sac  = (time_onset_sac - time_visual_sac) * 1e3;
            rt_bin_sac      = discretize(current_rt_sac, rt_edges);
            trials_in_bin_sac = (rt_bin_sac == iRT);

            idx_tuned_sac = istuned_sac(:) & has_cs_sac(:) & trials_in_bin_sac(:);

            if any(idx_tuned_sac)
                current_amp = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_amp(idx_tuned_sac);
                current_vel = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_vm_max(idx_tuned_sac);
                current_dec = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_offset(idx_tuned_sac) - ...
                               sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(idx_tuned_sac)) * 1e3;
                current_acc = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(idx_tuned_sac) - ...
                               sac_data_dir_sac(order_ang(counter_dir)).sac.time_onset(idx_tuned_sac)) * 1e3;
                current_raster_sac = sac_data_dir_sac(order_ang(counter_dir)).CS(:,idx_tuned_sac);
                current_vm_sac     = sac_data_dir_sac(order_ang(counter_dir)).eye_vm(:,idx_tuned_sac);

                poly_x = sac_crit.ms(1,:);
                poly_y = sac_crit.ms(2,:);
                [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);

                acc_dec_ = current_dec ./ current_acc;
                poly_x   = sac_crit.ma(1,:);
                poly_y   = sac_crit.ma(2,:);
                [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);

                rm_ind = ~(in_ms | on_ms) | ~(in_acc | on_acc) | (current_acc<=0) | (current_dec<=0);

                current_raster_sac(:,rm_ind) = [];
                current_vm_sac(:,rm_ind)     = [];

                cs_rate_tot_sac(:, counter_dir, cond_idx, iRT) = ESN_smooth(mean(current_raster_sac,2,'omitnan')) * 1e3;
                cs_vm_tot_sac(:,   counter_dir, cond_idx, iRT) = mean(current_vm_sac,2);
            else
                cs_rate_tot_sac(:, counter_dir, cond_idx, iRT) = nan;
                cs_vm_tot_sac(:,   counter_dir, cond_idx, iRT) = nan;
            end

        end % iRT
    end % cond
end % dir
res.rate_tot = cat(5, cs_rate_tot_vis, cs_rate_tot_sac); % [500 x 8 x n_conds x nRTbins x 2]
res.vm_tot   = cat(5, cs_vm_tot_vis, cs_vm_tot_sac);

end% function end
