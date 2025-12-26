function JDM_build_CS_prim_sac_stim_type(data_path)
JDM_params_funcs
load(fullfile(data_path,'population_data','behave_data'));
user = 'JDM';
loc  = 'ctx';

path_data_monkey_sorted = params.path_data_monkey_sorted;
animal_list  = params.animal_list;
session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);

num_animal = numel(session_list);

% --------------------------------------------------------
% STORE DATA PER ANIMAL (one cell per animal)
% --------------------------------------------------------
cs_on_data  = cell(num_animal,1);
cs_cell_ids = cell(num_animal,1);
cs_rate     = cell(num_animal,1);
cs_vm       = cell(num_animal,1);
cs_rt       = cell(num_animal,1);

for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path      = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess          = numel(current_sess_list);
    sub_name          = animal_list{counter_animal};

    cs_on_data_sess  = cell(num_sess,1);
    cs_cell_ids_sess = cell(num_sess,1);
    cs_rate_sess     = cell(num_sess,1);
    cs_vm_sess       = cell(num_sess,1);
    cs_rt_sess       = cell(num_sess,1); 

    current_sess_cr_ma = data_behav.sess_cr_ma{counter_animal};
    current_sess_cr_ms = data_behav.sess_cr_ms{counter_animal};

    % --------------------------------------------------------
    % LOOP OVER SESSIONS FOR THIS ANIMAL
    % --------------------------------------------------------
    parfor counter_sess = 1:num_sess
%     for counter_sess = 1:num_sess
        current_sess = current_sess_list{counter_sess};
        disp(['     ' current_sess]);

        current_cr_ma = current_sess_cr_ma{counter_sess};
        current_cr_ms = current_sess_cr_ms{counter_sess};

        % load eye_traces
        data_eye   = MAF_load_eye_traces(current_path,current_sess);
        units_info = MAF_extract_cell_metadata(current_path,current_sess);
        cell_types = units_info.cell_type;
        cell_ids   = units_info.cell_list;

        cs_cell_ids_ = cell_ids(ismember(cell_types,{'CS','PC'}));

        [cs_on_data_,cs_rate_,cs_vm_,cs_rt_] = par_loop_cs_on_data( ...
            current_path,current_sess,data_eye,cs_cell_ids_, ...
            current_cr_ma,current_cr_ms,sub_name);

        cs_cell_ids_sess{counter_sess} = cs_cell_ids_;
        cs_on_data_sess{counter_sess}  = cs_on_data_;
        cs_rate_sess{counter_sess}     = cs_rate_;
        cs_vm_sess{counter_sess}       = cs_vm_;
        cs_rt_sess{counter_sess}       = cs_rt_; 
    end

    % --------------------------------------------------------
    % COLLAPSE SESSIONS -> ONE BLOCK PER ANIMAL
    % --------------------------------------------------------
    cs_on_data{counter_animal}  = vertcat(cs_on_data_sess{:});
    cs_cell_ids{counter_animal} = vertcat(cs_cell_ids_sess{:});

    % cat over sessions along "cell" dimension
    cs_rate{counter_animal} = cat(1,cs_rate_sess{:});   % [nCells_animal x 500 x 8 x 4 x 2]
    cs_vm{counter_animal}   = cat(1,cs_vm_sess{:});
    cs_rt{counter_animal}   = cat(1,cs_rt_sess{:});
end

% --------------------------------------------------------
% SAVE EVERYTHING IN ONE FILE, BUT PER-ANIMAL SEPARATED
% --------------------------------------------------------
file_name = 'cs_on_rate_prim_sac_stim_type';
save(fullfile(data_path,'population_data',[file_name '.mat']), ...
     'cs_on_data','cs_cell_ids','cs_rate','cs_vm','cs_rt','animal_list');
end

%% par_loop_cs_on_data
function [cs_on_data,cs_rate,cs_vm,cs_rt] = par_loop_cs_on_data( ...
    current_path,current_sess,data_eye,cell_ids, ...
    current_cr_ma,current_cr_ms,sub_name)


    num_cs     = numel(cell_ids);
    cs_on_data = cell(num_cs,1);

    % dims: [cell x time x dir x 4 stim-types x 2 epochs(vis/sac)]
    cs_rate = nan(num_cs,500,8,4,2);
    cs_vm   = nan(num_cs,500,8,4,2);
    cs_rt   = nan(num_cs,8,4); 
    sac_crit.ma   = current_cr_ma;
    sac_crit.ms   = current_cr_ms;

    % session date (for 132F splitting)
    cutoff = datetime("2024-08-20");
    d      = datetime(current_sess);  
    % Decide how to treat stimuli based on subject
    if strcmp(sub_name,'132F')
        % date-based split: period 1 = before cutoff, 2 = after
        if d < cutoff
            stim_mode = '132F_period1';  % High1/Low1
        else
            stim_mode = '132F_period2';  % High2/Low2
        end

        stim_number = [];   % not used for 132F

    elseif strcmp(sub_name,'65F')
        % tgt_num-based split: (0,1) for high & low
        stim_mode   = '65F_tgtnum';
        stim_number = [0; 1; 0; 1];  % [High1, High2, Low1, Low2]

    else
        % default: lump everything (you can tune this later if needed)
        stim_mode   = 'default';
        stim_number = [];
    end

    % ----------------- loop over CS cells -----------------
    parfor counter_cs = 1:num_cs
%     for counter_cs = 1:num_cs
        current_cs = cell_ids{counter_cs};

        data            = MAF_load_cell(current_path,current_sess,current_cs);
        data_recordings = data.data_recordings;

        cs_on_data{counter_cs} = MAF_CS_on_analysis_fr( ...
            [data_recordings.eye], ...
            [data_recordings.Neural_Data], ...
            [1,4,6,7],[40,85],[-70,30]);

        if cs_on_data{counter_cs}.vis.rho_avg > cs_on_data{counter_cs}.sac.rho_avg
            cs_on_ang = cs_on_data{counter_cs}.vis.ang_avg*pi/180;
        else
            cs_on_ang = cs_on_data{counter_cs}.sac.ang_avg*pi/180;
        end 

    res_rate = load_cs_data_tuned(current_path, ...
            current_sess,current_cs,1:8,sac_crit,data_eye, ...
            stim_mode,stim_number);

        % res_rate.rate_tot: [500 x 8 x 4 x 2]
        cs_rate(counter_cs,:,:,:,:) = res_rate.rate_tot;
        cs_vm(counter_cs,:,:,:,:)   = res_rate.vm_tot;
        cs_rt(counter_cs,:,:)       = res_rate.rt_tot; 
    end
end

%% load_cs_data_tuned
function res = load_cs_data_tuned(current_path,current_sess,current_cell, ...
                                  order_ang,sac_crit,data_eye, ...
                                  stim_mode,stim_number)
JDM_params_funcs;
tags = 1;

% output 4 conditions: [High1, High2, Low1, Low2]
n_timepoints = 500;
n_dirs       = 8;
n_conds      = 4;   % H1, H2, L1, L2

cs_rate_tot_sac = nan(n_timepoints, n_dirs, n_conds);
cs_rate_tot_vis = nan(n_timepoints, n_dirs, n_conds);

cs_vm_tot_sac   = nan(n_timepoints, n_dirs, n_conds);
cs_vm_tot_vis   = nan(n_timepoints, n_dirs, n_conds);
cs_rt_tot       = nan(n_dirs, n_conds);

data       = MAF_load_cell(current_path,current_sess,current_cell);
eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);

neural_sac_data = JDM_combine_dataset_sac([data.data_recordings], ...
                                          [data.rec_info], ...
                                          {'visual','onset'},250,tags);
eye_vm = MAF_combine_vm_traces(eye_traces,{'visual','onset'},250,tags);

% ========================================================
% Data binned by angle
% ========================================================
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
sac_data_dir_sac        = MAF_bin_by_ang(neur_data_sac,sac_data_sac,'onset');

for counter_dir = 1:n_dirs
    for cond_idx = 1:n_conds

        % =====================================================
        % -------- VISUAL -------------
        % =====================================================
        task_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.task_cond(:);
        task_tag  = sac_data_dir_vis(order_ang(counter_dir)).sac.tag(:);
        tgt_cond  = sac_data_dir_vis(order_ang(counter_dir)).sac.tgt_cond(:); % high/low
        tgt_num   = sac_data_dir_vis(order_ang(counter_dir)).sac.tgt_num(:);  % 0/1

        % Decide which trials belong to this condition
        istuned_vis = select_trials(stim_mode,cond_idx, ...
            task_cond,task_tag,tgt_cond,tgt_num, ...
            stim_number);
        istuned_vis = istuned_vis(:);  % enforce column

        has_cs_vis  = sac_data_dir_vis(order_ang(counter_dir)).sac.has_cs(:);
        idx_tuned_vis = istuned_vis & has_cs_vis;

        % if nothing matches, skip (stay NaN)
        if any(idx_tuned_vis)
            current_amp = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_amp(idx_tuned_vis);
            current_vel = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_vm_max(idx_tuned_vis);
            current_dec = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_offset(idx_tuned_vis) - ...
                           sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(idx_tuned_vis))*1e3;
            current_acc = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(idx_tuned_vis) - ...
                           sac_data_dir_vis(order_ang(counter_dir)).sac.time_onset(idx_tuned_vis))*1e3;
        
            % ---- Reaction time (ms) ----

            time_onset_sac  = sac_data_dir_vis(order_ang(counter_dir)).sac.time_onset(idx_tuned_vis);
            time_visual_sac = sac_data_dir_vis(order_ang(counter_dir)).sac.time_visual(idx_tuned_vis);
            current_rt      = (time_onset_sac - time_visual_sac) * 1e3;

            current_raster_vis = sac_data_dir_vis(order_ang(counter_dir)).CS(:,idx_tuned_vis);
            current_vm_vis     = sac_data_dir_vis(order_ang(counter_dir)).eye_vm(:,idx_tuned_vis);
        
            % remove invalid sacs (amp–vel & acc/dec filters)
            poly_x = sac_crit.ms(1,:);
            poly_y = sac_crit.ms(2,:);
            [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);
        
            acc_dec_ = current_dec./current_acc;
            poly_x   = sac_crit.ma(1,:);
            poly_y   = sac_crit.ma(2,:);
            [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);
        
            rm_ind = ~(in_ms | on_ms) | ~(in_acc | on_acc) | ...
                     (current_acc<=0) | (current_dec<=0);
        
            current_raster_vis(:,rm_ind) = [];
            current_vm_vis(:,rm_ind)     = [];
            current_rt(rm_ind)           = [];   
        
            if ~isempty(current_raster_vis)
                cs_rate_tot_vis(:, counter_dir, cond_idx) = ESN_smooth(mean(current_raster_vis,2))*1e3;
                cs_vm_tot_vis(:,   counter_dir, cond_idx) = mean(current_vm_vis,2,'omitnan');
                cs_rt_tot(counter_dir, cond_idx)          = mean(current_rt,'omitnan'); 
            end
        end


        % =====================================================
        % -------- SACCADE -------------
        % =====================================================
        task_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.task_cond(:);
        tgt_cond  = sac_data_dir_sac(order_ang(counter_dir)).sac.tgt_cond(:);
        task_tag  = sac_data_dir_sac(order_ang(counter_dir)).sac.tag(:);
        tgt_num   = sac_data_dir_sac(order_ang(counter_dir)).sac.tgt_num(:);

        istuned_sac = select_trials(stim_mode,cond_idx, ...
            task_cond,task_tag,tgt_cond,tgt_num, ...
            stim_number);
        istuned_sac = istuned_sac(:);

        has_cs_sac  = sac_data_dir_sac(order_ang(counter_dir)).sac.has_cs(:);
        idx_tuned_sac = istuned_sac & has_cs_sac;

        if any(idx_tuned_sac)
            current_amp = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_amp(idx_tuned_sac);
            current_vel = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_vm_max(idx_tuned_sac);
            current_dec = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_offset(idx_tuned_sac) - ...
                           sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(idx_tuned_sac))*1e3;
            current_acc = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(idx_tuned_sac) - ...
                           sac_data_dir_sac(order_ang(counter_dir)).sac.time_onset(idx_tuned_sac))*1e3;
            current_raster_sac = sac_data_dir_sac(order_ang(counter_dir)).CS(:,idx_tuned_sac);
            current_vm_sac     = sac_data_dir_sac(order_ang(counter_dir)).eye_vm(:,idx_tuned_sac);

            poly_x = sac_crit.ms(1,:);
            poly_y = sac_crit.ms(2,:);
            [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);

            acc_dec_ = current_dec./current_acc;
            poly_x   = sac_crit.ma(1,:);
            poly_y   = sac_crit.ma(2,:);
            [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);

            rm_ind = ~(in_ms | on_ms) | ~(in_acc | on_acc) | ...
                     (current_acc<=0) | (current_dec<=0);

            current_raster_sac(:,rm_ind) = [];
            current_vm_sac(:,rm_ind)     = [];

            cs_rate_tot_sac(:, counter_dir, cond_idx) = ESN_smooth(mean(current_raster_sac,2))*1e3;
            cs_vm_tot_sac(:,   counter_dir, cond_idx) = mean(current_vm_sac,2,'omitnan');
        end
    end
end

% [time x dir x cond x epoch(1=vis,2=sac)]
res.rate_tot = cat(4, cs_rate_tot_vis, cs_rate_tot_sac);
res.vm_tot   = cat(4, cs_vm_tot_vis,  cs_vm_tot_sac);
res.rt_tot   = cs_rt_tot;
end

%%
function istuned = select_trials(stim_mode,cond_idx, ...
                                 task_cond,task_tag,tgt_cond,tgt_num, ...
                                 stim_number)
% cond_idx: 1=High1, 2=High2, 3=Low1, 4=Low2
%
% task_cond: usually 1 for your main task
% task_tag : usually 1 for stim trials
% tgt_cond : 1 = high, 0 = low
% tgt_num  : stimulus index (0 or 1)

    switch stim_mode

        case '132F_period1'   % Before 2024-08-20: High1/Low1
            switch cond_idx
                case 1 % High1
                    istuned = (task_cond==1 & task_tag==1 & tgt_cond==1);
                case 2 % High2 -> empty
                    istuned = false(size(task_cond));
                case 3 % Low1
                    istuned = (task_cond==1 & task_tag==1 & tgt_cond==0);
                case 4 % Low2 -> empty
                    istuned = false(size(task_cond));
            end

        case '132F_period2'   % On/after 2024-08-20: High2/Low2
            switch cond_idx
                case 1 % High1 -> empty
                    istuned = false(size(task_cond));
                case 2 % High2
                    istuned = (task_cond==1 & task_tag==1 & tgt_cond==1);
                case 3 % Low1 -> empty
                    istuned = false(size(task_cond));
                case 4 % Low2
                    istuned = (task_cond==1 & task_tag==1 & tgt_cond==0);
            end

        case '65F_tgtnum'
            % stim_number = [0;1;0;1] for H1, H2, L1, L2
            this_num = stim_number(cond_idx);

            switch cond_idx
                case {1,2}   % High1, High2
                    istuned = (task_cond==1 & task_tag==1 & ...
                               tgt_cond==1 & tgt_num==this_num);
                case {3,4}   % Low1, Low2
                    istuned = (task_cond==1 & task_tag==1 & ...
                               tgt_cond==0 & tgt_num==this_num);
            end

        otherwise
            % fallback: collapse all high/low, no split
            switch cond_idx
                case {1,2}   % both high
                    istuned = (task_cond==1 & task_tag==1 & tgt_cond==1);
                case {3,4}   % both low
                    istuned = (task_cond==1 & task_tag==1 & tgt_cond==0);
            end
    end
end




