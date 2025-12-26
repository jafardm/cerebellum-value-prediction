function JDM_build_cs_rate_velocity_bins(data_path)
JDM_params_funcs
load(fullfile(data_path,'population_data','behave_data'));
max_vlocities = load(fullfile(data_path, 'population_data', 'sac_kinematic_prim_sac.mat'),'tag1_max_vel');

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

for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    % ---- separate velocities into three tertiles (edges) ----
    max_vel_sub = max_vlocities.tag1_max_vel{counter_animal,1};  % cell of per-session vectors
    all_vel = vertcat(max_vel_sub{:});                           % concatenate
    
    % Get tertile thresholds; guard against pathological ties
    q = quantile(all_vel, [1/3, 2/3]);
    if ~issorted(q) || any(~isfinite(q))
        % Fallback if something weird happens
        q = quantile(all_vel, [0.3334, 0.6667]);
    end
    % Ensure strictly non-decreasing with tiny eps nudges
    if q(2) <= q(1)
        q(2) = q(1) + eps(q(1));
    end
    
    vel_edges = [0, q(1), q(2), 10000];  
    current_path = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess = numel(current_sess_list);
    cs_on_data_sess  = cell(num_sess,1);
    cs_cell_ids_sess = cell(num_sess,1);
    cs_rate_sess     = cell(num_sess,1);
    cs_vm_sess       = cell(num_sess,1);

    current_sess_cr_ma = data_behav.sess_cr_ma{counter_animal};
    current_sess_cr_ms = data_behav.sess_cr_ms{counter_animal};
    %parfor
    parfor counter_sess = 1:num_sess
        current_sess = current_sess_list{counter_sess};
        disp(['     ' current_sess]);
        current_cr_ma = current_sess_cr_ma{counter_sess};
        current_cr_ms = current_sess_cr_ms{counter_sess};
        % load eye_traces
        data_eye = MAF_load_eye_traces(current_path,current_sess);
        units_info = MAF_extract_cell_metadata(current_path,current_sess);
        cell_types = units_info.cell_type;
        cell_ids   = units_info.cell_list;
        cs_cell_ids_ = cell_ids(ismember(cell_types,{'CS','PC'}));

        [cs_on_data_,~,cs_rate_,cs_vm_] = par_loop_cs_on_data...
            (current_path,current_sess,data_eye,cs_cell_ids_,vel_edges,...
            current_cr_ma,current_cr_ms);
        cs_cell_ids_sess{counter_sess} = cs_cell_ids_;
        cs_on_data_sess{counter_sess}  = cs_on_data_;
        cs_rate_sess{counter_sess}     = cs_rate_;
        cs_vm_sess{counter_sess}       = cs_vm_;
    end
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
file_name = 'cs_on_rate_velocityBins_prim_sac';
save(fullfile(data_path,'population_data',sprintf('%s.mat',file_name)),...
   'cs_on_data','cs_cell_ids','cs_rate','cs_vm');
end

%% par_loop_cs_on_data
function [cs_on_data,perturb,cs_rate,cs_vm] = par_loop_cs_on_data...
    (current_path,current_sess,data_eye,cell_ids,vel_edges,...
    current_cr_ma,current_cr_ms)
    num_cs = numel(cell_ids);
    cs_on_data = cell(num_cs,1);
    perturb = cell(num_cs,1);

    n_vel_bins = numel(vel_edges) - 1;
    n_events   = 3;  % visual/onset/offset
    n_conds    = 2;  % from force_conds above
    n_dirs     = 8; 
    n_time     = 500;
    
    cs_rate = nan(num_cs, n_time, n_dirs, n_conds, n_vel_bins, n_events);
    cs_vm   = nan(num_cs, n_time, n_dirs, n_conds, n_vel_bins, n_events);

    sac_crit.ma   = current_cr_ma;
    sac_crit.ms   = current_cr_ms;
   %parfor
   parfor counter_cs = 1:num_cs
        current_cs = cell_ids{counter_cs};
    
        data = MAF_load_cell(current_path,current_sess,current_cs);
        data_recordings = data.data_recordings;
        cs_on_data{counter_cs} = MAF_CS_on_analysis_fr([data_recordings.eye],...
            [data_recordings.Neural_Data],[1,4,6,7],[40,85],[-70,30]);

        if cs_on_data{counter_cs}.vis.rho_avg > cs_on_data{counter_cs}.sac.rho_avg
            cs_on_ang = cs_on_data{counter_cs}.vis.ang_avg*pi/180;
        else
            cs_on_ang = cs_on_data{counter_cs}.sac.ang_avg*pi/180;
        end

        res_rate = load_cs_data_tuned(current_path,...
            current_sess,current_cs,1:8,sac_crit,data_eye,vel_edges);
  
        cs_rate(counter_cs,:,:,:,:,:) = res_rate.rate_tot;
        cs_vm(counter_cs,:,:,:,:,:)   = res_rate.vm_tot;
    end
end
%% load_cs_data_tuned

function res = load_cs_data_tuned(current_path,current_sess,current_cell,order_ang,sac_crit,data_eye,vel_bins)
% Load CS rate / eye velocity binned by direction, condition, and velocity.
% Returns:
%   res.rate_tot   [n_time x n_dirs x n_conds x n_vel_bins x n_events]
%   res.vm_tot     same shape
%   res.vel_bins   vector of bin edges used
%


JDM_params_funcs;
tags = 1;

% ------------ Force conditions ------------
force_conds = [1, 1, 1;   % H
               1, 1, 0];  % L

% ------------ Velocity bins ------------
n_vel_bins = numel(vel_bins) - 1;

% ------------ Dimensions / prealloc ------------
n_timepoints = 500;
n_dirs       = 8;
n_conds      = size(force_conds,1);

cs_rate_tot_vis    = nan(n_timepoints, n_dirs, n_conds, n_vel_bins);
cs_rate_tot_sac    = nan(n_timepoints, n_dirs, n_conds, n_vel_bins);
cs_rate_tot_offset = nan(n_timepoints, n_dirs, n_conds, n_vel_bins);

cs_vm_tot_vis    = nan(n_timepoints, n_dirs, n_conds, n_vel_bins);
cs_vm_tot_sac    = nan(n_timepoints, n_dirs, n_conds, n_vel_bins);
cs_vm_tot_offset = nan(n_timepoints, n_dirs, n_conds, n_vel_bins);

% ------------ Load cell data ------------
data       = MAF_load_cell(current_path,current_sess,current_cell);
eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);

neural_sac_data = JDM_combine_dataset_sac([data.data_recordings],[data.rec_info],{'visual','onset','offset'},250,tags);
eye_vm          = MAF_combine_vm_traces(eye_traces,{'visual','onset','offset'},250,tags);

% VISUAL
neur_data_vis.eye_vm    = eye_vm.visual.eye_vm;
neur_data_vis.CS        = neural_sac_data.visual.CS;
sac_data_vis.sac        = neural_sac_data.visual.sac;
sac_data_vis.fix        = neural_sac_data.visual.fix;
sac_data_vis.sac.has_ss = neural_sac_data.visual.has_ss;
sac_data_vis.sac.has_cs = neural_sac_data.visual.has_cs;
sac_data_dir_vis        = MAF_bin_by_ang(neur_data_vis,sac_data_vis,'visual');

% ONSET
neur_data_sac.eye_vm    = eye_vm.onset.eye_vm;
neur_data_sac.CS        = neural_sac_data.onset.CS;
sac_data_sac.sac        = neural_sac_data.onset.sac;
sac_data_sac.fix        = neural_sac_data.onset.fix;
sac_data_sac.sac.has_cs = neural_sac_data.onset.has_cs;
sac_data_dir_sac        = MAF_bin_by_ang(neur_data_sac,sac_data_sac,'onset');

% OFFSET
neur_data_offset.eye_vm    = eye_vm.offset.eye_vm;
neur_data_offset.CS        = neural_sac_data.offset.CS;
sac_data_offset.sac        = neural_sac_data.offset.sac;
sac_data_offset.fix        = neural_sac_data.offset.fix;
sac_data_offset.sac.has_cs = neural_sac_data.offset.has_cs;
sac_data_dir_offset        = MAF_bin_by_ang(neur_data_offset,sac_data_offset,'offset');

% ------------ Main loop ------------
for counter_dir = 1:n_dirs
    for cond_idx = 1:n_conds
        % Visual
        [cs_rate_tot_vis(:, counter_dir, cond_idx, :), ...
         cs_vm_tot_vis(:,   counter_dir, cond_idx, :)] = process_event_for_bins( ...
             sac_data_dir_vis, order_ang(counter_dir), cond_idx, ...
             force_conds, sac_crit, vel_bins, n_timepoints);

        % Onset
        [cs_rate_tot_sac(:, counter_dir, cond_idx, :), ...
         cs_vm_tot_sac(:,   counter_dir, cond_idx, :)] = process_event_for_bins( ...
             sac_data_dir_sac, order_ang(counter_dir), cond_idx, ...
             force_conds, sac_crit, vel_bins, n_timepoints);

        % Offset
        [cs_rate_tot_offset(:, counter_dir, cond_idx, :), ...
         cs_vm_tot_offset(:,   counter_dir, cond_idx, :)] = process_event_for_bins( ...
             sac_data_dir_offset, order_ang(counter_dir), cond_idx, ...
             force_conds, sac_crit, vel_bins, n_timepoints);
    end
end

% ------------ Pack results (5th dim = event) ------------
res.rate_tot = cat(5, cs_rate_tot_vis, cs_rate_tot_sac, cs_rate_tot_offset);
res.vm_tot   = cat(5, cs_vm_tot_vis,   cs_vm_tot_sac,   cs_vm_tot_offset);
res.vel_bins = vel_bins;

end  % <-- end of main function


% ============================================================
% Local function 
% ============================================================
function [rate_bin, vm_bin] = process_event_for_bins( ...
    sac_data_dir, dir_idx, cond_idx, ...
    force_conds, sac_crit, vel_edges, n_timepoints)
% Outputs:
%   rate_bin [n_timepoints x n_vel_bins]
%   vm_bin   [n_timepoints x n_vel_bins]

n_vel_bins = numel(vel_edges) - 1;

% ---- Pull per-trial fields (possibly unequal lengths) ----
task_cond = sac_data_dir(dir_idx).sac.task_cond(:);
task_tag  = sac_data_dir(dir_idx).sac.tag(:);
tgt_cond  = sac_data_dir(dir_idx).sac.tgt_cond(:);

amp   = sac_data_dir(dir_idx).sac.eye_amp(:);
vmax  = sac_data_dir(dir_idx).sac.eye_vm_max(:);
t_off = sac_data_dir(dir_idx).sac.time_offset(:);
t_max = sac_data_dir(dir_idx).sac.time_vmax(:);
t_on  = sac_data_dir(dir_idx).sac.time_onset(:);

raster_full = sac_data_dir(dir_idx).CS;      % [time x N?]
vm_full     = sac_data_dir(dir_idx).eye_vm;  % [time x N?]

% ---- Harmonize to a common number of trials across ALL arrays ----
lengths = [numel(task_cond), numel(task_tag), numel(tgt_cond), ...
           numel(amp), numel(vmax), numel(t_off), numel(t_max), numel(t_on), ...
           size(raster_full,2), size(vm_full,2)];
N_trials = min(lengths);  

% Truncate everything consistently to 1:N_trials
task_cond = task_cond(1:N_trials);
task_tag  = task_tag(1:N_trials);
tgt_cond  = tgt_cond(1:N_trials);
amp       = amp(1:N_trials);
vmax      = vmax(1:N_trials);
t_off     = t_off(1:N_trials);
t_max     = t_max(1:N_trials);
t_on      = t_on(1:N_trials);

raster = raster_full(:, 1:N_trials);
vm_tr  = vm_full(:,     1:N_trials);

% ---- Condition mask (High vs Low) ----
isHigh = ismember([task_cond task_tag tgt_cond], force_conds(1,:),'rows');
isLow  = ismember([task_cond task_tag tgt_cond], force_conds(2,:),'rows');
if     cond_idx == 1, in_cond = isHigh;
elseif cond_idx == 2, in_cond = isLow;
else, error('cond_idx must be 1 (High) or 2 (Low)'); 
end

% ---- Validity (MS/MA polygons, timing) ----
dec_ms = (t_off - t_max) * 1e3;
acc_ms = (t_max - t_on) * 1e3;

[in_ms,on_ms] = inpolygon(log10(amp), log10(vmax),    sac_crit.ms(1,:), sac_crit.ms(2,:));
acc_dec       = dec_ms ./ max(acc_ms, eps);
[in_ma,on_ma] = inpolygon(log10(amp), log10(acc_dec), sac_crit.ma(1,:), sac_crit.ma(2,:));

valid = (in_ms | on_ms) & (in_ma | on_ma) & (acc_ms > 0) & (dec_ms > 0);

% Keep only valid trials in the  condition
keep = valid & in_cond;

vmax_k = vmax(keep);
raster = raster(:, keep);
vm_tr  = vm_tr(:,   keep);

% ---- Outputs ----
rate_bin = nan(n_timepoints, n_vel_bins);
vm_bin   = nan(n_timepoints, n_vel_bins);

% ---- Weighted averages per velocity bin  ----
for vb = 1:n_vel_bins
    vmin = vel_edges(vb); vmax_edge = vel_edges(vb+1);
    if vb < n_vel_bins
        in_bin = (vmax_k >= vmin) & (vmax_k <  vmax_edge);
    else
        in_bin = (vmax_k >= vmin) & (vmax_k <= vmax_edge);  % right-inclusive
    end
    ntr = sum(in_bin);
    if ntr == 0, continue; end

    w = ones(1, ntr) / ntr;                

    r_bin = raster(:, in_bin);             % [time x ntr]
    vm_b  = vm_tr(:, in_bin);

    rate_bin(:, vb) = ESN_smooth(r_bin * w.') * 1e3;  % Hz
    vm_bin(:,   vb) =             (vm_b  * w.' );
end
end

