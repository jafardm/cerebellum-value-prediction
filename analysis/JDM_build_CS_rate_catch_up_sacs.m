function JDM_build_CS_rate_catch_up_sacs(data_path)
JDM_params_funcs
load(fullfile(data_path,'population_data','behave_data'));
user = 'JDM'; loc = 'ctx';

path_data_monkey_sorted = params.path_data_monkey_sorted;
animal_list   = params.animal_list;
session_list  = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);
num_animal    = numel(session_list);

cs_on_data  = cell(0);
cs_cell_ids = cell(0);
cs_rate     = [];
cs_vm       = [];

N_selected_all = cell(num_animal,1);  
N_total_all    = cell(num_animal,1);   
N_amp_all      = cell(num_animal,1); 
N_tags_all     = cell(num_animal,1); 
N_lat_all      = cell(num_animal,1);   
N_lat_cond_all = cell(num_animal,1);  

for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path      = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess          = numel(current_sess_list);

    cs_on_data_sess  = cell(num_sess,1);
    cs_cell_ids_sess = cell(num_sess,1);
    cs_rate_sess     = cell(num_sess,1);
    cs_vm_sess       = cell(num_sess,1);

    N_selected_sess  = cell(num_sess,1);   
    N_total_sess     = cell(num_sess,1);   % <-- NEW
    N_amp_sess       = cell(num_sess,1);   
    N_tags_sess      = cell(num_sess,1); 

    current_sess_cr_ma = data_behav.sess_cr_ma{counter_animal};
    current_sess_cr_ms = data_behav.sess_cr_ms{counter_animal};

    parfor counter_sess = 1:num_sess
        current_sess = current_sess_list{counter_sess};
        disp(['     ' current_sess]);

        current_cr_ma = current_sess_cr_ma{counter_sess};
        current_cr_ms = current_sess_cr_ms{counter_sess};

        data_eye   = MAF_load_eye_traces(current_path,current_sess);
        units_info = MAF_extract_cell_metadata(current_path,current_sess);
        cell_types = units_info.cell_type;
        cell_ids   = units_info.cell_list;
        cs_cell_ids_ = cell_ids(ismember(cell_types,{'CS','PC'}));

        % --- UPDATED: return N_total as well
         [cs_on_data_,cs_rate_,cs_vm_,Nsel,Ntot,amp_selected,tags_selected, ...
         lat_selected,lat_cond] = par_loop_cs_on_data(...
         current_path,current_sess,data_eye,cs_cell_ids_, ...
            current_cr_ma,current_cr_ms);

        cs_cell_ids_sess{counter_sess} = cs_cell_ids_;
        cs_on_data_sess{counter_sess}  = cs_on_data_;
        cs_rate_sess{counter_sess}     = cs_rate_;
        cs_vm_sess{counter_sess}       = cs_vm_;
        N_selected_sess{counter_sess}  = Nsel;  
        N_total_sess{counter_sess}     = Ntot;   
        N_amp_sess{counter_sess}       = amp_selected; 
        N_tags_sess{counter_sess}      = tags_selected;
        N_lat_sess{counter_sess}       = lat_selected;  
        N_lat_cond_sess{counter_sess}  = lat_cond;       
    end

    cs_on_data  = [cs_on_data;  vertcat(cs_on_data_sess{:})];
    cs_cell_ids = [cs_cell_ids; vertcat(cs_cell_ids_sess{:})];

    if counter_animal == 1
        cs_rate = cat(1, cs_rate_sess{:});
        cs_vm   = cat(1, cs_vm_sess{:});
    else
        cs_rate = cat(1, cs_rate, cat(1, cs_rate_sess{:}));
        cs_vm   = cat(1, cs_vm,   cat(1, cs_vm_sess{:}));
    end

    N_selected_all{counter_animal} = N_selected_sess;  
    N_total_all{counter_animal}    = N_total_sess;     
    N_amp_all{counter_animal}      = N_amp_sess;
    N_tags_all{counter_animal}     = N_tags_sess;
    N_lat_all{counter_animal}      = N_lat_sess;       % <-- NEW
    N_lat_cond_all{counter_animal} = N_lat_cond_sess;  % <-- NEW
end

% ---------- META ----------
file_name = 'cs_on_rate_catch_up_sac_sac';
meta = struct();
meta.file_name  = file_name;
meta.created_at = datestr(now,'yyyy-mm-dd HH:MM:SS');
meta.user = user; meta.loc = loc;
meta.animals      = animal_list;
meta.session_list = session_list;
meta.n_cells      = numel(cs_cell_ids);
meta.n_timepoints = 500; meta.n_dirs = 8;
meta.force_conds  = params.reward_tag1_conds(1:4,:);  
meta.order_ang    = 1:8;
meta.selection.ang_thresh_deg = 30;
meta.selection.amp_thresh_deg = 2.0;
meta.selection.pos_thresh_deg = 0.5;

N_selected_sessions = struct();
N_selected_sessions.animals      = animal_list;
N_selected_sessions.session_list = session_list;
N_selected_sessions.N            = N_selected_all;   
N_selected_sessions.N_total      = N_total_all;       % <-- NEW

% ---------- SAVE ----------
out_path = fullfile(data_path,'population_data');
if ~exist(out_path,'dir'), mkdir(out_path); end

save(fullfile(out_path, sprintf('%s.mat',file_name)), ...
     'cs_on_data','cs_cell_ids','cs_rate','cs_vm','meta', ...
     'N_selected_sessions','N_amp_all','N_tags_all', ...
     'N_lat_all','N_lat_cond_all','-v7.3');   % <-- NEW

save(fullfile(out_path, sprintf('%s_meta.mat',file_name)), ...
     'meta','N_selected_sessions');
end


%% par_loop_cs_on_data
function [cs_on_data,cs_rate,cs_vm,N_selected_sess,N_total_sess, ...
          amp_selected_sess,tag_selected_sess,lat_selected_sess,lat_cond_sess] = ...
    par_loop_cs_on_data(current_path,current_sess,data_eye,cell_ids, ...
                        current_cr_ma,current_cr_ms)

num_cs = numel(cell_ids);
cs_on_data = cell(num_cs,1);

cs_rate = nan(num_cs,500,8,4);
cs_vm   = nan(num_cs,500,8,4);

N_selected_candidates   = cell(num_cs,1);
N_total_candidates      = cell(num_cs,1);   
amp_selected_candidates = cell(num_cs,1);   
tag_selected_candidates = cell(num_cs,1);   
lat_selected_candidates = cell(num_cs,1);   % <-- NEW
lat_cond_candidates     = cell(num_cs,1);   % <-- NEW

sac_crit.ma = current_cr_ma;
sac_crit.ms = current_cr_ms;

parfor counter_cs = 1:num_cs
    current_cs = cell_ids{counter_cs};
    data = MAF_load_cell(current_path,current_sess,current_cs);
    data_recordings = data.data_recordings;
    cs_on_data{counter_cs} = MAF_CS_on_analysis_fr([data_recordings.eye], ...
        [data_recordings.Neural_Data],[1,4,6,7],[40,85],[-70,30]);

    res_rate = load_cs_data_tuned(current_path,current_sess,current_cs,1:8,data_eye);

    cs_rate(counter_cs,:,:,:) = res_rate.cs_rate_sac;  
    cs_vm(counter_cs,:,:,:)   = res_rate.cs_vm_sac;

    if isfield(res_rate,'N_selected')
        N_selected_candidates{counter_cs} = res_rate.N_selected;
    end
    if isfield(res_rate,'N_total')
        N_total_candidates{counter_cs} = res_rate.N_total;
    end
    if isfield(res_rate,'amp_selected') && ~isempty(res_rate.amp_selected)
        amp_selected_candidates{counter_cs} = res_rate.amp_selected;
    end
    if isfield(res_rate,'amp_tags') && ~isempty(res_rate.amp_tags)
        tag_selected_candidates{counter_cs} = res_rate.amp_tags;
    end
    if isfield(res_rate,'catchup_latency') && ~isempty(res_rate.catchup_latency)
        lat_selected_candidates{counter_cs} = res_rate.catchup_latency;
        lat_cond_candidates{counter_cs}     = res_rate.catchup_latency_cond;
    end
end

% --- Post-process ---
N_selected_sess = get_first_nonempty(N_selected_candidates,[8,4]);
N_total_sess    = get_first_nonempty(N_total_candidates,[8,4]);

if any(~cellfun(@isempty,amp_selected_candidates))
    amp_selected_sess = vertcat(amp_selected_candidates{:});
else
    amp_selected_sess = [];
end

if any(~cellfun(@isempty,tag_selected_candidates))
    tag_selected_sess = vertcat(tag_selected_candidates{:});
else
    tag_selected_sess = [];
end

if any(~cellfun(@isempty,lat_selected_candidates))
    lat_selected_sess = vertcat(lat_selected_candidates{:});
    lat_cond_sess     = vertcat(lat_cond_candidates{:});
else
    lat_selected_sess = [];
    lat_cond_sess     = [];
end
end

function out = get_first_nonempty(cellarr,sz)
if any(~cellfun(@isempty,cellarr))
    idx = find(~cellfun(@isempty,cellarr),1,'first');
    out = cellarr{idx};
else
    out = zeros(sz);
end
end

%% load_cs_data_tuned
function res = load_cs_data_tuned(current_path,current_sess,current_cell,order_ang,data_eye)
JDM_params_funcs;

force_conds  = params.reward_tag1_conds(1:4,:);  
n_timepoints = 500;
n_dirs       = 8;
n_conds      = 4;

cs_rate_tot_sac = nan(n_timepoints, n_dirs, n_conds);
cs_vm_tot_sac   = nan(n_timepoints, n_dirs, n_conds);

% --- trial counters ---
N_catchup = zeros(n_dirs, n_conds);   % catch-up trials (as before)
N_total   = zeros(n_dirs, n_conds);   % all tuned trials (NEW)

% --- collectors ---
all_amp_deg   = [];     
all_cond_idx  = [];     
all_dir_idx   = [];     
all_tags      = [];     
all_latencies = [];     
all_lat_cond  = [];     

data       = MAF_load_cell(current_path,current_sess,current_cell);
eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);

neural_sac_data = JDM_combine_dataset_sac([data.data_recordings],[data.rec_info],{'onset'},250,1:10);
eye_vm          = MAF_combine_vm_traces(eye_traces,{'onset'},250,1:10);

neur_data_sac.eye_vm    = eye_vm.onset.eye_vm;
neur_data_sac.CS        = neural_sac_data.onset.CS;   
sac_data_sac.sac        = neural_sac_data.onset.sac;
sac_data_sac.fix        = neural_sac_data.onset.fix;
sac_data_sac.sac.has_cs = neural_sac_data.onset.has_cs;

sac_data_dir_sac = MAF_bin_by_ang(neur_data_sac, sac_data_sac, 'onset');

ang_thresh_deg = 30;
amp_thresh_deg = 2.2;
pos_thresh     = 0.5;
letency_tol    = 500;

res.after_tuned_small(n_dirs, n_conds).pairs = [];  

for counter_dir = 1:n_dirs
    dir_idx = order_ang(counter_dir);

    CS_mat = sac_data_dir_sac(dir_idx).CS; 
    VM_mat = sac_data_dir_sac(dir_idx).eye_vm;
    S      = sac_data_dir_sac(dir_idx).sac;

    st  = (S.eye_px_onset(:)  + 1i*S.eye_py_onset(:));
    en  = (S.eye_px_offset(:) + 1i*S.eye_py_offset(:));

    for cond_idx = 1:n_conds
        task_cond = S.task_cond;  tgt_cond = S.tgt_cond;  task_tag = S.tag;
        jmp_cond  = S.jump_cond;  rwd_cond = S.rew_cond;

        if cond_idx <= 2
            px_ = S.cue_x_high_rew;  py_ = S.cue_y_high_rew;
        else
            px_ = S.cue_x_low_rew;   py_ = S.cue_y_low_rew;
        end
        tgt = px_(:) + 1i*py_(:);

        sacs_tuned = ismember([task_cond' task_tag' tgt_cond' jmp_cond' rwd_cond'], ...
                              force_conds(cond_idx,:), 'rows');
        idx_tuned = find(sacs_tuned(:));

        % --- count all tuned trials ---
        N_total(counter_dir, cond_idx) = numel(idx_tuned);

        keep_idx     = [];
        pairs        = [];
        next_amp_deg = [];
        ang_err_deg  = [];
        dist_improve = [];

        for k = idx_tuned.'
            cand_rel = find(abs(st(k+1:end) - en(k)) < pos_thresh, 1, 'first');
            if isempty(cand_rel), continue; end
            m = k + cand_rel;

            next_vec   = en(m) - st(m);
            to_tgt_vec = tgt(m) - st(m);
            if ~isfinite(next_vec) || ~isfinite(to_tgt_vec) || next_vec==0 || to_tgt_vec==0
                continue
            end

            amp   = abs(next_vec);
            c     = real(next_vec*conj(to_tgt_vec)) / (abs(next_vec)*abs(to_tgt_vec));
            c     = max(-1, min(1, c));
            theta = acosd(c);

            d_before = abs(tgt(m) - st(m));
            d_after  = abs(tgt(m) - en(m));
            moved_closer = d_after < d_before;

            if (theta <= ang_thresh_deg) && moved_closer && (amp < amp_thresh_deg) ...
                           && ismember(task_tag(m),[9 10])
                if  ~isnan(S.time_onset(m)) && ~isnan(S.time_offset(k))
                    latency = (S.time_onset(m) - S.time_offset(k))*1e3; % ms
                    
                    if latency <= letency_tol
                        keep_idx     = [keep_idx; m];
                        pairs        = [pairs; [k m]];
                        next_amp_deg = [next_amp_deg; amp];
                        ang_err_deg  = [ang_err_deg; theta];
                        dist_improve = [dist_improve; d_before - d_after];
                        
                        all_latencies = [all_latencies; latency];
                        all_lat_cond  = [all_lat_cond; cond_idx];
                        
                        all_amp_deg  = [all_amp_deg;  amp];
                        all_cond_idx = [all_cond_idx; cond_idx];
                        all_dir_idx  = [all_dir_idx;  counter_dir];
                        all_tags     = [all_tags;     task_tag(m)];
                    end
                end
            end

        end

        % --- count catch-up trials ---
        N_catchup(counter_dir, cond_idx) = numel(keep_idx);

        % --- save condition-specific results ---
        if ~isempty(keep_idx)
            cs_mean_raster = mean(double(CS_mat(:, keep_idx)), 2, 'omitnan'); 
            cs_rate_hz     = ESN_smooth(cs_mean_raster) * 1e3;         
            cs_rate_tot_sac(:, counter_dir, cond_idx) = cs_rate_hz;
            cs_vm_tot_sac(:,  counter_dir, cond_idx) = mean(VM_mat(:, keep_idx), 2, 'omitnan');
        end

        res.after_tuned_small(counter_dir,cond_idx).pairs        = pairs;
        res.after_tuned_small(counter_dir,cond_idx).next_idx     = keep_idx;
        res.after_tuned_small(counter_dir,cond_idx).next_amp_deg = next_amp_deg;
        res.after_tuned_small(counter_dir,cond_idx).ang_err_deg  = ang_err_deg;
        res.after_tuned_small(counter_dir,cond_idx).dist_improve = dist_improve;
        res.after_tuned_small(counter_dir,cond_idx).N            = numel(keep_idx);
    end
end

% --- pack outputs ---
res.cs_rate_sac          = cs_rate_tot_sac;     
res.cs_vm_sac            = cs_vm_tot_sac;       
res.N_selected           = N_catchup;    % catch-up only (as before)
res.N_total              = N_total;      %  all tuned trials
res.catchup_latency      = all_latencies;
res.catchup_latency_cond = all_lat_cond;
res.amp_selected         = all_amp_deg;
res.amp_cond_idx         = all_cond_idx;
res.amp_dir_idx          = all_dir_idx;
res.amp_tags             = all_tags;

res.meta.n_timepoints    = n_timepoints;
res.meta.ang_thresh_deg  = ang_thresh_deg;
res.meta.amp_thresh_deg  = amp_thresh_deg;
res.meta.pos_thresh_deg  = pos_thresh;
res.meta.order_ang       = order_ang;
res.meta.force_conds     = force_conds;
end

