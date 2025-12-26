function JDM_build_CS_rate_by_trial (data_path)
JDM_params_funcs
load(fullfile(data_path,'population_data','behave_data'));
user = 'JDM';
loc  = 'ctx';

path_data_monkey_sorted = params.path_data_monkey_sorted;
animal_list = params.animal_list;
session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);

num_animal = numel(session_list);

cs_cell_ids = cell(0);
cs_rate     = [];


for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess = numel(current_sess_list);


    cs_cell_ids_sess = cell(num_sess,1);
    cs_rate_sess     = cell(num_sess,1);

    %parfor
   parfor counter_sess = 1:num_sess
        current_sess = current_sess_list{counter_sess};
        disp(['     ' current_sess]);
        % load eye_traces

        units_info = MAF_extract_cell_metadata(current_path,current_sess);
        cell_types = units_info.cell_type;
        cell_ids   = units_info.cell_list;
        cs_cell_ids_ = cell_ids(ismember(cell_types,{'CS','PC'}));

        cs_rate_ = par_loop_cs_on_data(current_path,current_sess,cs_cell_ids_);
        cs_cell_ids_sess{counter_sess} = cs_cell_ids_;
        cs_rate_sess{counter_sess}     = cs_rate_;

    end

    cs_cell_ids = [cs_cell_ids; vertcat(cs_cell_ids_sess{:})];

    if counter_animal == 1
        cs_rate = cat(1,cs_rate_sess{:});
    else
        cs_rate = cat(1,cs_rate, cat(1,cs_rate_sess{:}));
    end
end
file_name = 'CS_rate_by_trial';
save(fullfile(data_path,'population_data',sprintf('%s.mat',file_name)),...
   'cs_cell_ids','cs_rate');
end

%% par_loop_cs_on_data
function cs_rate= par_loop_cs_on_data(current_path,current_sess,cell_ids)
    num_cs = numel(cell_ids);

    cs_rate = nan(num_cs,500,4,2);
   %parfor
   parfor counter_cs = 1:num_cs
        current_cs = cell_ids{counter_cs};    
        res_rate = csrate_trail_data_tuned(current_path,current_sess,current_cs);
  
        cs_rate(counter_cs,:,:,:)   = res_rate.rate_tot;
    end
end
%% load_cs_data_tuned
function res = csrate_trail_data_tuned(current_path, current_sess, current_cell)
   JDM_params_funcs;
   tags = 1;
   n_timepoints = 500;

   data = MAF_load_cell(current_path, current_sess, current_cell);
   neural_sac_data = JDM_combine_dataset_sac([data.data_recordings],[data.rec_info],{'visual','onset'},250,tags);
    
    % --- Visual Condition ---
    task_cond_vis = neural_sac_data.visual.sac.task_cond';
    tgt_cond_vis  = neural_sac_data.visual.sac.tgt_cond';
    rwd_cond_vis  = neural_sac_data.visual.sac.rew_cond';
    conditions_vis = [task_cond_vis, tgt_cond_vis, rwd_cond_vis];

    trial_labels_vis = nan(size(conditions_vis, 1), 1);
    trial_labels_vis(ismember(conditions_vis, [1,1,1], 'rows')) = 1; % HH
    trial_labels_vis(ismember(conditions_vis, [1,1,0], 'rows')) = 2; % HL
    trial_labels_vis(ismember(conditions_vis, [1,0,0], 'rows')) = 3; % LL
    trial_labels_vis(ismember(conditions_vis, [1,0,1], 'rows')) = 4; % LH

    has_cs_vis   = neural_sac_data.visual.has_cs';
    cs_raster_vis = neural_sac_data.visual.CS;  % [time x trials]

    rate_diff_vis = nan(n_timepoints,4);  % averaged trial-wise rate differences

    transitions = {
        1, [1 2];  % HH → HH|HL
        2, [1 2];  % HL → HH|HL
        3, [3 4];  % LL → LL|LH
        4, [3 4];  % LH → LL|LH
    };

    for t = 1:size(transitions, 1)
        from_label = transitions{t,1};
        to_labels  = transitions{t,2};

        prev = trial_labels_vis(1:end-1);
        curr = trial_labels_vis(2:end);

        valid_trans = find(prev == from_label & ismember(curr, to_labels));

        if isempty(valid_trans)
            continue;
        end

        diffs = [];

        for idx = valid_trans'
            t1 = idx;
            t2 = idx + 1;

            if has_cs_vis(t1) && has_cs_vis(t2)
                cs1 = cs_raster_vis(:,t1);
                cs2 = cs_raster_vis(:,t2);
                rate_diff = ESN_smooth(cs1 - cs2) * 1e3;
                diffs = [diffs, rate_diff];
            end
        end

        rate_diff_vis(:,t) = mean(diffs,2);
    end
    
    outlayers = MAF_find_sac_outliers(neural_sac_data.onset.sac);
    idx_valid = not(outlayers.rm_ind)';
    % --- Saccade Condition ---
    task_cond = neural_sac_data.onset.sac.task_cond';
    tgt_cond  = neural_sac_data.onset.sac.tgt_cond';
    rwd_cond  = neural_sac_data.onset.sac.rew_cond';

    conditions = [task_cond, tgt_cond, rwd_cond];
    conditions = conditions(idx_valid,:);

    trial_labels = nan(size(conditions, 1), 1);
    trial_labels(ismember(conditions, [1,1,1], 'rows')) = 1; % HH
    trial_labels(ismember(conditions, [1,1,0], 'rows')) = 2; % HL
    trial_labels(ismember(conditions, [1,0,0], 'rows')) = 3; % LL
    trial_labels(ismember(conditions, [1,0,1], 'rows')) = 4; % LH

    has_cs    = neural_sac_data.onset.has_cs(idx_valid)';
    cs_raster = neural_sac_data.onset.CS(:,idx_valid);
    
    noftransit = size(transitions, 1);
    rate_diff_sac = nan(n_timepoints,4);
 
    for t = 1:noftransit
        from_label = transitions{t,1};   % e.g., HH
        to_labels  = transitions{t,2};   % e.g., HH or HL
    
        prev = trial_labels(1:end-1);    % previous trial labels
        curr = trial_labels(2:end);      % current trial labels
    
        valid_trans = find(prev == from_label & ismember(curr, to_labels));

        diffs = [];  
        for idx = valid_trans'
            t1 = idx;
            t2 = idx + 1;
    
            if has_cs(t1) && has_cs(t2)
                cs1 = cs_raster(:,t1);                   % CS trace of trial 1
                cs2 = cs_raster(:,t2);                   % CS trace of trial 2
                rate_diff = ESN_smooth(cs1 - cs2) * 1e3; % smooth and convert to Hz
                diffs = [diffs, rate_diff];              % store in [time x trials] matrix 

            end
        end
    
        rate_diff_sac(:,t) = mean(diffs,2); % average across trials for this transition
    
    end

    % --- Output ---
    res.rate_tot = cat(3, rate_diff_vis, rate_diff_sac);  % [time x 4 conditions x 2 epochs] 
end