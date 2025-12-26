function JDM_build_Behave_by_trial
JDM_params_funcs
user = 'JDM';
loc  = 'ctx';
tags = 1:10;

path_data_monkey_sorted = params.path_data_monkey_sorted;
animal_list = params.animal_list;
session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);

num_animal = numel(session_list);
transitions = {
        1, [1 2];  % HH → HH|HL
        2, [1 2];  % HL → HH|HL
        3, [3 4];  % LL → LL|LH
        4, [3 4];  % LH → LL|LH
    };

data_behav.sac_rt    = cell(num_animal,1);
data_behav.vmax      = cell(num_animal,1);
data_behav.sac_error = cell(num_animal,1);

for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess = numel(current_sess_list);

    sess_sac_rt   = cell(num_sess,1);
    sess_vmax_sac = cell(num_sess,1);
    sac_error     = cell(num_sess,1);

      for counter_sess = 1:num_sess
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
          tags_tr    = [sacs_tr.tag]';
          data_rec   = data.data_recordings;   
          eye        = [data_rec.eye];
          sac        = [eye.sac];
          sac_tag_   = [sac.tag]';
          sac_vmax_  = [sac.eye_vm_max];
          sac_amp_   = [sac.eye_amp]';

           assert(all(tags_tr == sac_tag_));
    
           tgt_cond_  = [sac.tgt_cond]';
           task_cond_ = [sac.task_cond]';
           rwd_cond_  = [sac.rew_cond]';
           sac_rt_    = (([sac.time_onset] - [sac.time_visual])*1e3)';
           sac_end_x_ = [sac.eye_px_offset];
           sac_end_y_ = [sac.eye_py_offset];
    
           sac_st_x_  = [sac.eye_px_onset];
           sac_st_y_  = [sac.eye_py_onset];
    
           vis_end_x_ = [sac.vis_px_offset];
           vis_end_y_ = [sac.vis_py_offset];
           sac_vec_   = (sac_end_x_ - sac_st_x_) + 1j*(sac_end_y_ - sac_st_y_);
           gl_vec_    = (vis_end_x_ - sac_st_x_) + 1j*(vis_end_y_ - sac_st_y_);
           idx_rxt    = sac_rt_ > 0 & sac_rt_ < 600;
           tag_1      = sac_tag_ == 1;
           sac_amp_idx = sac_amp_ > 3.5;
           assert(length(sac_rt_) == length(sac_tag_));
    
            res = MAF_find_sac_outliers(sac);
            % get the outlayers indices
            ind_selected = (ismember(sac_tag_,tags) & not(res.rm_ind)' & idx_rxt); % not outlies + tag 1-10
            idx_tuned = ind_selected & tag_1 & sac_amp_idx;
        
            sess_task_cond = task_cond_(idx_tuned);
            sess_tag = tag_1(idx_tuned);
            sess_tgt = tgt_cond_(idx_tuned);
            sess_rwd = rwd_cond_(idx_tuned);
            sess_vmax = sac_vmax_(idx_tuned);
            sess_sac_rt_ = sac_rt_(idx_tuned);
            [px_error, py_error] = compute_enpoint_error(gl_vec_(idx_tuned), sac_vec_(idx_tuned));
            sess_sac_error = sqrt(px_error.^2 + py_error.^2 );
    
            conditions = [sess_task_cond, sess_tag sess_tgt, sess_rwd ];

        
            trial_labels = nan(size(conditions, 1), 1);
            trial_labels(ismember(conditions, [1,1,1,1], 'rows')) = 1; % HH
            trial_labels(ismember(conditions, [1,1,1,0], 'rows')) = 2; % HL
            trial_labels(ismember(conditions, [1,1,0,0], 'rows')) = 3; % LL
            trial_labels(ismember(conditions, [1,1,0,1], 'rows')) = 4; % LH
    
    
            noftransit = size(transitions, 1);
            rt_sac_    = cell(noftransit,1);
            vmax_sac_  = cell(noftransit,1);
            sac_error_ = cell(noftransit,1);
    
    
         for t = 1:noftransit
            from_label = transitions{t,1};   % e.g., HH
            to_labels  = transitions{t,2};   % e.g., HH or HL
        
            prev = trial_labels(1:end-1);    % previous trial labels
            curr = trial_labels(2:end);      % current trial labels
        
            valid_trans = find(prev == from_label & ismember(curr, to_labels));
    
            sac_rt_temp = [];
            vmax_temp   = [];
            error_temp  = [];
          
            for idx = valid_trans'
                t1 = idx;
                t2 = idx + 1;
                sac_rt_temp  = [sac_rt_temp; sess_sac_rt_(t1) - sess_sac_rt_(t2)];
                vmax_temp    = [vmax_temp; sess_vmax(t1) - sess_vmax(t2)];
                error_temp   = [error_temp; sess_sac_error(t1) - sess_sac_error(t2)];

            end 
            rt_sac_{t,1}    = sac_rt_temp;
            vmax_sac_{t,1}  = vmax_temp;
            sac_error_{t,1} = error_temp;
    
         end  
         sess_sac_rt{counter_sess}   = rt_sac_;
         sess_vmax_sac{counter_sess} = vmax_sac_;
         sac_error{counter_sess}     = sac_error_;
       end
% Concatenate across sessions per transition (i.e., across rows)
data_behav.sac_rt{counter_animal} = ...
    cellfun(@(varargin) vertcat(varargin{:}), sess_sac_rt{:}, 'UniformOutput', false);
data_behav.vmax{counter_animal}   =...
    cellfun(@(varargin) vertcat(varargin{:}), sess_vmax_sac{:}, 'UniformOutput', false);
data_behav.sac_error{counter_animal} =...
    cellfun(@(varargin) vertcat(varargin{:}), sac_error{:}, 'UniformOutput', false);
              
end

save_path = 'C:\Users\Jafar\Documents\reward\';
file_name = 'Behave_trial_by_trial';
save(fullfile(save_path,'population_data',sprintf('%s.mat',file_name)), 'data_behav', '-v7.3');
end

%% compute_enpoint_error
function [px_error, py_error] = compute_enpoint_error(gl_vec, sac_vec)
    angles_deg_wrapped = wrapTo360(rad2deg(angle(gl_vec)));
    N = length(angles_deg_wrapped);
    rot_eye_px_offset = zeros(1, N);
    rot_eye_py_offset = zeros(1, N);
    rot_goal_x = zeros(1, N);
    rot_goal_y = zeros(1, N);

    for ii = 1:N
        theta = angles_deg_wrapped(ii);
        R = [cosd(-theta), -sind(-theta); sind(-theta), cosd(-theta)];
        rot_sac  = R * [real(sac_vec(ii)); imag(sac_vec(ii))];
        rot_goal = R * [real(gl_vec(ii)); imag(gl_vec(ii))];
        rot_eye_px_offset(ii) = rot_sac(1);
        rot_eye_py_offset(ii) = rot_sac(2);
        rot_goal_x(ii) = rot_goal(1);
        rot_goal_y(ii) = rot_goal(2);
    end

    px_error = rot_eye_px_offset - rot_goal_x;
    py_error = rot_eye_py_offset - rot_goal_y;
end
