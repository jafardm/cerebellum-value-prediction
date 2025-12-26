%% Build behave dataset
function JDM_build_stats_endpoint_error_corr_sac
    JDM_params_funcs
    user = 'JDM';
    loc  = 'ctx';
    tags = 1:10;

    % 4 reward conditions
    sac_conds = [
      1, 4, 1, 1, 1, 1; % HH
      1, 4, 1, 1, 0, 1; % HL
      1, 4, 0, 1, 0, 1; % LL
      1, 4, 0, 1, 1, 1];% LH

    path_data_monkey_sorted = params.path_data_monkey_sorted;
    animal_list = params.animal_list;
    session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);

    num_animal = numel(session_list);

    data_behav.var_x = cell(num_animal,1);
    data_behav.var_y = cell(num_animal,1);
    data_behav.sac_conds = sac_conds; % store condition definitions

    num_cond = size(sac_conds,1);

    for counter_animal = 1:num_animal
        disp(animal_list{counter_animal});
        current_path = path_data_monkey_sorted{counter_animal};
        current_sess_list = session_list{counter_animal};
        num_sess = numel(current_sess_list);

        % one column per condition
        var_x_   = nan(num_sess,num_cond);
        var_y_   = nan(num_sess,num_cond);

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

            tgt_cond_  = [sac.tgt_cond]';
            task_cond_ = [sac.task_cond]';
            jump_cond_ = [sac.jump_cond]';
            rwd_cond_  = [sac.rew_cond]';

            sac_rt_    = (([sac.time_onset] - [sac.time_visual])*1e3)';

            assert(all(tags_tr == sac_tag_));

            sac_end_x_ = [sac.eye_px_offset];
            sac_end_y_ = [sac.eye_py_offset];

            sac_st_x_  = [sac.eye_px_onset];
            sac_st_y_  = [sac.eye_py_onset];

            vis_end_x_ = [sac.vis_px_offset];
            vis_end_y_ = [sac.vis_py_offset];

            sac_vec_   = (sac_end_x_ - sac_st_x_) + 1j*(sac_end_y_ - sac_st_y_);
            gl_vec_    = (vis_end_x_ - sac_st_x_) + 1j*(vis_end_y_ - sac_st_y_);

            idx_rxt = sac_rt_ > 0 & sac_rt_ < 600;

            assert(length(sac_rt_) == length(sac_tag_));

            res = MAF_find_sac_outliers(sac);
            ind_selected = (ismember(sac_tag_,tags)' & ~res.rm_ind)'; 

            % ---- loop over conditions ----
            for c = 1:num_cond
                cond_mask = ismember([task_cond_ sac_tag_ tgt_cond_ jump_cond_ rwd_cond_ double(idx_rxt)], ...
                                     sac_conds(c,:), 'rows') & ind_selected;

                if any(cond_mask)
                    [px_error, py_error] = compute_enpoint_error(gl_vec_(cond_mask), sac_vec_(cond_mask));
                    var_x_(counter_sess,c) = std(px_error,'omitnan');
                    var_y_(counter_sess,c) = std(py_error,'omitnan');
                end
            end
        end

        % store results for this animal
        data_behav.var_x{counter_animal}  = var_x_; % [num_sessions × 4]
        data_behav.var_y{counter_animal}  = var_y_;
    end

    save_path = 'C:\Users\Jafar\Documents\reward\';
    file_name = 'stat_data_enderror_corr_sac';
    save(fullfile(save_path,'population_data',sprintf('%s.mat',file_name)), ...
        'data_behav', '-v7.3');
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
