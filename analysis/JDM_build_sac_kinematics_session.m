function JDM_build_sac_kinematics_session(data_path)
    JDM_params_funcs
    fprintf('\nLoading Data..........\n');
    load(fullfile(data_path,'population_data','behave_data_session'));  

    % [idx_tsk ,idx_tag, idx_tgt, idx_jmp, idx_rew, rxt]
    tag4_conds   = [params.reward_tag4_conds ones(4,1)];
    % [idx_tsk ,idx_tag, idx_tgt, idx_rew, rxt]
    tag1_conds   = [params.reward_tag1_conds(:,[1,2,3,5]) ones(4,1)];
    % [idx_tsk ,idx_ch, idx_tag, rxt]
    choice_conds = [params.reward_choice_conds ones(2,1) ones(2,1)];

    n_forc_cond  = size(tag1_conds,1);
    n_conds_tag1 = n_forc_cond + size(choice_conds,1);
    n_conds_tag4 = size(tag4_conds,1);
    subj_names   = params.animal_list;
    num_subjects = size(subj_names,2);

    % Top-level: {num_subjects,1}, then {nSess,1}, then {nConds,1}
    tag1_vigor            = cell(num_subjects, 1);
    tag1_sac_RT           = cell(num_subjects, 1);
    tag1_max_vel          = cell(num_subjects, 1);
    tag1_end_point_error  = cell(num_subjects, 1);
    tag1_eye_vm_vis       = cell(num_subjects, 1);
    tag1_eye_vm_sac       = cell(num_subjects, 1);

    tag4_vigor            = cell(num_subjects, 1);
    tag4_sac_RT           = cell(num_subjects, 1);
    tag4_max_vel          = cell(num_subjects, 1);
    tag4_end_point_error  = cell(num_subjects, 1);
    tag4_eye_vm_vis       = cell(num_subjects, 1);
    tag4_eye_vm_sac       = cell(num_subjects, 1);

    for ii = 1:num_subjects
        fprintf('Processing subject %d of %d...\n', ii, num_subjects);

        % Each of these is nSess×1 cell for this subject
        tot_sac_amp    = data_behav.sac_amp{ii};
        tot_sac_vel    = data_behav.sac_vel{ii};
        tot_sac_tag    = data_behav.sac_tag{ii};
        tot_tsk_cond   = data_behav.task_cond{ii};
        tot_tgt_cond   = data_behav.tgt_cond{ii};
        tot_jmp_cond   = data_behav.jump_cond{ii};
        tot_rwd_cond   = data_behav.rew_cond{ii};
        tot_rt         = data_behav.sac_rt{ii};
        tot_choice     = data_behav.choice{ii};
        tot_sac_vec    = data_behav.sac_vec{ii};
        tot_gl_vec     = data_behav.gl_vec{ii};
        tot_eye_vm_vis = data_behav.eye_vm_vis{ii};
        tot_eye_vm_sac = data_behav.eye_vm_sac{ii};

        n_sess = numel(tot_rt);   % e.g. 52

        % Per-subject containers: {nSess,1}
        tag1_vigor_subj           = cell(n_sess,1);
        tag1_sac_RT_subj          = cell(n_sess,1);
        tag1_max_vel_subj         = cell(n_sess,1);
        tag1_end_point_error_subj = cell(n_sess,1);
        tag1_eye_vm_vis_subj      = cell(n_sess,1);
        tag1_eye_vm_sac_subj      = cell(n_sess,1);

        tag4_vigor_subj           = cell(n_sess,1);
        tag4_sac_RT_subj          = cell(n_sess,1);
        tag4_max_vel_subj         = cell(n_sess,1);
        tag4_end_point_error_subj = cell(n_sess,1);
        tag4_eye_vm_vis_subj      = cell(n_sess,1);
        tag4_eye_vm_sac_subj      = cell(n_sess,1);

        % ==================== LOOP OVER SESSIONS ====================
        for ss = 1:n_sess
            fprintf('  Session %d of %d\n', ss, n_sess);

            % Pull this session’s vectors
            sac_amp    = tot_sac_amp{ss};
            sac_vel    = tot_sac_vel{ss};
            sac_tag    = tot_sac_tag{ss};
            tsk_cond   = tot_tsk_cond{ss};
            tgt_cond   = tot_tgt_cond{ss};
            jmp_cond   = tot_jmp_cond{ss};
            rwd_cond   = tot_rwd_cond{ss};
            rt         = tot_rt{ss};
            choice     = tot_choice{ss};
            sac_vec    = tot_sac_vec{ss};
            gl_vec     = tot_gl_vec{ss};
            eye_vm_vis = tot_eye_vm_vis{ss};  % [T × Ntr]
            eye_vm_sac = tot_eye_vm_sac{ss};  % [T × Ntr]

            % ---------- filters INSIDE session ----------
            idx_rxt    = rt > 0 & rt < 600;
            targ_tags  = ismember(sac_tag, [1 4 6 7]);
            fit_mask   = targ_tags & idx_rxt;

            if nnz(fit_mask) >= 5
                fit_obj = fitlm( ...
                    log10(sac_amp(fit_mask)), ...
                    log10(sac_vel(fit_mask)) );
            else
                % fallback: use all tagged trials even if RT is weird
                fallback_mask = targ_tags;
                fit_obj = fitlm( ...
                    log10(sac_amp(fallback_mask)), ...
                    log10(sac_vel(fallback_mask)) );
            end

            %% ---------------- Tag 1 (Primary) ----------------
            vigor         = cell(n_conds_tag1, 1);
            sac_RT        = cell(n_conds_tag1, 1);
            max_vel       = cell(n_conds_tag1, 1);
            end_point_err = cell(n_conds_tag1, 1);
            eye_vm_vis_out = cell(n_conds_tag1, 1);
            eye_vm_sac_out = cell(n_conds_tag1, 1);

            for cond_idx = 1:n_conds_tag1
                fprintf('    Tag1 Condition %d of %d\n', cond_idx, n_conds_tag1);

                if cond_idx <= n_forc_cond
                    cond_mat = [tsk_cond sac_tag tgt_cond rwd_cond double(idx_rxt)];
                    base_match = ismember(cond_mat, tag1_conds(cond_idx,:), 'rows');
                else
                    cond_mat = [tsk_cond choice sac_tag double(idx_rxt)];
                    base_match = ismember(cond_mat, ...
                                          choice_conds(cond_idx-n_forc_cond,:), 'rows');
                end

                istuned = base_match;

                amps_temp = log10(sac_amp(istuned));
                vels_temp = log10(sac_vel(istuned));

                if ~isempty(amps_temp)
                    log_pred_vel = predict(fit_obj, ...
                        table(amps_temp, 'VariableNames', {'x1'}));

                    vigor{cond_idx}          = vels_temp ./ log_pred_vel;
                    sac_RT{cond_idx}         = rt(istuned);
                    max_vel{cond_idx}        = sac_vel(istuned);
                    eye_vm_vis_out{cond_idx} = eye_vm_vis(:,istuned);
                    eye_vm_sac_out{cond_idx} = eye_vm_sac(:,istuned);

                    gl   = gl_vec(istuned);
                    sac  = sac_vec(istuned);
                    [px_err, py_err] = compute_enpoint_error(gl, sac);
                    end_point_err{cond_idx} = [px_err' py_err'];
                else
                    vigor{cond_idx}          = [];
                    sac_RT{cond_idx}         = [];
                    max_vel{cond_idx}        = [];
                    eye_vm_vis_out{cond_idx} = [];
                    eye_vm_sac_out{cond_idx} = [];
                    end_point_err{cond_idx}  = [];
                end
            end

            tag1_vigor_subj{ss}           = vigor;
            tag1_sac_RT_subj{ss}          = sac_RT;
            tag1_max_vel_subj{ss}         = max_vel;
            tag1_end_point_error_subj{ss} = end_point_err;
            tag1_eye_vm_vis_subj{ss}      = eye_vm_vis_out;
            tag1_eye_vm_sac_subj{ss}      = eye_vm_sac_out;

            %% ---------------- Tag 4 (Corrective) ----------------
            vigor         = cell(n_conds_tag4, 1);
            sac_RT        = cell(n_conds_tag4, 1);
            max_vel       = cell(n_conds_tag4, 1);
            end_point_err = cell(n_conds_tag4, 1);
            eye_vm_vis_out = cell(n_conds_tag4, 1);
            eye_vm_sac_out = cell(n_conds_tag4, 1);

            for cond_idx = 1:n_conds_tag4
                fprintf('    Tag4 Condition %d of %d\n', cond_idx, n_conds_tag4);

                cond_mat = [tsk_cond sac_tag tgt_cond jmp_cond rwd_cond double(idx_rxt)];
                base_match = ismember(cond_mat, tag4_conds(cond_idx,:), 'rows');
                istuned = base_match;

                amps_temp = log10(sac_amp(istuned));
                vels_temp = log10(sac_vel(istuned));

                if ~isempty(amps_temp)
                    log_pred_vel = predict(fit_obj, ...
                        table(amps_temp, 'VariableNames', {'x1'}));

                    vigor{cond_idx}          = vels_temp ./ log_pred_vel;
                    sac_RT{cond_idx}         = rt(istuned);
                    max_vel{cond_idx}        = sac_vel(istuned);
                    eye_vm_vis_out{cond_idx} = eye_vm_vis(:,istuned);
                    eye_vm_sac_out{cond_idx} = eye_vm_sac(:,istuned);

                    gl   = gl_vec(istuned);
                    sac  = sac_vec(istuned);
                    [px_err, py_err] = compute_enpoint_error(gl, sac);
                    end_point_err{cond_idx} = [px_err' py_err'];
                else
                    vigor{cond_idx}          = [];
                    sac_RT{cond_idx}         = [];
                    max_vel{cond_idx}        = [];
                    eye_vm_vis_out{cond_idx} = [];
                    eye_vm_sac_out{cond_idx} = [];
                    end_point_err{cond_idx}  = [];
                end
            end

            tag4_vigor_subj{ss}           = vigor;
            tag4_sac_RT_subj{ss}          = sac_RT;
            tag4_max_vel_subj{ss}         = max_vel;
            tag4_end_point_error_subj{ss} = end_point_err;
            tag4_eye_vm_vis_subj{ss}      = eye_vm_vis_out;
            tag4_eye_vm_sac_subj{ss}      = eye_vm_sac_out;

        end % session loop

        % Put subject-level into top-level
        tag1_vigor{ii}           = tag1_vigor_subj;
        tag1_sac_RT{ii}          = tag1_sac_RT_subj;
        tag1_max_vel{ii}         = tag1_max_vel_subj;
        tag1_end_point_error{ii} = tag1_end_point_error_subj;
        tag1_eye_vm_vis{ii}      = tag1_eye_vm_vis_subj;
        tag1_eye_vm_sac{ii}      = tag1_eye_vm_sac_subj;

        tag4_vigor{ii}           = tag4_vigor_subj;
        tag4_sac_RT{ii}          = tag4_sac_RT_subj;
        tag4_max_vel{ii}         = tag4_max_vel_subj;
        tag4_end_point_error{ii} = tag4_end_point_error_subj;
        tag4_eye_vm_vis{ii}      = tag4_eye_vm_vis_subj;
        tag4_eye_vm_sac{ii}      = tag4_eye_vm_sac_subj;

    end % subject loop

    %% -------- Save everything to file --------
    fprintf('\nSaving Data..........\n');
    save(fullfile(data_path, 'population_data','sac_kinematic_prim_sac_session.mat' ), ...
        'tag1_end_point_error', 'tag1_vigor', 'tag1_sac_RT', 'tag1_max_vel', ...
        'tag1_eye_vm_vis','tag1_eye_vm_sac','-v7.3');

    save(fullfile(data_path, 'population_data', 'sac_kinematic_corr_sac_session.mat'), ...
        'tag4_end_point_error', 'tag4_vigor', 'tag4_sac_RT', 'tag4_max_vel', ...
        'tag4_eye_vm_vis','tag4_eye_vm_sac','-v7.3');

    fprintf('\nDone.\n');
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
        rot_goal_x(ii)        = rot_goal(1);
        rot_goal_y(ii)        = rot_goal(2);
    end

    px_error = rot_eye_px_offset - rot_goal_x;
    py_error = rot_eye_py_offset - rot_goal_y;
end
