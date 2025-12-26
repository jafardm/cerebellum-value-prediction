function JDM_build_sac_kinematics(data_path)
    JDM_params_funcs
    fprintf('\nLoading Data..........\n');
    load(fullfile(data_path,'population_data','behave_data'));

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

    % Each variable will be a cell of size {num_subjects, 1}, where each element is {n_conds, 1}
    tag1_vigor = cell(num_subjects, 1);
    tag1_sac_RT = cell(num_subjects, 1);
    tag1_max_vel = cell(num_subjects, 1);
    tag1_end_point_error = cell(num_subjects, 1);
    tag1_eye_vm_vis = cell(num_subjects, 1);
    tag1_eye_vm_sac  = cell(num_subjects, 1);

    tag4_vigor = cell(num_subjects, 1);
    tag4_sac_RT = cell(num_subjects, 1);
    tag4_max_vel = cell(num_subjects, 1);
    tag4_end_point_error = cell(num_subjects, 1);
    tag4_eye_vm_vis = cell(num_subjects, 1);
    tag4_eye_vm_sac  = cell(num_subjects, 1);

    for ii = 1:num_subjects
        fprintf('Processing subject %d of %d...\n', ii, num_subjects);
        tot_sac_amp = data_behav.sac_amp{ii}; 
        tot_sac_vel = data_behav.sac_vel{ii}; 
        tot_sac_tag = data_behav.sac_tag{ii}; 
        tot_tsk_cond = data_behav.task_cond{ii};
        tot_tgt_cond = data_behav.tgt_cond{ii};
        tot_jmp_cond = data_behav.jump_cond{ii};
        tot_rwd_cond = data_behav.rew_cond{ii};
        tot_rt       = data_behav.sac_rt{ii};
        tot_choice   = data_behav.choice{ii};
        tot_sac_vec  = data_behav.sac_vec{ii};
        tot_gl_vec   = data_behav.gl_vec{ii};
        tot_eye_vm_vis   = data_behav.eye_vm_vis{ii};
        tot_eye_vm_sac   = data_behav.eye_vm_sac{ii};

        idx_rxt = tot_rt > 0 & tot_rt < 600;
        targ_tags = ismember(tot_sac_tag, [1 4 6 7]);

        fit_obj = fitlm(log10(tot_sac_amp(targ_tags)), log10(tot_sac_vel(targ_tags)));
        %% ---------------- Tag 1 (Primary) ----------------
        vigor = cell(n_conds_tag1, 1);
        sac_RT = cell(n_conds_tag1, 1);
        max_vel = cell(n_conds_tag1, 1);
        end_point_error = cell(n_conds_tag1, 1);
        eye_vm_vis = cell(n_conds_tag1, 1);
        eye_vm_sac = cell(n_conds_tag1, 1);

        for cond_idx = 1:n_conds_tag1
             fprintf('  Tag1 Condition %d of %d\n', cond_idx, n_conds_tag1);
            if cond_idx <= n_forc_cond
                istuned = ismember([tot_tsk_cond tot_sac_tag tot_tgt_cond tot_rwd_cond double(idx_rxt)], ...
                                   tag1_conds(cond_idx,:), 'rows');
            else
                istuned = ismember([tot_tsk_cond tot_choice tot_sac_tag double(idx_rxt)], ...
                                   choice_conds(cond_idx-n_forc_cond,:), 'rows');
            end

            amps_temp = log10(tot_sac_amp(istuned));
            vels_temp = log10(tot_sac_vel(istuned));
            
            % Wrap the predictor in a table with the correct variable name
            log_pred_vel = predict(fit_obj, table(amps_temp, 'VariableNames', {'x1'}));
        
            vigor{cond_idx} = vels_temp ./ log_pred_vel;
            sac_RT{cond_idx} = tot_rt(istuned);
            max_vel{cond_idx} = tot_sac_vel(istuned);
            eye_vm_vis{cond_idx} = tot_eye_vm_vis(:,istuned);
            eye_vm_sac{cond_idx} = tot_eye_vm_sac(:,istuned);
            gl_vec = tot_gl_vec(istuned);
            sac_vec = tot_sac_vec(istuned);
            [px_error, py_error] = compute_enpoint_error(gl_vec, sac_vec);
            end_point_error{cond_idx} = [px_error' py_error'];
           
        end

        % Store in top-level cell
        tag1_vigor{ii} = vigor;
        tag1_sac_RT{ii} = sac_RT;
        tag1_max_vel{ii} = max_vel;
        tag1_end_point_error{ii} = end_point_error;
        tag1_eye_vm_vis{ii} = eye_vm_vis;
        tag1_eye_vm_sac{ii}  = eye_vm_sac;
        
        %% ---------------- Tag 4 (Corrective) ----------------
        vigor = cell(n_conds_tag4, 1);
        sac_RT = cell(n_conds_tag4, 1);
        max_vel = cell(n_conds_tag4, 1);
        end_point_error = cell(n_conds_tag4, 1);
        eye_vm_vis = cell(n_conds_tag4, 1);
        eye_vm_sac = cell(n_conds_tag4, 1);

        for cond_idx = 1:n_conds_tag4
            fprintf('  Tag4 Condition %d of %d\n', cond_idx, n_conds_tag4);
             istuned = ismember([tot_tsk_cond tot_sac_tag tot_tgt_cond...
                       tot_jmp_cond tot_rwd_cond double(idx_rxt)],...
                       tag4_conds(cond_idx,:), 'rows');

            amps_temp = log10(tot_sac_amp(istuned));
            vels_temp = log10(tot_sac_vel(istuned));
            log_pred_vel = predict(fit_obj, amps_temp);

            vigor{cond_idx} = vels_temp ./ log_pred_vel;
            sac_RT{cond_idx} = tot_rt(istuned);
            max_vel{cond_idx} = tot_sac_vel(istuned);
            eye_vm_vis{cond_idx} = tot_eye_vm_vis(:,istuned);
            eye_vm_sac{cond_idx} = tot_eye_vm_sac(:,istuned);

            gl_vec = tot_gl_vec(istuned);
            sac_vec = tot_sac_vec(istuned);
            [px_error, py_error] = compute_enpoint_error(gl_vec, sac_vec);
            end_point_error{cond_idx} = [px_error' py_error'];
        end

        % Store in top-level cell
        tag4_vigor{ii} = vigor;
        tag4_sac_RT{ii} = sac_RT;
        tag4_max_vel{ii} = max_vel;
        tag4_end_point_error{ii} = end_point_error;
        tag4_eye_vm_vis{ii} = eye_vm_vis;
        tag4_eye_vm_sac{ii}  = eye_vm_sac;

    
    end
            %% -------- Save everything to file --------
        fprintf('\nSaving Data..........\n');
        save(fullfile(data_path, 'population_data','sac_kinematic_prim_sac.mat' ), ...
            'tag1_end_point_error', 'tag1_vigor', 'tag1_sac_RT', 'tag1_max_vel',...
            'tag1_eye_vm_vis','tag1_eye_vm_sac','-v7.3');
        
        save(fullfile(data_path, 'population_data', 'sac_kinematic_corr_sac.mat'), ...
            'tag4_end_point_error', 'tag4_vigor', 'tag4_sac_RT', 'tag4_max_vel',...
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
        rot_goal_x(ii) = rot_goal(1);
        rot_goal_y(ii) = rot_goal(2);
    end

    px_error = rot_eye_px_offset - rot_goal_x;
    py_error = rot_eye_py_offset - rot_goal_y;
end
