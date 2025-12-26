function [sac_data, meta_data] = JDM_Sac_Sorter(trials_data, meta_data,UNEYE_DATA,CONVERTED_DATA,event_data)
%% Define saccade tags
meta_data.sac_tag_list = { ...
    'prim_success', ... % tag 1
    'prim_attempt', ... % tag 2
    'prim_fail', ... % tag 3
    'corr_success', ... % tag 4
    'corr_fail', ... % tag 5
    'back_center_success', ... % tag 6
    'back_center_prim', ... % tag 7
    'back_center_irrelev', ... % tag 8
    'target_irrelev', ... % tag 9
    'other_irrelev', ... % tag 10
    'prim_no_corr',... % tag 11; prim. sac. that is not followed by corr. sac.
    'db_corr_success',... % tag 12; sac. that follows first corr. sac., back to 2nd jumped cue
    'corr_no_db_corr',... % tag 13; corr. sac. that is not followed by another corr. sac.
    'other_irrelev_visual',... % tag 14; like tag 10, but visual ang. based on the pursuit target present on the screen
    'back_center_irrelev_visual',... % tag 15; like tag 10, but visual start position based on offset of previous saccade
    ... % Add more tags here, do not reorder or change the tags defined above.
    };

%% Analyze trials
num_trials = length(trials_data.time_start);
counter_valid_trial = 1; % keeps track of number of valid trials
sampling_freq   = 1e3;
for counter_trial = 1 : num_trials
    %% Clear variables
    clearvars -except meta_data trials_data counter_trial num_trials SACS_ALL counter_valid_trial trial_num ...
                      sampling_freq UNEYE_DATA CONVERTED_DATA event_data
    
    %% Extract Trial Varibales
    TRIAL.start_x       = trials_data.start_x(1,counter_trial);
    TRIAL.start_y       = trials_data.start_y(1,counter_trial);
    
    TRIAL.cue_x         = trials_data.cue_x(1,counter_trial);
    TRIAL.cue_y         = trials_data.cue_y(1,counter_trial);
    TRIAL.cue_x_high_rew   = trials_data.cue_x_high_rew(  1,counter_trial);     % # only applicable under choice condition
    TRIAL.cue_x_low_rew    = trials_data.cue_x_low_rew(  1,counter_trial);      % # only applicable under choice condition
    TRIAL.cue_y_high_rew   = trials_data.cue_y_high_rew(  1,counter_trial);     % # only applicable under choice condition
    TRIAL.cue_y_low_rew    = trials_data.cue_y_low_rew(  1,counter_trial);      % # only applicable under choice condition

    TRIAL.end_x         = trials_data.end_x(  1,counter_trial);
    TRIAL.end_y         = trials_data.end_y(  1,counter_trial);
    TRIAL.iss_x         = trials_data.iss_x(  1,counter_trial);
    TRIAL.iss_y         = trials_data.iss_y(  1,counter_trial);

    TRIAL.tgt_num          = trials_data.tgt_num(  1,counter_trial);            % # tgt. num, under forced condition
    TRIAL.high_rew_tgt_num = trials_data.high_rew_tgt_num(  1,counter_trial);   % # tgt. num. for high reward, if under high reward condition
    TRIAL.low_rew_tgt_num  = trials_data.low_rew_tgt_num(  1,counter_trial);    % # tgt. num for low reward, if under low reward condition
    TRIAL.task_cond        = trials_data.task_cond(  1,counter_trial);
    TRIAL.choice           = trials_data.choice(  1,counter_trial);
    TRIAL.rew_cond         = trials_data.rew_cond(  1,counter_trial);
    TRIAL.jump_cond        = trials_data.jump_cond(  1,counter_trial);
    TRIAL.tgt_cond         = trials_data.tgt_cond(  1,counter_trial);

    TRIAL.time_start    = trials_data.time_start(1,counter_trial);
    TRIAL.time_end      = trials_data.time_end(1,counter_trial);
    variable_list = {'time_state_str_pursuit',...
        'time_state_str_present',...
        'time_state_str_fixation',...
        'time_state_cue_present',...
        'time_state_sac_detect_on',...
        'time_state_sac_onset',...
        'time_state_sac_detect_off',...
        'time_state_reward',...
        'time_state_end_fixation',...
         };

    for counter_variable = 1 : length(variable_list)
        variable_name = variable_list{counter_variable};
        if isa(trials_data.(variable_name),'double')
            TRIAL.(variable_name) = trials_data.(variable_name)(1,counter_trial);
        elseif isa(trials_data.(variable_name),'cell')
            TRIAL.(variable_name) = trials_data.(variable_name){1,counter_trial};
        end
    end
    TRIAL.time_iti        = trials_data.time_iti(1,counter_trial);
    TRIAL.time_punishment = trials_data.time_punishment(1,counter_trial);
    TRIAL.time_fixation   = trials_data.time_fixation(1,counter_trial);
    TRIAL.time_pursuit    = trials_data.time_pursuit(1,counter_trial);
    TRIAL.time_pursuit(isnan(TRIAL.time_pursuit)) = 0.200;
    
    time_device = trials_data.time_1K{1,counter_trial};

    TRIAL.tgt_px        = trials_data.tgt_px{       1,counter_trial};
    TRIAL.tgt_py        = trials_data.tgt_py{       1,counter_trial};
    
    TRIAL.eye_px_filt = trials_data.eye_px_filt{1,counter_trial};
    TRIAL.eye_py_filt = trials_data.eye_py_filt{1,counter_trial};
    TRIAL.eye_vx_filt = trials_data.eye_vx_filt{1,counter_trial};
    TRIAL.eye_vy_filt = trials_data.eye_vy_filt{1,counter_trial};
    TRIAL.eye_vm_filt = trials_data.eye_vm_filt{1,counter_trial};
    
    % Backstep exp. parameters
    if isfield(trials_data, 'time_state_dwell_on')
        TRIAL.time_state_dwell_on = trials_data.time_state_dwell_on{1,counter_trial};
        TRIAL.time_dwell = trials_data.time_dwell(1,counter_trial);
        % If dwell time is too long, the trial might have ended w/o tgt. jumping back to cue, so 
        % time state dwell off may not exist
        if isfield(trials_data, 'time_state_dwell_off')
            TRIAL.time_state_dwell_off = trials_data.time_state_dwell_off{1,counter_trial};
        end
    end
    %% Extract Saccades
    clearvars -except meta_data trials_data counter_trial num_trials TRIAL SACS_ALL counter_valid_trial trial_num time_device...
                      sampling_freq r_eye_tracked l_eye_tracked UNEYE_DATA CONVERTED_DATA event_data

    all_sac_ind_onset_temp  = UNEYE_DATA.Sac.(sprintf('trial_%d',counter_trial)).all_sac_ind_onset;
    all_sac_ind_vmax_temp   = UNEYE_DATA.Sac.(sprintf('trial_%d',counter_trial)).all_sac_ind_vmax;
    all_sac_ind_offset_temp = UNEYE_DATA.Sac.(sprintf('trial_%d',counter_trial)).all_sac_ind_offset;
    all_sac_validity_temp   = UNEYE_DATA.Sac.(sprintf('trial_%d',counter_trial)).all_sac_validity;
    % Skip the trial if no saccade valid
    if isempty(all_sac_validity_temp)
        continue
    end
    all_sac_ind_onset_temp = reshape(all_sac_ind_onset_temp, 1, length(all_sac_ind_onset_temp));
    all_sac_ind_vmax_temp = reshape(all_sac_ind_vmax_temp, 1, length(all_sac_ind_vmax_temp));
    all_sac_ind_offset_temp = reshape(all_sac_ind_offset_temp, 1, length(all_sac_ind_offset_temp));
    all_sac_validity_temp = reshape(all_sac_validity_temp, 1, length(all_sac_validity_temp));
    max_ind = length(CONVERTED_DATA.(sprintf('trial_%d',counter_trial)).time);
    all_sac_ind_onset_temp(all_sac_ind_onset_temp<1) = 1;
    all_sac_ind_vmax_temp(all_sac_ind_vmax_temp<1) = 1;
    all_sac_ind_offset_temp(all_sac_ind_offset_temp<1) = 1;
    all_sac_ind_onset_temp(all_sac_ind_onset_temp>max_ind) = max_ind;
    all_sac_ind_vmax_temp(all_sac_ind_vmax_temp>max_ind) = max_ind;
    all_sac_ind_offset_temp(all_sac_ind_offset_temp>max_ind) = max_ind;
    all_sac_time_onset = CONVERTED_DATA.(sprintf('trial_%d',counter_trial)).time(all_sac_ind_onset_temp);
    all_sac_time_vmax = CONVERTED_DATA.(sprintf('trial_%d',counter_trial)).time(all_sac_ind_vmax_temp);
    all_sac_time_offset = CONVERTED_DATA.(sprintf('trial_%d',counter_trial)).time(all_sac_ind_offset_temp);
    temp = repmat(time_device(:), [1 length(all_sac_time_onset)]);
    [~, all_sac_ind_onset] = min(abs(temp-all_sac_time_onset));
    [~, all_sac_ind_vmax] = min(abs(temp-all_sac_time_vmax));
    [~, all_sac_ind_offset] = min(abs(temp-all_sac_time_offset));
    all_sac_validity  = all_sac_validity_temp;
    
    %% Init SACS_ALL and add common params
    SACS_ALL_TRIAL = struct;
    SACS_ALL_TRIAL.validity    = reshape(all_sac_validity,1, []);
    SACS_ALL_TRIAL.trial_num   = reshape(ones(size(all_sac_validity))*counter_trial, 1, []);
    SACS_ALL_TRIAL.tag         = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.count       = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.flag_last_cue = false(size(all_sac_validity));
    SACS_ALL_TRIAL.time_onset  = reshape(time_device(all_sac_ind_onset), 1, []);
    SACS_ALL_TRIAL.time_vmax   = reshape(time_device(all_sac_ind_vmax),  1, []);
    SACS_ALL_TRIAL.time_offset = reshape(time_device(all_sac_ind_offset),1, []);
    SACS_ALL_TRIAL.start_x = ones(size(all_sac_ind_onset))*TRIAL.start_x;
    SACS_ALL_TRIAL.start_y = ones(size(all_sac_ind_onset))*TRIAL.start_y;
    SACS_ALL_TRIAL.cue_x = ones(size(all_sac_ind_onset))*TRIAL.cue_x;
    SACS_ALL_TRIAL.cue_y = ones(size(all_sac_ind_onset))*TRIAL.cue_y;
    SACS_ALL_TRIAL.end_x = ones(size(all_sac_ind_onset))*TRIAL.end_x;
    SACS_ALL_TRIAL.end_y = ones(size(all_sac_ind_onset))*TRIAL.end_y;
    
    % REWARD TASK PARAMS
    SACS_ALL_TRIAL.tgt_num          = ones(size(all_sac_ind_onset))*TRIAL.tgt_num;          
    SACS_ALL_TRIAL.cue_x_high_rew   = ones(size(all_sac_ind_onset))*TRIAL.cue_x_high_rew;     
    SACS_ALL_TRIAL.cue_x_low_rew    = ones(size(all_sac_ind_onset))*TRIAL.cue_x_low_rew;       
    SACS_ALL_TRIAL.cue_y_high_rew   = ones(size(all_sac_ind_onset))*TRIAL.cue_y_high_rew;     
    SACS_ALL_TRIAL.cue_y_low_rew    = ones(size(all_sac_ind_onset))*TRIAL.cue_y_low_rew;       
    SACS_ALL_TRIAL.high_rew_tgt_num = ones(size(all_sac_ind_onset))*TRIAL.high_rew_tgt_num;  
    SACS_ALL_TRIAL.low_rew_tgt_num  = ones(size(all_sac_ind_onset))*TRIAL.low_rew_tgt_num;   
    SACS_ALL_TRIAL.task_cond       = ones(size(all_sac_ind_onset))*TRIAL.task_cond;
    SACS_ALL_TRIAL.choice          = ones(size(all_sac_ind_onset))*TRIAL.choice;
    SACS_ALL_TRIAL.rew_cond        = ones(size(all_sac_ind_onset))*TRIAL.rew_cond;
    SACS_ALL_TRIAL.jump_cond        = ones(size(all_sac_ind_onset))*TRIAL.jump_cond;
    SACS_ALL_TRIAL.tgt_cond        = ones(size(all_sac_ind_onset))*TRIAL.tgt_cond;

    SACS_ALL_TRIAL.time_visual      = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.time_auditory    = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_px_onset  = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_px_offset = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_py_onset  = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_py_offset = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.reaction         = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_amp_x     = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_amp_y     = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_amp_m     = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.visual_ang       = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.diff_start       = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.diff_finish      = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.diff_ang         = nan(size(all_sac_validity));


    SACS_ALL_TRIAL.duration         = reshape((SACS_ALL_TRIAL.time_offset - SACS_ALL_TRIAL.time_onset) * 1000.0, 1, []);

    SACS_ALL_TRIAL.eye_vm_max     = reshape(TRIAL.eye_vm_filt(all_sac_ind_vmax),   1, []);
    SACS_ALL_TRIAL.eye_px_onset   = reshape(TRIAL.eye_px_filt(all_sac_ind_onset),  1, []);
    SACS_ALL_TRIAL.eye_px_offset  = reshape(TRIAL.eye_px_filt(all_sac_ind_offset), 1, []);
    SACS_ALL_TRIAL.eye_py_onset   = reshape(TRIAL.eye_py_filt(all_sac_ind_onset),  1, []);
    SACS_ALL_TRIAL.eye_py_offset  = reshape(TRIAL.eye_py_filt(all_sac_ind_offset), 1, []);
    SACS_ALL_TRIAL.eye_amp_x      = reshape(SACS_ALL_TRIAL.eye_px_offset - SACS_ALL_TRIAL.eye_px_onset, 1, []);
    SACS_ALL_TRIAL.eye_amp_y      = reshape(SACS_ALL_TRIAL.eye_py_offset - SACS_ALL_TRIAL.eye_py_onset, 1, []);
    SACS_ALL_TRIAL.eye_amp_m      = reshape(sqrt(SACS_ALL_TRIAL.eye_amp_x.^2 + SACS_ALL_TRIAL.eye_amp_y.^2), 1, []);
    SACS_ALL_TRIAL.eye_ang        = reshape(atan2d(SACS_ALL_TRIAL.eye_amp_y, SACS_ALL_TRIAL.eye_amp_x), 1, []);
    SACS_ALL_TRIAL.tgt_px_onset   = reshape(TRIAL.tgt_px(all_sac_ind_onset),  1, []);
    SACS_ALL_TRIAL.tgt_px_offset  = reshape(TRIAL.tgt_px(all_sac_ind_offset), 1, []);
    SACS_ALL_TRIAL.tgt_py_onset   = reshape(TRIAL.tgt_py(all_sac_ind_onset),  1, []);
    SACS_ALL_TRIAL.tgt_py_offset  = reshape(TRIAL.tgt_py(all_sac_ind_offset), 1, []);
    % Backstep exp. parameters
    if isfield(TRIAL, 'time_state_dwell_on')
        SACS_ALL_TRIAL.time_state_dwell_on = ones(size(all_sac_ind_onset))*TRIAL.time_state_dwell_on(end);
        SACS_ALL_TRIAL.time_dwell = ones(size(all_sac_ind_onset))*TRIAL.time_dwell; 
        % If dwell time is too long, the trial might have ended w/o tgt. jumping back to cue, so 
        % time state dwell off may not exist
        if isfield(TRIAL, 'time_state_dwell_off')
            SACS_ALL_TRIAL.time_state_dwell_off = ones(size(all_sac_ind_onset))*TRIAL.time_state_dwell_off(end);
        else
            SACS_ALL_TRIAL.time_state_dwell_off = ones(size(all_sac_ind_onset))*nan;
        end
    end
    %% Tag saccades
    threshold_pos = 2.2; % default 1.5deg; for backstep 1.7deg
    threshold_ang = 45.0; % deg
    threshold_micro_saccade = 1.0; % deg
    
    % Tag the primary saccades
    length_cue_presentation = length(TRIAL.time_state_cue_present);
    for counter_cue_pres = 1 : length_cue_presentation
        time_start_search  = TRIAL.time_state_cue_present(counter_cue_pres);
        % For primSac. that comes after last cue present., search the period after the last cue_present till the next_trial
        if counter_cue_pres == length_cue_presentation
            time_finish_search = TRIAL.time_end;
        % For the rest, search the period after the cue_present. till the next str_present
        else
            time_finish_search = TRIAL.time_state_str_fixation( find(TRIAL.time_state_str_fixation > time_start_search, 1, 'first') );
        end
        idx_sac = find((SACS_ALL_TRIAL.time_onset > time_start_search) & (SACS_ALL_TRIAL.time_onset < time_finish_search), 1, 'first');
        idx_sac_end = find((SACS_ALL_TRIAL.time_onset > time_start_search) & (SACS_ALL_TRIAL.time_onset < time_finish_search), 1, 'last');
        if ~isempty(idx_sac) && ~isempty(idx_sac_end)
            while (idx_sac <= idx_sac_end) && (SACS_ALL_TRIAL.eye_amp_m( idx_sac) < threshold_micro_saccade)
                idx_sac = idx_sac + 1;
            end
            if idx_sac == idx_sac_end + 1
                continue;
            end
            SACS_ALL_TRIAL.time_visual(idx_sac)      = time_start_search;
            SACS_ALL_TRIAL.time_auditory(idx_sac)    = time_start_search; % neutral beep at the cue presentation
            SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = TRIAL.start_x;
            SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = TRIAL.start_y;
            SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.cue_x;
            SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.cue_y;

            visual_low_px_offset  = TRIAL.cue_x_low_rew;
            visual_low_py_offset  = TRIAL.cue_y_low_rew;
            visual_high_px_offset = TRIAL.cue_x_high_rew;
            visual_high_py_offset = TRIAL.cue_y_high_rew;

            SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
            SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
            SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
            SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
            SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
            visual_amp_x_ = SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.eye_px_onset(idx_sac);
            visual_amp_y_ = SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.eye_py_onset(idx_sac);
            visual_amp_m_ = sqrt(visual_amp_x_.^2 + visual_amp_y_.^2);

            visual_amp_x_low  = visual_low_px_offset - SACS_ALL_TRIAL.eye_px_onset(idx_sac);
            visual_amp_y_low  = visual_low_py_offset - SACS_ALL_TRIAL.eye_py_onset(idx_sac);
            visual_amp_m_low  = sqrt(visual_amp_x_low.^2 + visual_amp_y_low.^2);
            visual_amp_x_high = visual_high_px_offset - SACS_ALL_TRIAL.eye_px_onset(idx_sac);
            visual_amp_y_high = visual_high_py_offset - SACS_ALL_TRIAL.eye_py_onset(idx_sac);
            visual_amp_m_high = sqrt(visual_amp_x_high.^2 + visual_amp_y_high.^2);

            eye_amp_x_ = SACS_ALL_TRIAL.eye_amp_x(idx_sac);
            eye_amp_y_ = SACS_ALL_TRIAL.eye_amp_y(idx_sac);
            eye_amp_m_ = SACS_ALL_TRIAL.eye_amp_m(idx_sac);
            SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
                ((SACS_ALL_TRIAL.eye_px_onset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac)).^2) + ...
                ((SACS_ALL_TRIAL.eye_py_onset( idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac)).^2));
            SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2));
            
            diff_finish_low  = sqrt( ...
                ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - visual_low_px_offset).^2) + ...
                ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - visual_low_py_offset).^2));
            diff_finish_high = sqrt( ...
                ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - visual_high_px_offset).^2) + ...
                ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - visual_high_py_offset).^2));
            
            diff_finish_tgt_end = sqrt( ...
                ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - TRIAL.end_x).^2) + ...
                ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - TRIAL.end_y).^2) );
            SACS_ALL_TRIAL.diff_ang(idx_sac) = abs(acosd(... % x dot y = |x||y|cos(diff_ang); angle btwn vector from eye start pos. to cue and vector of saccade trajectory
                (visual_amp_x_.* eye_amp_x_ + visual_amp_y_.* eye_amp_y_)./visual_amp_m_./eye_amp_m_) );

            diff_ang_low  = abs(acosd(... % x dot y = |x||y|cos(diff_ang); angle btwn vector from eye start pos. to cue and vector of saccade trajectory
                (visual_amp_x_low.* eye_amp_x_ + visual_amp_y_low.* eye_amp_y_)./visual_amp_m_low./eye_amp_m_) );
            diff_ang_high = abs(acosd(... % x dot y = |x||y|cos(diff_ang); angle btwn vector from eye start pos. to cue and vector of saccade trajectory
                (visual_amp_x_high.* eye_amp_x_ + visual_amp_y_high.* eye_amp_y_)./visual_amp_m_high./eye_amp_m_) );

            if counter_cue_pres == length_cue_presentation
                flag_last_cue = true;
                SACS_ALL_TRIAL.flag_last_cue(idx_sac) = flag_last_cue;
            end
            
            % If saccade starts near start tgt. && moves twrd. cue tgt. && lands near cue or end tgt., then 'prim_success' tag 1
            % Note: allowing primSac. landing near end tgt. bc saccade gets hypometric as exp. progresses (only for last cue presentation),
            % so saccade in trial with end tgt. that steps back from cue tgt. is counted
            % Note: saccade could be considered successful primSac. and even have successful corrSac. that follows but not comes
            % after the last cue present. in a trial, bc. the calibration during exp. may be off and thus counting "successful" saccade 
            % as incorrect and repeating cue present.
            
            % Forced
            if TRIAL.task_cond == 1
                if ((SACS_ALL_TRIAL.diff_start( idx_sac) < threshold_pos) && (SACS_ALL_TRIAL.diff_ang(idx_sac) < threshold_ang) &&...
                        ((SACS_ALL_TRIAL.diff_finish(idx_sac) < threshold_pos) || ((diff_finish_tgt_end < threshold_pos)&&(counter_cue_pres == length_cue_presentation))))
                    SACS_ALL_TRIAL.tag(idx_sac) = 1; 
                % If saccade starts near start tgt. && moves twrd. cue tgt. but not lands near cue nor end tgt., then 'prim_attempt' tag 2
                elseif ((SACS_ALL_TRIAL.diff_start( idx_sac) < threshold_pos) &&(SACS_ALL_TRIAL.diff_ang(idx_sac) < threshold_ang ) &&...
                        ~((SACS_ALL_TRIAL.diff_finish(idx_sac) < threshold_pos) || (diff_finish_tgt_end < threshold_pos)))
                    SACS_ALL_TRIAL.tag(idx_sac) = 2; 
                % If saccade starts near start tgt. but not moves twrd. cue tgt. (implying not landing near cue nor end tgt.), then 'prim_fail' tag 3
                elseif ((SACS_ALL_TRIAL.diff_start( idx_sac) < threshold_pos) && (SACS_ALL_TRIAL.diff_ang(idx_sac) >= threshold_ang))
                    SACS_ALL_TRIAL.tag(idx_sac) = 3; 
                end
            % Choice
            elseif TRIAL.task_cond == 0
                if SACS_ALL_TRIAL.diff_start( idx_sac) < threshold_pos
                    % Chose low reward
                    if (diff_ang_low <= threshold_ang) && (diff_finish_low <= threshold_pos)
                        SACS_ALL_TRIAL.visual_px_offset(idx_sac) = visual_low_px_offset;
                        SACS_ALL_TRIAL.visual_py_offset(idx_sac) = visual_low_py_offset;
                        SACS_ALL_TRIAL.tag(idx_sac) = 1;
                    % Chose high reward
                    elseif (diff_ang_high <= threshold_ang) && (diff_finish_high <= threshold_pos)
                        SACS_ALL_TRIAL.visual_px_offset(idx_sac) = visual_high_px_offset;
                        SACS_ALL_TRIAL.visual_py_offset(idx_sac) = visual_high_py_offset;
                        SACS_ALL_TRIAL.tag(idx_sac) = 1;
                    % Attempted low reward
                    elseif (diff_ang_low <= threshold_ang) && (diff_finish_low > threshold_pos)
                        SACS_ALL_TRIAL.visual_px_offset(idx_sac) = visual_low_px_offset;
                        SACS_ALL_TRIAL.visual_py_offset(idx_sac) = visual_low_py_offset;
                        SACS_ALL_TRIAL.tag(idx_sac) = 2;
                    % Attempted high reward
                    elseif (diff_ang_high <= threshold_ang) && (diff_finish_high > threshold_pos)
                        SACS_ALL_TRIAL.visual_px_offset(idx_sac) = visual_high_px_offset;
                        SACS_ALL_TRIAL.visual_py_offset(idx_sac) = visual_high_py_offset;
                        SACS_ALL_TRIAL.tag(idx_sac) = 2;
                    % Failed saccade
                    elseif (diff_ang_low > threshold_ang) && (diff_ang_high > threshold_ang) 
                        SACS_ALL_TRIAL.tag(idx_sac) = 3;
                    end
                    
                
                end
                SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
                SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
                visual_amp_x_ = SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.eye_px_onset(idx_sac);
                visual_amp_y_ = SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.eye_py_onset(idx_sac);
                visual_amp_m_ = sqrt(visual_amp_x_.^2 + visual_amp_y_.^2);
                SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2));
                
                SACS_ALL_TRIAL.diff_ang(idx_sac) = abs(acosd(... % x dot y = |x||y|cos(diff_ang); angle btwn vector from eye start pos. to cue and vector of saccade trajectory
                    (visual_amp_x_.* eye_amp_x_ + visual_amp_y_.* eye_amp_y_)./visual_amp_m_./eye_amp_m_) );
            end
            
        end
    end
    
    % Applies to jump and forced condition only
    % Tag the corr_success and corr_fail
    % search the period after a success sac_prim offset till the next str_present
    % or till the next_trial
    idx_prim_success = find(SACS_ALL_TRIAL.tag == 1); % find prim_success
    if ~isempty(idx_prim_success) && (TRIAL.jump_cond == 1) && (TRIAL.task_cond == 1)
        for counter_prim = 1 : length(idx_prim_success)
            idx_prim = idx_prim_success(counter_prim);
            time_last_cue_pres = TRIAL.time_state_cue_present(end);
            time_start_search = SACS_ALL_TRIAL.time_offset(idx_prim);
            if time_start_search < time_last_cue_pres % refer to 'prim_success' tag 1 comments to see how 'corr_success' could occur but exp. doesn't proceed to next trial
                % search from prim sac offset till the next str presentation
                time_finish_search = TRIAL.time_state_str_fixation( find(TRIAL.time_state_str_fixation > time_start_search, 1, 'first') );
                flag_last_cue = false;
            else
                % search from prim sac offset till end of trial
                time_finish_search = TRIAL.time_end;
                flag_last_cue = true;
            end
            idx_sac = find((SACS_ALL_TRIAL.time_onset > time_start_search) & (SACS_ALL_TRIAL.time_onset < time_finish_search), 1, 'first');
            if ~isempty(idx_sac)
                SACS_ALL_TRIAL.time_visual(idx_sac)      = time_start_search;
                SACS_ALL_TRIAL.time_auditory(idx_sac)    = time_start_search; % good beep at the end of prim success
                % we will use tgt_cue and tgt_end for visual params, but we will use
                % the prim offset for diff values
                SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = TRIAL.cue_x;
                SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = TRIAL.cue_y;
                SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.end_x;
                SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.end_y;
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
                SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
                SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
                % we will use the prim offset for diff values
                visual_amp_x_prim = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.eye_px_offset(idx_prim));
                visual_amp_y_prim = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.eye_py_offset(idx_prim));
                visual_amp_m_prim = sqrt((visual_amp_x_prim.^2) + (visual_amp_y_prim.^2));
                SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_onset(idx_sac) - SACS_ALL_TRIAL.eye_px_offset(idx_prim)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_onset(idx_sac) - SACS_ALL_TRIAL.eye_py_offset(idx_prim)).^2) );
                SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_ang(idx_sac) = abs(acosd( ...
                    ( (SACS_ALL_TRIAL.eye_amp_x(idx_sac) .* visual_amp_x_prim) + ...
                      (SACS_ALL_TRIAL.eye_amp_y(idx_sac) .* visual_amp_y_prim) ) ...
                    ./ (SACS_ALL_TRIAL.eye_amp_m(idx_sac)) ./ (visual_amp_m_prim) ));
           
                % If end tgt. is different from cue tgt. && saccade starts near where eye was @ previous successful primSac. offset
                % && moves twrd. and lands near end tgt., then 'corr_success' tag 4
                if ((SACS_ALL_TRIAL.visual_amp_m(idx_sac) > eps)&&(SACS_ALL_TRIAL.diff_start(idx_sac) < threshold_pos)...
                        &&(SACS_ALL_TRIAL.diff_ang(idx_sac) < threshold_ang)&&(SACS_ALL_TRIAL.diff_finish(idx_sac) < threshold_pos))
                    SACS_ALL_TRIAL.tag(idx_sac) = 4; 
                    SACS_ALL_TRIAL.flag_last_cue(idx_sac) = flag_last_cue;
                % If saccade doesn't meet reqs. for 'corr_success' but still starts near where eye was @ previous successful primSac. offset,
                % then 'corr_fail' tag 5
                elseif (SACS_ALL_TRIAL.diff_start( idx_sac) < threshold_pos)
                    SACS_ALL_TRIAL.tag(idx_sac) = 5; 
                    SACS_ALL_TRIAL.flag_last_cue(idx_sac) = flag_last_cue;
                end
            end
        end
    end
    
    % Tag back_center
    % search the period before cue_present till the previous str_present (str_pursuit)
    % tag them back_center_irrelev, later we will revisit these tags.
    for counter_cue_pres = 1 : length(TRIAL.time_state_cue_present)
        time_finish_search  = TRIAL.time_state_cue_present(counter_cue_pres);
        if counter_cue_pres==1
            time_start_search = TRIAL.time_start;
        else
            time_start_search = TRIAL.time_state_str_fixation( find(TRIAL.time_state_str_fixation < time_finish_search, 1, 'last') ) - TRIAL.time_pursuit; % str_pursuit is earlier than str_present by pursuit duration (usually 200 ms)
            % Not simply using 'time_state_str_pursuit,' as finite state machine of the exp., can change the state from 'DETECT_SACCADE_START' to
            % 'INCORRECT_SACCADE' to 'STR_TARGET_PURSUIT' and we recorded 'time_state_str_pursuit' when state changed from 'INIT' to 'STR_TARGET_PURSUIT,'
            % so in the case of this sequence of states, 'time_state_str_pursuit' is not when the pursuit target, which results in the back-to-center saccade
            % of the current 'counter_cue_pres' loop, appears
        end
        idx_sac = find((SACS_ALL_TRIAL.time_onset > time_start_search) & (SACS_ALL_TRIAL.time_onset < time_finish_search), 1, 'last');
        if ~isempty(idx_sac)
            if (SACS_ALL_TRIAL.tag(idx_sac) >= 1) && (SACS_ALL_TRIAL.tag(idx_sac) <= 3)
                % 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3
                % if the saccade is tagged as a prim sac, do not re-tag it
                % and skip to next saccade.
                continue;
            end
            if idx_sac ~= 1
                SACS_ALL_TRIAL.time_visual(idx_sac) = SACS_ALL_TRIAL.time_offset(idx_sac-1);
            else
                SACS_ALL_TRIAL.time_visual(idx_sac) = TRIAL.time_state_str_pursuit(idx_sac);
            end
            SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0; % counting the rxn. time from presentation of moving pursuit target
            SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = SACS_ALL_TRIAL.eye_px_onset(idx_sac);
            SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = SACS_ALL_TRIAL.eye_py_onset(idx_sac);
            SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.start_x; % used for CS-on later
            SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.start_y; % used for CS-on later
            SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
            SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
            SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
            SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
            SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
                ((SACS_ALL_TRIAL.eye_px_onset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac)).^2) + ...
                ((SACS_ALL_TRIAL.eye_py_onset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac)).^2) );
            SACS_ALL_TRIAL.diff_ang(idx_sac) = abs(acosd( ...
                ( (SACS_ALL_TRIAL.eye_amp_x(idx_sac) .* SACS_ALL_TRIAL.visual_amp_x(idx_sac)) + ...
                  (SACS_ALL_TRIAL.eye_amp_y(idx_sac) .* SACS_ALL_TRIAL.visual_amp_y(idx_sac)) ) ...
                ./ (SACS_ALL_TRIAL.eye_amp_m(idx_sac)) ./ (SACS_ALL_TRIAL.visual_amp_m(idx_sac)) ));
            SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );

            % If saccade lands near start tgt., then 'back_center_irrelev' tag 8
            % Visual direction based on sac. onset pos. and start tgt. pos.
            if (SACS_ALL_TRIAL.diff_finish(idx_sac) < threshold_pos)
                SACS_ALL_TRIAL.tag(idx_sac) = 8; 
            % If saccade doesn't land near start tgt., then 'other_irrelev' tag 10  
            elseif (SACS_ALL_TRIAL.diff_finish(idx_sac) >= threshold_pos)
                SACS_ALL_TRIAL.tag(idx_sac) = 10; 
            end
        end
    end
    
    % Re-tag the back_center_irrelev
    % review tags and based on the events before these tags, re-tag them as
    % back_center_success or back_center_prim
    idx_back_center = find(SACS_ALL_TRIAL.tag == 8); % find back_center_irrelev
    if ~isempty(idx_back_center)
        for counter_back_center = 1 : length(idx_back_center)
            idx_sac = idx_back_center(counter_back_center);
            if (idx_sac == 1) 
                % this is the 1st sac of the trial, i.e., back_center_success
                SACS_ALL_TRIAL.time_visual(idx_sac)      = TRIAL.time_state_str_pursuit(idx_sac);
                SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = SACS_ALL_TRIAL.eye_px_onset(idx_sac);
                SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = SACS_ALL_TRIAL.eye_py_onset(idx_sac);
                SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.start_x;
                SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.start_y;
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
                SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
                SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
                SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_onset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_onset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_ang(idx_sac) = abs(acosd( ...
                    ( (SACS_ALL_TRIAL.eye_amp_x(idx_sac) .* SACS_ALL_TRIAL.visual_amp_x(idx_sac)) + ...
                      (SACS_ALL_TRIAL.eye_amp_y(idx_sac) .* SACS_ALL_TRIAL.visual_amp_y(idx_sac)) ) ...
                    ./ (SACS_ALL_TRIAL.eye_amp_m(idx_sac)) ./ (SACS_ALL_TRIAL.visual_amp_m(idx_sac))));
                SACS_ALL_TRIAL.tag(idx_sac) = 6; % 'back_center_success' tag 6
                continue;
            end
            if (SACS_ALL_TRIAL.tag(idx_sac-1) >= 1) && (SACS_ALL_TRIAL.tag(idx_sac-1) <= 5)
                % 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3 % 'corr_success' tag 4 % 'cord_fail' tag 5
                % If the previous sac is a prim or corr saccade, then this is a back_center_prim.
                % The time_visual in case of back_center_prim, can be considered as when the animal observes the target
                % after the str_pursuit which is the time_visual. However, another definition can be when the boop got played 
                % at the end of prim saccade. To cover this scenario we added time_auditory.
                SACS_ALL_TRIAL.time_auditory(idx_sac)      = SACS_ALL_TRIAL.time_offset(idx_sac-1); % bad boop at the the end of prim attempt/fail
                % In case of corrective saccades the time_visual should be set to the offset of corr sac and not str_pursuit
                if (SACS_ALL_TRIAL.tag(idx_sac-1) >= 4) && (SACS_ALL_TRIAL.tag(idx_sac-1) <= 5)
                    % 'corr_success' tag 4 % 'cord_fail' tag 5
                    SACS_ALL_TRIAL.time_visual(idx_sac)      = SACS_ALL_TRIAL.time_offset(idx_sac-1);
                end
                SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = SACS_ALL_TRIAL.eye_px_offset(idx_sac-1);
                SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = SACS_ALL_TRIAL.eye_py_offset(idx_sac-1);
                SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.start_x;
                SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.start_y;
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
                SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
                SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
                SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_onset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_onset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_ang(idx_sac) = abs(acosd( ...
                    ((SACS_ALL_TRIAL.eye_amp_x(idx_sac) .* SACS_ALL_TRIAL.visual_amp_x(idx_sac)) + ...
                    (SACS_ALL_TRIAL.eye_amp_y(idx_sac) .* SACS_ALL_TRIAL.visual_amp_y(idx_sac))) ...
                    ./ (SACS_ALL_TRIAL.eye_amp_m(idx_sac)) ./ (SACS_ALL_TRIAL.visual_amp_m(idx_sac))));
                SACS_ALL_TRIAL.tag(idx_sac) = 7; % 'back_center_prim' tag 7
                continue;
            end
        end
    end
    
    % Tag the potential back_center_prim which happen after prim or corr
    % but did not trigger the "time_state_cue_present"
    idx_prim_corr = find((SACS_ALL_TRIAL.tag >= 1) & (SACS_ALL_TRIAL.tag <= 5));
    % find 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3 % 'corr_success' tag 4 % 'corr_fail' tag 5
    if ~isempty(idx_prim_corr)
        for counter_prim_corr = 1 : length(idx_prim_corr)
            idx_sac = idx_prim_corr(counter_prim_corr)+1;
            if idx_sac > length(SACS_ALL_TRIAL.tag)
                % the idx is out of range. break the for loop.
                break;
            end
            if ~isnan(SACS_ALL_TRIAL.tag(idx_sac))
                % the saccade has been tagged before, go to the next
                % saccade
                continue;
            end
            diff_finish = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - TRIAL.start_x).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - TRIAL.start_y).^2) );
            if diff_finish < threshold_pos
                % this is a missed back_center_prim and should be tagged.
                
                % The time_visual in case of back_center_prim, can be considered as when the animal observes the target
                % after the str_pursuit which is the time_visual. However, another definition can be when the boop got played 
                % at the end of prim saccade. To cover this scenario we added time_auditory.
                SACS_ALL_TRIAL.time_visual(idx_sac)      = SACS_ALL_TRIAL.time_offset(idx_sac-1);
                SACS_ALL_TRIAL.time_auditory(idx_sac)    = SACS_ALL_TRIAL.time_offset(idx_sac-1); % bad boop at the the end of prim attempt/fail
                if (SACS_ALL_TRIAL.tag(idx_sac-1) >= 1) && (SACS_ALL_TRIAL.tag(idx_sac-1) <= 3)
                    % 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3
                    SACS_ALL_TRIAL.time_visual(idx_sac)      = SACS_ALL_TRIAL.time_offset(idx_sac-1) + TRIAL.time_punishment;
                end
                SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = SACS_ALL_TRIAL.eye_px_offset(idx_sac-1);
                SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = SACS_ALL_TRIAL.eye_py_offset(idx_sac-1);
                SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.start_x;
                SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.start_y;
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
                SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
                SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
                SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
                SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_onset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_onset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
                SACS_ALL_TRIAL.diff_ang(idx_sac) = abs(acosd( ...
                    ( (SACS_ALL_TRIAL.eye_amp_x(idx_sac) .* SACS_ALL_TRIAL.visual_amp_x(idx_sac)) + ...
                    (SACS_ALL_TRIAL.eye_amp_y(idx_sac) .* SACS_ALL_TRIAL.visual_amp_y(idx_sac)) ) ...
                    ./ (SACS_ALL_TRIAL.eye_amp_m(idx_sac)) ./ (SACS_ALL_TRIAL.visual_amp_m(idx_sac)) ));
                SACS_ALL_TRIAL.tag(idx_sac) = 7; % 'back_center_prim' tag 7
            end
        end
    end
    
    % Tag target_irrelev. These are the remaining saccades that the saccade offset has landed on the target.
    diff_finish = sqrt( ...
        ((SACS_ALL_TRIAL.eye_px_offset - SACS_ALL_TRIAL.tgt_px_offset).^2) + ...
        ((SACS_ALL_TRIAL.eye_py_offset - SACS_ALL_TRIAL.tgt_py_offset).^2) );
    idx_target_irrelev = isnan(SACS_ALL_TRIAL.tag) & (diff_finish < threshold_pos);
    SACS_ALL_TRIAL.visual_px_onset( idx_target_irrelev) = SACS_ALL_TRIAL.eye_px_onset(idx_target_irrelev);
    SACS_ALL_TRIAL.visual_py_onset( idx_target_irrelev) = SACS_ALL_TRIAL.eye_py_onset(idx_target_irrelev); 
    SACS_ALL_TRIAL.visual_px_offset(idx_target_irrelev) = SACS_ALL_TRIAL.tgt_px_offset( idx_target_irrelev);
    SACS_ALL_TRIAL.visual_py_offset(idx_target_irrelev) = SACS_ALL_TRIAL.tgt_py_offset( idx_target_irrelev);
    SACS_ALL_TRIAL.diff_finish(     idx_target_irrelev) = diff_finish(idx_target_irrelev);
    SACS_ALL_TRIAL.tag(             idx_target_irrelev) = 9; % 'target_irrelev' tag 9
    
    % Tag the remaining saccades as other_irrelev
    idx_other_irrelev = isnan(SACS_ALL_TRIAL.tag);
    SACS_ALL_TRIAL.visual_px_onset( idx_other_irrelev) = SACS_ALL_TRIAL.eye_px_onset( idx_other_irrelev); % this entry is for the sake of filling up the nan values
    SACS_ALL_TRIAL.visual_py_onset( idx_other_irrelev) = SACS_ALL_TRIAL.eye_py_onset( idx_other_irrelev); % this entry is for the sake of filling up the nan values
    SACS_ALL_TRIAL.visual_px_offset(idx_other_irrelev) = SACS_ALL_TRIAL.eye_px_offset(idx_other_irrelev); % this entry is for the sake of filling up the nan values
    SACS_ALL_TRIAL.visual_py_offset(idx_other_irrelev) = SACS_ALL_TRIAL.eye_py_offset(idx_other_irrelev); % this entry is for the sake of filling up the nan values
    SACS_ALL_TRIAL.diff_finish(     idx_other_irrelev) = 0.0; % this entry is for the sake of filling up the nan values

    SACS_ALL_TRIAL.tag(             idx_other_irrelev) = 10; % 'other_irrelev' tag 10
    
    % Tag the prim. sac. that is not followed by corr. sac.
    % Sac. is prim_success, but visual info. is what occurs at its offset
    % Onset time is either the end of trial or onset of next sac.
    % Offset time is the same as that of prim_success: offset of sac.
    % Visual time the offset time
    num_saccades = length(SACS_ALL_TRIAL.validity); % number of detected saccades; this should not change
    idx_prim_success = find(SACS_ALL_TRIAL.tag == 1); % find prim_success
    if ~isempty(idx_prim_success)
        for counter_prim = 1 : length(idx_prim_success)
            idx_prim = idx_prim_success(counter_prim);
            next_idx_prim = idx_prim + 1; % index after the prim_success
            diff_finish = sqrt( ...
                ((SACS_ALL_TRIAL.tgt_px_offset( idx_prim) - SACS_ALL_TRIAL.eye_px_offset( idx_prim)).^2) + ...
                ((SACS_ALL_TRIAL.tgt_py_offset( idx_prim) - SACS_ALL_TRIAL.eye_py_offset( idx_prim)).^2) );
            % If (the prim. is the last saccade of trial || (the next saccade exists && isn't corr_success)) &&
            % @ prim. sac. offset, tgt. is at diff. pos. from where it lands, 
            % then it is counted as 'prim_no_corr'
            if (((next_idx_prim > num_saccades) || ...
                    ((next_idx_prim <= num_saccades) && ~(SACS_ALL_TRIAL.tag(next_idx_prim) == 4))) && ...
                    (diff_finish > threshold_pos))
                % Append the saccade data, mostly using the the prim. data
                SACS_ALL_TRIAL.validity(end+1) = SACS_ALL_TRIAL.validity(idx_prim);
                SACS_ALL_TRIAL.trial_num(end+1) = SACS_ALL_TRIAL.trial_num(idx_prim);
                SACS_ALL_TRIAL.tag(end+1) = 11;
                SACS_ALL_TRIAL.count(end+1) = SACS_ALL_TRIAL.count(idx_prim);
                SACS_ALL_TRIAL.flag_last_cue(end+1) = SACS_ALL_TRIAL.flag_last_cue(idx_prim);
                SACS_ALL_TRIAL.start_x(end+1) = SACS_ALL_TRIAL.start_x(idx_prim);
                SACS_ALL_TRIAL.start_y(end+1) = SACS_ALL_TRIAL.start_y(idx_prim);
                SACS_ALL_TRIAL.cue_x(end+1) = SACS_ALL_TRIAL.cue_x(idx_prim);
                SACS_ALL_TRIAL.cue_y(end+1) = SACS_ALL_TRIAL.cue_y(idx_prim);
                SACS_ALL_TRIAL.end_x(end+1) = SACS_ALL_TRIAL.end_x(idx_prim);
                SACS_ALL_TRIAL.end_y(end+1) = SACS_ALL_TRIAL.end_y(idx_prim);
                SACS_ALL_TRIAL.time_auditory(end+1) = SACS_ALL_TRIAL.time_auditory(idx_prim);
                SACS_ALL_TRIAL.time_vmax(end+1) = SACS_ALL_TRIAL.time_vmax(idx_prim);
                if (next_idx_prim > num_saccades)
                    SACS_ALL_TRIAL.time_onset(end+1) = TRIAL.time_end;
                else
                    SACS_ALL_TRIAL.time_onset(end+1) = SACS_ALL_TRIAL.time_onset(next_idx_prim);
                end
                SACS_ALL_TRIAL.time_visual(end+1) = SACS_ALL_TRIAL.time_offset(idx_prim);
                SACS_ALL_TRIAL.time_offset(end+1) = SACS_ALL_TRIAL.time_offset(idx_prim);
                SACS_ALL_TRIAL.visual_px_onset(end+1) = SACS_ALL_TRIAL.eye_px_offset(idx_prim); 
                SACS_ALL_TRIAL.visual_py_onset(end+1) = SACS_ALL_TRIAL.eye_py_offset(idx_prim);
                SACS_ALL_TRIAL.visual_px_offset(end+1) = SACS_ALL_TRIAL.tgt_px_offset(idx_prim);
                SACS_ALL_TRIAL.visual_py_offset(end+1) = SACS_ALL_TRIAL.tgt_py_offset(idx_prim);
                SACS_ALL_TRIAL.reaction(end+1) = SACS_ALL_TRIAL.reaction(idx_prim);
                SACS_ALL_TRIAL.visual_amp_x(end+1) = (SACS_ALL_TRIAL.visual_px_offset(next_idx_prim) - SACS_ALL_TRIAL.visual_px_onset(next_idx_prim));
                SACS_ALL_TRIAL.visual_amp_y(end+1) = (SACS_ALL_TRIAL.visual_py_offset(next_idx_prim) - SACS_ALL_TRIAL.visual_py_onset(next_idx_prim));
                SACS_ALL_TRIAL.visual_amp_m(end+1) = sqrt((SACS_ALL_TRIAL.visual_amp_x(next_idx_prim).^2) + (SACS_ALL_TRIAL.visual_amp_y(next_idx_prim).^2));
                SACS_ALL_TRIAL.visual_ang(end+1) = atan2d(SACS_ALL_TRIAL.visual_amp_y(next_idx_prim), SACS_ALL_TRIAL.visual_amp_x(next_idx_prim));
                SACS_ALL_TRIAL.diff_start(end+1) = SACS_ALL_TRIAL.diff_start(idx_prim);
                SACS_ALL_TRIAL.diff_finish(end+1) = SACS_ALL_TRIAL.diff_finish(idx_prim);
                SACS_ALL_TRIAL.diff_ang(end+1) = SACS_ALL_TRIAL.diff_ang(idx_prim);
                SACS_ALL_TRIAL.duration(end+1) = SACS_ALL_TRIAL.duration(idx_prim);
                SACS_ALL_TRIAL.eye_vm_max(end+1) = SACS_ALL_TRIAL.eye_vm_max(idx_prim);
                % Offset of prim. sac. as 'eye onset'; 'eye offset' the same 
                SACS_ALL_TRIAL.eye_px_onset(end+1) = SACS_ALL_TRIAL.eye_px_offset(idx_prim);
                SACS_ALL_TRIAL.eye_py_onset(end+1) = SACS_ALL_TRIAL.eye_py_offset(idx_prim);
                SACS_ALL_TRIAL.eye_px_offset(end+1) = SACS_ALL_TRIAL.eye_px_offset(idx_prim);
                SACS_ALL_TRIAL.eye_py_offset(end+1) = SACS_ALL_TRIAL.eye_py_offset(idx_prim);
                SACS_ALL_TRIAL.eye_amp_x(end+1) = SACS_ALL_TRIAL.eye_amp_x(idx_prim);
                SACS_ALL_TRIAL.eye_amp_y(end+1) = SACS_ALL_TRIAL.eye_amp_y(idx_prim);
                SACS_ALL_TRIAL.eye_amp_m(end+1) = SACS_ALL_TRIAL.eye_amp_m(idx_prim);
                SACS_ALL_TRIAL.eye_ang(end+1) = SACS_ALL_TRIAL.eye_ang(idx_prim);                
                SACS_ALL_TRIAL.tgt_px_onset(end+1) = SACS_ALL_TRIAL.tgt_px_onset(idx_prim);
                SACS_ALL_TRIAL.tgt_px_offset(end+1) = SACS_ALL_TRIAL.tgt_px_offset(idx_prim);
                SACS_ALL_TRIAL.tgt_py_onset(end+1) = SACS_ALL_TRIAL.tgt_py_onset(idx_prim);
                SACS_ALL_TRIAL.tgt_py_offset(end+1) = SACS_ALL_TRIAL.tgt_py_offset(idx_prim); 
                % REWARD PARAMS
                SACS_ALL_TRIAL.tgt_num(end+1) = SACS_ALL_TRIAL.tgt_num(idx_prim); 
                SACS_ALL_TRIAL.cue_x_high_rew(end+1) = SACS_ALL_TRIAL.cue_x_high_rew(idx_prim); 
                SACS_ALL_TRIAL.cue_y_high_rew(end+1) = SACS_ALL_TRIAL.cue_y_high_rew(idx_prim);
                SACS_ALL_TRIAL.cue_x_low_rew(end+1) = SACS_ALL_TRIAL.cue_x_low_rew(idx_prim);
                SACS_ALL_TRIAL.cue_y_low_rew(end+1) = SACS_ALL_TRIAL.cue_y_low_rew(idx_prim);
                SACS_ALL_TRIAL.high_rew_tgt_num(end+1) = SACS_ALL_TRIAL.high_rew_tgt_num(idx_prim);
                SACS_ALL_TRIAL.low_rew_tgt_num(end+1) = SACS_ALL_TRIAL.low_rew_tgt_num(idx_prim);
                SACS_ALL_TRIAL.task_cond(end+1) = SACS_ALL_TRIAL.task_cond(idx_prim);
                SACS_ALL_TRIAL.choice(end+1) = SACS_ALL_TRIAL.choice(idx_prim);
                SACS_ALL_TRIAL.rew_cond(end+1) = SACS_ALL_TRIAL.rew_cond(idx_prim);
                SACS_ALL_TRIAL.jump_cond(end+1) = SACS_ALL_TRIAL.jump_cond(idx_prim);
                SACS_ALL_TRIAL.tgt_cond(end+1) = SACS_ALL_TRIAL.tgt_cond(idx_prim);
                % Backstep exp. parameters
                if isfield(TRIAL, 'time_state_dwell_on')
                    SACS_ALL_TRIAL.time_state_dwell_on(end+1) = SACS_ALL_TRIAL.time_state_dwell_on(idx_prim);
                    SACS_ALL_TRIAL.time_dwell(end+1) = SACS_ALL_TRIAL.time_dwell(idx_prim);
                    SACS_ALL_TRIAL.time_state_dwell_off(end+1) = SACS_ALL_TRIAL.time_state_dwell_off(idx_prim);
                end
            end
        end
    end
      

    % Fill up the nan values
    idx_nan_visual_values = isnan(SACS_ALL_TRIAL.visual_amp_x);
    SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values)     = (SACS_ALL_TRIAL.visual_px_offset(idx_nan_visual_values) - SACS_ALL_TRIAL.visual_px_onset(idx_nan_visual_values));
    SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values)     = (SACS_ALL_TRIAL.visual_py_offset(idx_nan_visual_values) - SACS_ALL_TRIAL.visual_py_onset(idx_nan_visual_values));
    SACS_ALL_TRIAL.visual_amp_m(idx_nan_visual_values)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values).^2));
    SACS_ALL_TRIAL.visual_ang(idx_nan_visual_values)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values), SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values));
    SACS_ALL_TRIAL.diff_start(idx_nan_visual_values)  = sqrt( ...
        ((SACS_ALL_TRIAL.eye_px_onset( idx_nan_visual_values) - SACS_ALL_TRIAL.visual_px_onset( idx_nan_visual_values)).^2) + ...
        ((SACS_ALL_TRIAL.eye_py_onset( idx_nan_visual_values) - SACS_ALL_TRIAL.visual_py_onset( idx_nan_visual_values)).^2) );
    SACS_ALL_TRIAL.diff_ang(idx_nan_visual_values) = abs(acosd( ...
        ( (SACS_ALL_TRIAL.eye_amp_x(idx_nan_visual_values) .* SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values)) + ...
          (SACS_ALL_TRIAL.eye_amp_y(idx_nan_visual_values) .* SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values)) ) ...
        ./ (SACS_ALL_TRIAL.eye_amp_m(idx_nan_visual_values)) ./ (SACS_ALL_TRIAL.visual_amp_m(idx_nan_visual_values)) ));
    
    
    % Set the reaction time of target_irrelev and other_irrelev as the time
    % from the offset of the previous saccade.
    idx_irrelev = find((SACS_ALL_TRIAL.tag==9) | (SACS_ALL_TRIAL.tag==10));
     % 'target_irrelev' tag 9 % 'other_irrelev' tag 10
    if ~isempty(idx_irrelev)
        for counter_target_irrelev = 1 : length(idx_irrelev)
            idx_sac = idx_irrelev(counter_target_irrelev);
            if idx_sac == 1
                SACS_ALL_TRIAL.time_visual(idx_sac)      = TRIAL.time_state_str_pursuit(idx_sac);
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
            else
                SACS_ALL_TRIAL.time_visual(idx_sac)      = SACS_ALL_TRIAL.time_offset(idx_sac-1);
                SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
            end
        end
    end
    
    % Re-tag the irrelevant sac. possibly as db_corr_success; most are target_irrelevant sac. but some were found to be classified as 
    % irrelevant sac. bc. it comes late in the trial, so the tgt. has jumped back to start by the time db_corr sac. lands at cue tgt. 
    % @ the previous corr_success offset, tgt. should have jumped to a different location from where corr_success was aiming;
    % dwell time could be long enough that tgt. doesn't jump back to cue bf. the first corr_success ends, in which case,
    % there is no "corrective" saccade to be made
    idx_irrelev = find((SACS_ALL_TRIAL.tag==9) | (SACS_ALL_TRIAL.tag==10));
    if ~isempty(idx_irrelev)
        for counter_irrelev = 1 : length(idx_irrelev)
            idx_sac = idx_irrelev(counter_irrelev);
            prev_idx_sac = idx_sac - 1; % index before the target_irrelevant sac.
            % Make sure there is a previous saccade
            if (prev_idx_sac >= 1) 
                start_pos_drift = sqrt((SACS_ALL_TRIAL.eye_px_onset(idx_sac)-SACS_ALL_TRIAL.eye_px_offset(prev_idx_sac)).^2 + ...
                    (SACS_ALL_TRIAL.eye_py_onset(idx_sac)-SACS_ALL_TRIAL.eye_py_offset(prev_idx_sac)).^2); 
                sac_offset_dist_to_prev_sac_tgt = sqrt((SACS_ALL_TRIAL.eye_px_offset(idx_sac)-SACS_ALL_TRIAL.tgt_px_offset(prev_idx_sac)).^2 + ...
                    (SACS_ALL_TRIAL.eye_py_offset(idx_sac)-SACS_ALL_TRIAL.tgt_py_offset(prev_idx_sac)).^2);
                diff_finish = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_offset(prev_idx_sac) - SACS_ALL_TRIAL.tgt_px_offset(prev_idx_sac)).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_offset(prev_idx_sac) - SACS_ALL_TRIAL.tgt_py_offset(prev_idx_sac)).^2) );
                % If the previous sac. is corr_success && current sac. lands where tgt. was @ prev. sac. offset &&
                % eye didn't drift away from where it landed @ prev. sac. offset when it starts &&
                % tgt. @ prev. sac. offset is away from where eye lands, then tag 12
                if  ((SACS_ALL_TRIAL.tag(prev_idx_sac) == 4) && (sac_offset_dist_to_prev_sac_tgt < threshold_pos) && ...
                        (start_pos_drift < 0.5) && (diff_finish > 1.0))
                    SACS_ALL_TRIAL.tag(idx_sac) = 12;
                end
            end     
        end
    end
   
    % Tag the corr. sac. that is not followed by another corr. sac.
    % Sac. is corr_success, but visual info. is what occurs at its offset
    % Onset time is either the end of the trial or onset of next saccade
    % Offset time is the same as that of corr_success: offset of sac.
    % Visual time is the offset time
    % @ its offset, tgt. should have jumped to a different location from where corr_success was aiming;
    % dwell time could be long enough that tgt. doesn't jump back to cue bf. the first corr_success ends, in which case,
    % there is no "corrective" saccade to be made
    idx_corr_success = find(SACS_ALL_TRIAL.tag == 4); % find corr_success
    if ~isempty(idx_corr_success)
        for counter_corr = 1 : length(idx_corr_success)
            idx_corr = idx_corr_success(counter_corr);
            next_idx_corr = idx_corr + 1; % index after the corr_success
            diff_finish = sqrt( ...
                ((SACS_ALL_TRIAL.tgt_px_offset( idx_corr) - SACS_ALL_TRIAL.eye_px_offset( idx_corr)).^2) + ...
                ((SACS_ALL_TRIAL.tgt_py_offset( idx_corr) - SACS_ALL_TRIAL.eye_py_offset( idx_corr)).^2) );
            % If (the corr. is the last saccade of trial || (the next saccade exists && isn't db_corr_success)) &&
            % && tgt. @ prev. sac. offset is away from where eye lands, then count as tag 13
            if (((next_idx_corr > num_saccades) || ((next_idx_corr <= num_saccades) && (SACS_ALL_TRIAL.tag(next_idx_corr) ~= 12))) && ...
                    (diff_finish > 1.0))
                SACS_ALL_TRIAL.validity(end+1) = SACS_ALL_TRIAL.validity(idx_corr);
                SACS_ALL_TRIAL.trial_num(end+1) = SACS_ALL_TRIAL.trial_num(idx_corr);
                SACS_ALL_TRIAL.tag(end+1) = 13;
                SACS_ALL_TRIAL.count(end+1) = SACS_ALL_TRIAL.count(idx_corr);
                SACS_ALL_TRIAL.flag_last_cue(end+1) = SACS_ALL_TRIAL.flag_last_cue(idx_corr);
                SACS_ALL_TRIAL.start_x(end+1) = SACS_ALL_TRIAL.start_x(idx_corr);
                SACS_ALL_TRIAL.start_y(end+1) = SACS_ALL_TRIAL.start_y(idx_corr);
                SACS_ALL_TRIAL.cue_x(end+1) = SACS_ALL_TRIAL.cue_x(idx_corr);
                SACS_ALL_TRIAL.cue_y(end+1) = SACS_ALL_TRIAL.cue_y(idx_corr);
                SACS_ALL_TRIAL.end_x(end+1) = SACS_ALL_TRIAL.end_x(idx_corr);
                SACS_ALL_TRIAL.end_y(end+1) = SACS_ALL_TRIAL.end_y(idx_corr);
                SACS_ALL_TRIAL.time_auditory(end+1) = SACS_ALL_TRIAL.time_auditory(idx_corr); 
                SACS_ALL_TRIAL.time_vmax(end+1) = SACS_ALL_TRIAL.time_vmax(idx_corr);
                if (next_idx_corr > num_saccades)
                    SACS_ALL_TRIAL.time_onset(end+1) = TRIAL.time_end;
                else
                    SACS_ALL_TRIAL.time_onset(end+1) = SACS_ALL_TRIAL.time_onset(next_idx_corr);
                end
                SACS_ALL_TRIAL.time_visual(end+1) = SACS_ALL_TRIAL.time_offset(idx_corr);
                SACS_ALL_TRIAL.time_offset(end+1) = SACS_ALL_TRIAL.time_offset(idx_corr);
                
                SACS_ALL_TRIAL.visual_px_onset(end+1) = SACS_ALL_TRIAL.eye_px_offset(idx_corr); 
                SACS_ALL_TRIAL.visual_py_onset(end+1) = SACS_ALL_TRIAL.eye_py_offset(idx_corr);
                SACS_ALL_TRIAL.visual_px_offset(end+1) = SACS_ALL_TRIAL.tgt_px_offset(idx_corr);
                SACS_ALL_TRIAL.visual_py_offset(end+1) = SACS_ALL_TRIAL.tgt_py_offset(idx_corr);
                SACS_ALL_TRIAL.reaction(end+1) = SACS_ALL_TRIAL.reaction(idx_corr);
                SACS_ALL_TRIAL.visual_amp_x(end+1) = (SACS_ALL_TRIAL.visual_px_offset(next_idx_corr) - SACS_ALL_TRIAL.visual_px_onset(next_idx_corr));
                SACS_ALL_TRIAL.visual_amp_y(end+1) = (SACS_ALL_TRIAL.visual_py_offset(next_idx_corr) - SACS_ALL_TRIAL.visual_py_onset(next_idx_corr));
                SACS_ALL_TRIAL.visual_amp_m(end+1) = sqrt((SACS_ALL_TRIAL.visual_amp_x(next_idx_corr).^2) + (SACS_ALL_TRIAL.visual_amp_y(next_idx_corr).^2));
                SACS_ALL_TRIAL.visual_ang(end+1) = atan2d(SACS_ALL_TRIAL.visual_amp_y(next_idx_corr), SACS_ALL_TRIAL.visual_amp_x(next_idx_corr));
                SACS_ALL_TRIAL.diff_start(end+1) = SACS_ALL_TRIAL.diff_start(idx_corr);
                SACS_ALL_TRIAL.diff_finish(end+1) = SACS_ALL_TRIAL.diff_finish(idx_corr);
                SACS_ALL_TRIAL.diff_ang(end+1) = SACS_ALL_TRIAL.diff_ang(idx_corr);
                SACS_ALL_TRIAL.duration(end+1) = SACS_ALL_TRIAL.duration(idx_corr);
                SACS_ALL_TRIAL.eye_vm_max(end+1) = SACS_ALL_TRIAL.eye_vm_max(idx_corr);
                % Offset of corr. sac. as 'eye onset'; 'eye offset' the same 
                SACS_ALL_TRIAL.eye_px_onset(end+1) = SACS_ALL_TRIAL.eye_px_offset(idx_corr);
                SACS_ALL_TRIAL.eye_py_onset(end+1) = SACS_ALL_TRIAL.eye_py_offset(idx_corr);
                SACS_ALL_TRIAL.eye_px_offset(end+1) = SACS_ALL_TRIAL.eye_px_offset(idx_corr);
                SACS_ALL_TRIAL.eye_py_offset(end+1) = SACS_ALL_TRIAL.eye_py_offset(idx_corr);
                SACS_ALL_TRIAL.eye_amp_x(end+1) = SACS_ALL_TRIAL.eye_amp_x(idx_corr);
                SACS_ALL_TRIAL.eye_amp_y(end+1) = SACS_ALL_TRIAL.eye_amp_y(idx_corr);
                SACS_ALL_TRIAL.eye_amp_m(end+1) = SACS_ALL_TRIAL.eye_amp_m(idx_corr);
                SACS_ALL_TRIAL.eye_ang(end+1) = SACS_ALL_TRIAL.eye_ang(idx_corr);                
                SACS_ALL_TRIAL.tgt_px_onset(end+1) = SACS_ALL_TRIAL.tgt_px_onset(idx_corr);
                SACS_ALL_TRIAL.tgt_px_offset(end+1) = SACS_ALL_TRIAL.tgt_px_offset(idx_corr);
                SACS_ALL_TRIAL.tgt_py_onset(end+1) = SACS_ALL_TRIAL.tgt_py_onset(idx_corr);
                SACS_ALL_TRIAL.tgt_py_offset(end+1) = SACS_ALL_TRIAL.tgt_py_offset(idx_corr);
                % reward params
                SACS_ALL_TRIAL.tgt_num(end+1) = SACS_ALL_TRIAL.tgt_num(idx_corr); 
                SACS_ALL_TRIAL.cue_x_high_rew(end+1) = SACS_ALL_TRIAL.cue_x_high_rew(idx_corr); 
                SACS_ALL_TRIAL.cue_y_high_rew(end+1) = SACS_ALL_TRIAL.cue_y_high_rew(idx_corr);
                SACS_ALL_TRIAL.cue_x_low_rew(end+1) = SACS_ALL_TRIAL.cue_x_low_rew(idx_corr);
                SACS_ALL_TRIAL.cue_y_low_rew(end+1) = SACS_ALL_TRIAL.cue_y_low_rew(idx_corr);
                SACS_ALL_TRIAL.high_rew_tgt_num(end+1) = SACS_ALL_TRIAL.high_rew_tgt_num(idx_corr);
                SACS_ALL_TRIAL.low_rew_tgt_num(end+1) = SACS_ALL_TRIAL.low_rew_tgt_num(idx_corr);
                SACS_ALL_TRIAL.task_cond(end+1) = SACS_ALL_TRIAL.task_cond(idx_corr);
                SACS_ALL_TRIAL.choice(end+1) = SACS_ALL_TRIAL.choice(idx_corr);
                SACS_ALL_TRIAL.rew_cond(end+1) = SACS_ALL_TRIAL.rew_cond(idx_corr);
                SACS_ALL_TRIAL.jump_cond(end+1) = SACS_ALL_TRIAL.jump_cond(idx_corr);
                SACS_ALL_TRIAL.tgt_cond(end+1) = SACS_ALL_TRIAL.tgt_cond(idx_corr);
                % Backstep exp. parameters
                if isfield(TRIAL, 'time_state_dwell_on')
                    SACS_ALL_TRIAL.time_state_dwell_on(end+1) = SACS_ALL_TRIAL.time_state_dwell_on(idx_corr);
                    SACS_ALL_TRIAL.time_dwell(end+1) = SACS_ALL_TRIAL.time_dwell(idx_corr);
                    SACS_ALL_TRIAL.time_state_dwell_off(end+1) = SACS_ALL_TRIAL.time_state_dwell_off(idx_corr);
                end
            end
        end
    end
      
    % Re-tag the 1st saccade as back_center_success if it landed at start
    % Handles case when there is 'back_center_success' but eye makes irrelevant saccades before the saccade twrd. target that initiates cue present.
    idx_sac = 1;
    diff_finish  = sqrt( ...
        ((SACS_ALL_TRIAL.eye_px_offset( idx_sac) - TRIAL.start_x).^2) + ...
        ((SACS_ALL_TRIAL.eye_py_offset( idx_sac) - TRIAL.start_y).^2) );
    if(diff_finish < threshold_pos)
        
        SACS_ALL_TRIAL.time_visual(idx_sac)      = TRIAL.time_state_str_pursuit(idx_sac);
        SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = SACS_ALL_TRIAL.eye_px_onset(idx_sac);
        SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = SACS_ALL_TRIAL.eye_py_onset(idx_sac);
        SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.start_x;
        SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.start_y;
        SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
        SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = (SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac));
        SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = (SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac));
        SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = sqrt((SACS_ALL_TRIAL.visual_amp_x(idx_sac).^2) + (SACS_ALL_TRIAL.visual_amp_y(idx_sac).^2));
        SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));
        SACS_ALL_TRIAL.diff_start(idx_sac)  = sqrt( ...
            ((SACS_ALL_TRIAL.eye_px_onset( idx_sac) - SACS_ALL_TRIAL.visual_px_onset( idx_sac)).^2) + ...
            ((SACS_ALL_TRIAL.eye_py_onset( idx_sac) - SACS_ALL_TRIAL.visual_py_onset( idx_sac)).^2) );
        SACS_ALL_TRIAL.diff_finish(idx_sac) = sqrt( ...
            ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac)).^2) + ...
            ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac)).^2) );
        SACS_ALL_TRIAL.diff_ang(idx_sac) = abs(acosd( ...
            ( (SACS_ALL_TRIAL.eye_amp_x(idx_sac) .* SACS_ALL_TRIAL.visual_amp_x(idx_sac)) + ...
              (SACS_ALL_TRIAL.eye_amp_y(idx_sac) .* SACS_ALL_TRIAL.visual_amp_y(idx_sac)) ) ...
            ./ (SACS_ALL_TRIAL.eye_amp_m(idx_sac)) ./ (SACS_ALL_TRIAL.visual_amp_m(idx_sac)) ));
        SACS_ALL_TRIAL.tag(idx_sac) = 6; % 'back_center_success' tag 6
    end
    
    
    
    % Fill up the count values
    prim_tags = (SACS_ALL_TRIAL.tag==1) | (SACS_ALL_TRIAL.tag==2); % 'prim_success' tag 1 % 'prim_attempt' tag 2
    if sum(prim_tags) > 0
        SACS_ALL_TRIAL.count(prim_tags) = 1 : sum(prim_tags);
    end
    prim_no_corr_tags = SACS_ALL_TRIAL.tag == 11; % 'prim_no_corr' tag 11 
    if sum(prim_tags) > 0
        SACS_ALL_TRIAL.count(prim_no_corr_tags) = 1 : sum(prim_no_corr_tags);
    end
    corr_tags = (SACS_ALL_TRIAL.tag==4) | (SACS_ALL_TRIAL.tag==5); % 'corr_success' tag 4 % 'cord_fail' tag 5
    if sum(corr_tags) > 0
        SACS_ALL_TRIAL.count(corr_tags) = 1 : sum(corr_tags);
    end
    db_corr_tags = SACS_ALL_TRIAL.tag == 12; % 'db_corr_sac' tag 12
    if sum(db_corr_tags) > 0
        SACS_ALL_TRIAL.count(db_corr_tags) = 1 : sum(db_corr_tags);
    end 
    back_center_tags = (SACS_ALL_TRIAL.tag==6) | (SACS_ALL_TRIAL.tag==7); % 'back_center_success' tag 6 % 'back_center_prim' tag 7
    if sum(back_center_tags) > 0
        SACS_ALL_TRIAL.count(back_center_tags) = 1 : sum(back_center_tags);
    end
    irrelev_tags = (SACS_ALL_TRIAL.tag==3) | (SACS_ALL_TRIAL.tag==8) | (SACS_ALL_TRIAL.tag==9) | (SACS_ALL_TRIAL.tag==10);
    % 'prim_fail' tag 3 % 'back_center_irrelev' tag 8 % 'target_irrelev' tag 9 % 'other_irrelev' tag 10
    if sum(irrelev_tags) > 0
        SACS_ALL_TRIAL.count(irrelev_tags) = 1 : sum(irrelev_tags);
    end
    
    % Make another tag for 'other_irrelev' where the visual ang. is 
    % based on direction of center target from the offset pos. of previous saccade
    SACS_fields = fieldnames(SACS_ALL_TRIAL);
    idx_other_irrelev = find(SACS_ALL_TRIAL.tag == 10);
    if ~isempty(idx_other_irrelev)
        for counter_other_irrelev = 1 : length(idx_other_irrelev)
            sac_idx = idx_other_irrelev(counter_other_irrelev);
            next_sac_idx = sac_idx + 1;
            % Copy the data first
            for counter_field = 1 : length(SACS_fields)
                field_name = SACS_fields{counter_field};
                SACS_ALL_TRIAL.(field_name)(end+1) = SACS_ALL_TRIAL.(field_name)(sac_idx);
            end
            % Change the data specific to this tag
            new_tag_idx = length(SACS_ALL_TRIAL.(field_name));
            SACS_ALL_TRIAL.tag(new_tag_idx) = 14;
            SACS_ALL_TRIAL.time_visual(new_tag_idx) = SACS_ALL_TRIAL.time_offset(sac_idx);
            SACS_ALL_TRIAL.time_onset(new_tag_idx) = SACS_ALL_TRIAL.time_onset(next_sac_idx);
            SACS_ALL_TRIAL.visual_px_offset(new_tag_idx) = SACS_ALL_TRIAL.tgt_px_offset(sac_idx); % for CS on
            SACS_ALL_TRIAL.visual_py_offset(new_tag_idx) = SACS_ALL_TRIAL.tgt_py_offset(sac_idx); % for CS on
            SACS_ALL_TRIAL.visual_px_onset(new_tag_idx) = SACS_ALL_TRIAL.eye_px_offset(sac_idx);
            SACS_ALL_TRIAL.visual_py_onset(new_tag_idx) = SACS_ALL_TRIAL.eye_py_offset(sac_idx);
            SACS_ALL_TRIAL.visual_amp_x(new_tag_idx) = SACS_ALL_TRIAL.visual_px_offset(new_tag_idx) - SACS_ALL_TRIAL.visual_px_onset(new_tag_idx);
            SACS_ALL_TRIAL.visual_amp_y(new_tag_idx) = SACS_ALL_TRIAL.visual_py_offset(new_tag_idx) - SACS_ALL_TRIAL.visual_py_onset(new_tag_idx);
            SACS_ALL_TRIAL.visual_amp_m(new_tag_idx) = sqrt(SACS_ALL_TRIAL.visual_amp_x(new_tag_idx).^2 + SACS_ALL_TRIAL.visual_amp_y(new_tag_idx).^2);
            SACS_ALL_TRIAL.visual_ang(new_tag_idx) = atan2d(SACS_ALL_TRIAL.visual_amp_y(new_tag_idx), SACS_ALL_TRIAL.visual_amp_x(new_tag_idx));
            SACS_ALL_TRIAL.eye_px_onset(new_tag_idx) = SACS_ALL_TRIAL.eye_px_offset(sac_idx); % for CS on
            SACS_ALL_TRIAL.eye_py_onset(new_tag_idx) = SACS_ALL_TRIAL.eye_py_offset(sac_idx); % for CS on
            SACS_ALL_TRIAL.reaction(new_tag_idx) = (SACS_ALL_TRIAL.time_onset(sac_idx+1) - SACS_ALL_TRIAL.time_offset(sac_idx))*1000;
            SACS_ALL_TRIAL.diff_start(new_tag_idx)  = sqrt( ...
                ((SACS_ALL_TRIAL.tgt_px_onset( sac_idx) - SACS_ALL_TRIAL.eye_px_onset( sac_idx)).^2) + ...
                ((SACS_ALL_TRIAL.tgt_py_onset( sac_idx) - SACS_ALL_TRIAL.eye_py_onset( sac_idx)).^2) );
            SACS_ALL_TRIAL.diff_finish(new_tag_idx) = sqrt( ...
                ((SACS_ALL_TRIAL.tgt_px_offset(sac_idx) - SACS_ALL_TRIAL.eye_px_offset(sac_idx)).^2) + ...
                ((SACS_ALL_TRIAL.tgt_py_offset(sac_idx) - SACS_ALL_TRIAL.eye_py_offset(sac_idx)).^2) );
            SACS_ALL_TRIAL.diff_ang(new_tag_idx) = 0;
            SACS_ALL_TRIAL.eye_amp_x(new_tag_idx) = 0;
            SACS_ALL_TRIAL.eye_amp_y(new_tag_idx) = 0;
            SACS_ALL_TRIAL.eye_amp_m(new_tag_idx) = 0;
            SACS_ALL_TRIAL.eye_ang(new_tag_idx) = 0;
        end
    end
    
    % Make another tag for 'back_center_irrelev' where the visual ang. is 
    % based on direction of center target from the offset pos. of previous saccade
    % Visual time is the same as that of tag 8
    SACS_fields = fieldnames(SACS_ALL_TRIAL);
    idx_back_center_irrelev = find(SACS_ALL_TRIAL.tag == 8);
    if ~isempty(idx_back_center_irrelev)
        for counter_back_center_irrelev = 1 : length(idx_back_center_irrelev)
            sac_idx = idx_back_center_irrelev(counter_back_center_irrelev);
            prev_sac_idx = sac_idx - 1;
            % Copy the data first
            for counter_field = 1 : length(SACS_fields)
                field_name = SACS_fields{counter_field};
                SACS_ALL_TRIAL.(field_name)(end+1) = SACS_ALL_TRIAL.(field_name)(sac_idx);
            end
            % Change the data specific to this tag
            new_tag_idx = length(SACS_ALL_TRIAL.(field_name));
            SACS_ALL_TRIAL.tag(new_tag_idx) = 15;
            SACS_ALL_TRIAL.visual_px_onset(new_tag_idx)  = SACS_ALL_TRIAL.eye_px_offset(prev_sac_idx);
            SACS_ALL_TRIAL.visual_py_onset(new_tag_idx)  = SACS_ALL_TRIAL.eye_py_offset(prev_sac_idx);
            SACS_ALL_TRIAL.visual_px_offset(new_tag_idx) = TRIAL.start_x; % for CS on
            SACS_ALL_TRIAL.visual_py_offset(new_tag_idx) = TRIAL.start_y; % for CS on
            SACS_ALL_TRIAL.visual_amp_x(new_tag_idx) = SACS_ALL_TRIAL.visual_px_offset(new_tag_idx) - SACS_ALL_TRIAL.visual_px_onset(new_tag_idx);
            SACS_ALL_TRIAL.visual_amp_y(new_tag_idx) = SACS_ALL_TRIAL.visual_py_offset(new_tag_idx) - SACS_ALL_TRIAL.visual_py_onset(new_tag_idx);
            SACS_ALL_TRIAL.visual_amp_m(new_tag_idx) = sqrt(SACS_ALL_TRIAL.visual_amp_x(new_tag_idx).^2 + SACS_ALL_TRIAL.visual_amp_y(new_tag_idx).^2);
            SACS_ALL_TRIAL.visual_ang(new_tag_idx) = atan2d(SACS_ALL_TRIAL.visual_amp_y(new_tag_idx), SACS_ALL_TRIAL.visual_amp_x(new_tag_idx));
            SACS_ALL_TRIAL.eye_px_onset(new_tag_idx) = SACS_ALL_TRIAL.eye_px_offset(prev_sac_idx); % for CS on
            SACS_ALL_TRIAL.eye_py_onset(new_tag_idx) = SACS_ALL_TRIAL.eye_py_offset(prev_sac_idx); % for CS on
            SACS_ALL_TRIAL.diff_start(new_tag_idx)  = sqrt( ...
                ((SACS_ALL_TRIAL.tgt_px_onset( sac_idx) - SACS_ALL_TRIAL.eye_px_onset( sac_idx)).^2) + ...
                ((SACS_ALL_TRIAL.tgt_py_onset( sac_idx) - SACS_ALL_TRIAL.eye_py_onset( sac_idx)).^2) );
            SACS_ALL_TRIAL.eye_amp_x(new_tag_idx) = 0; % fill-in value
            SACS_ALL_TRIAL.eye_amp_y(new_tag_idx) = 0; % fill-in value
            SACS_ALL_TRIAL.eye_amp_m(new_tag_idx) = 0; % fill-in value
            SACS_ALL_TRIAL.eye_ang(new_tag_idx) = 0;   % fill-in value
            SACS_ALL_TRIAL.diff_ang(new_tag_idx) = 0;  % fill-in value
        end
    end
    
    
    %% Debugging
    flag_plot_trial = 0;
    %     flag_plot_trial = 1; % Jinyu Temp
    if flag_plot_trial %&& (counter_trial == trial_num)
        %trial_num = counter_trial;

%         eye_px_filt = TRIAL.eye_px_filt;
%         eye_py_filt = TRIAL.eye_py_filt;
%         eye_vm_filt = TRIAL.eye_vm_filt;
% 
%         hAx_(1) = subplot(3,1,1);
%         hold on;
%         plot(time_device, eye_px_filt)
%         plot(time_device, TRIAL.tgt_px)
%         %     plot(time, X, 'r') % Jinyu Temp
%         %     plot(t_raw, X_raw, 'm') % Jinyu Temp
%         plot(time_device(all_sac_ind_onset(all_sac_validity)), eye_px_filt(all_sac_ind_onset(all_sac_validity)), 'ok')
%         plot(time_device(all_sac_ind_offset(all_sac_validity)), eye_px_filt(all_sac_ind_offset(all_sac_validity)), 'ob')
%         plot(TRIAL.time_state_str_fixation'-TRIAL.time_pursuit, zeros(size(TRIAL.time_state_str_fixation')), '*k')
%         plot(TRIAL.time_state_sac_detect_on', zeros(size(TRIAL.time_state_sac_detect_on')), '*r')
%         plot(TRIAL.time_state_sac_detect_off', zeros(size(TRIAL.time_state_sac_detect_off')), '*m')
%         plot(TRIAL.time_state_reward', zeros(size(TRIAL.time_state_reward')), '*b')
%         %     plot(TRIAL.time_state_iti', zeros(size(TRIAL.time_state_iti')), '*g')
%         plot(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)', zeros(size(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)')), 'sk')
%         ylim([-12 12])
%         %     xlabel('Time (s)')
%         %     set(gca,'xtick',[])
%         ylabel('Horz. Eye (deg)')
%         %     title(['Trial no.: ' num2str(trial_num) ', (191111-142741)'])
% 
%         hAx_(2) = subplot(3,1,2);
%         hold on;
%         plot(time_device, eye_py_filt)
%         plot(time_device, TRIAL.tgt_py)
%         %     plot(time, Y, 'r') % Jinyu Temp
%         %     plot(t_raw, Y_raw, 'm') % Jinyu Temp
%         plot(time_device(all_sac_ind_onset(all_sac_validity)), eye_py_filt(all_sac_ind_onset(all_sac_validity)), 'ok')
%         plot(time_device(all_sac_ind_offset(all_sac_validity)), eye_py_filt(all_sac_ind_offset(all_sac_validity)), 'ob')
%         plot(TRIAL.time_state_str_fixation'-TRIAL.time_pursuit, zeros(size(TRIAL.time_state_str_fixation')), '*k')
%         plot(TRIAL.time_state_sac_detect_on', zeros(size(TRIAL.time_state_sac_detect_on')), '*r')
%         plot(TRIAL.time_state_sac_detect_off', zeros(size(TRIAL.time_state_sac_detect_off')), '*m')
%         plot(TRIAL.time_state_reward', zeros(size(TRIAL.time_state_reward')), '*b')
%         %     plot(TRIAL.time_state_iti', zeros(size(TRIAL.time_state_iti')), '*g')
%         plot(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)', zeros(size(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)')), 'sk')
%         ylim([-12 12])
%         %     xlabel('Time (s)')
%         %     set(gca,'xtick',[])
%         ylabel('Vert. Eye (deg)')
% 
%         hAx_(3) = subplot(3,1,3);
%         hold on;
%         plot(time_device, eye_vm_filt)
%         plot(time_device(all_sac_ind_onset(all_sac_validity)), eye_vm_filt(all_sac_ind_onset(all_sac_validity)), 'ok')
%         plot(time_device(all_sac_ind_vmax(all_sac_validity)), eye_vm_filt(all_sac_ind_vmax(all_sac_validity)), 'or')
%         plot(time_device(all_sac_ind_offset(all_sac_validity)), eye_vm_filt(all_sac_ind_offset(all_sac_validity)), 'ob')
%         plot(TRIAL.time_state_str_fixation'-TRIAL.time_pursuit, zeros(size(TRIAL.time_state_str_fixation')), '*k')
%         plot(TRIAL.time_state_sac_detect_on', zeros(size(TRIAL.time_state_sac_detect_on')), '*r')
%         plot(TRIAL.time_state_sac_detect_off', zeros(size(TRIAL.time_state_sac_detect_off')), '*m')
%         plot(TRIAL.time_state_reward', zeros(size(TRIAL.time_state_reward')), '*b')
%         %     plot(TRIAL.time_state_iti', zeros(size(TRIAL.time_state_iti')), '*g')
%         plot(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)', zeros(size(SACS_ALL_TRIAL.time_visual(SACS_ALL_TRIAL.tag==7)')), 'sk')
%         ylim([0 1000])
%         xlabel('Time (s)')
%         ylabel('Eye Vel. (deg/s)')
% 
%         for counter_sac = 1 : length(SACS_ALL_TRIAL.validity)
%             tag_ = SACS_ALL_TRIAL.tag(counter_sac);
%             if     tag_==1
%                 color_ = 'r';
%             elseif tag_==2
%                 color_ = [0.6350    0.0780    0.1840]; % 'r';
%             elseif tag_==3
%                 color_ = [0.4 0.4 0.4];
%             elseif tag_==4
%                 color_ = 'b';
%             elseif tag_==5
%                 color_ = [0    0.4470    0.7410]; % 'b';
%             elseif tag_==6
%                 color_ = [0.4940    0.1840    0.5560]; % 'm';
%             elseif tag_==7
%                 color_ = [0.4940    0.1840    0.5560]; % 'm';
%             elseif tag_==8
%                 color_ = [0.6 0.6 0.6];
%             elseif tag_==9
%                 color_ = [0.6 0.6 0.6];
%             elseif tag_==10
%                 color_ = [0.6 0.6 0.6];
%             elseif tag_ == 11
%                 color_ = [0.4660 0.6740 0.1880]; % greenish
%             elseif tag_ == 13
%                 color_ = [0.3010 0.7450 0.9330];
%             elseif ((tag_ == 14)) % for tag 14, ignore
%                 continue;
%             else
%                 color_ = [0.6 0.6 0.6];
%             end
%             label_ = [meta_data.sac_tag_list{tag_} '_' num2str(SACS_ALL_TRIAL.count(counter_sac))];
%             if ((tag_ == 11) || (tag_== 13))
%                 if tag_ == 11
%                     legend_text_ = 'prim no corr';
%                 elseif tag_ == 13
%                     legend_text_ = 'corr no db corr';
%                 end
%                 x_1 = xline(hAx_(1), SACS_ALL_TRIAL.time_offset(counter_sac), '-','LineWidth',2 ,'color', color_);
%                 legend(x_1,legend_text_);
%                 legend(hAx_(1),'boxoff');
%                 x_2 = xline(hAx_(2), SACS_ALL_TRIAL.time_offset(counter_sac), '-','LineWidth',2 ,'color', color_);
%                 legend(x_2,legend_text_);
%                 legend(hAx_(2),'boxoff');
%                 x_3 = xline(hAx_(3), SACS_ALL_TRIAL.time_offset(counter_sac), '-','LineWidth',2 ,'color', color_);
%                 legend(x_3,legend_text_);
%                 legend(hAx_(3),'boxoff');
% 
%             else
%                 xline(hAx_(1), SACS_ALL_TRIAL.time_onset(counter_sac), '-', label_, 'color', color_, 'interpret', 'none', 'FontSize', 12)
%                 xline(hAx_(2), SACS_ALL_TRIAL.time_onset(counter_sac), '-', label_, 'color', color_, 'interpret', 'none', 'FontSize', 12)
%                 xline(hAx_(3), SACS_ALL_TRIAL.time_onset(counter_sac), '-', label_, 'color', color_, 'interpret', 'none', 'FontSize', 12)
%             end
% 
% 
%             %         linkaxes(hAx_,'x')
%         end
%         %     w = waitforbuttonpress;
%         linkaxes(hAx_,'x');
%         ESN_Beautify_Plot(hFig_, [8,4], 12)
%         pause;
%         
        hFig_ = figure(1);
        clf(hFig_);
        disp(counter_trial)

        eye_px_filt = TRIAL.eye_px_filt;
        eye_py_filt = TRIAL.eye_py_filt;
        eye_vx_filt = TRIAL.eye_vx_filt;
        eye_vy_filt = TRIAL.eye_vy_filt;
        eye_vm_filt = TRIAL.eye_vm_filt;
        tgt_px      = TRIAL.tgt_px;
        tgt_py      = TRIAL.tgt_py;
        dist_to_cue = sqrt((TRIAL.cue_x - eye_px_filt).^2 + (TRIAL.cue_y - eye_py_filt).^2);
        dist_to_end = sqrt((TRIAL.end_x - eye_px_filt).^2 + (TRIAL.end_y - eye_py_filt).^2);
        dist_to_cue_low_rew  = sqrt((TRIAL.cue_x_low_rew - eye_px_filt).^2 + (TRIAL.cue_y_low_rew - eye_py_filt).^2);
        dist_to_cue_high_rew = sqrt((TRIAL.cue_x_high_rew - eye_px_filt).^2 + (TRIAL.cue_y_high_rew - eye_py_filt).^2);
        dist_from_start_tgt = sqrt((TRIAL.start_x - eye_px_filt).^2 + (TRIAL.start_y - eye_py_filt).^2);
        hAx_(1) = subplot(2,1,1);
        hold on;
        p_eye_px = plot(time_device, eye_px_filt, '-','Color','r','LineWidth',1);
        p_eye_py = plot(time_device, eye_py_filt, '-','Color','b','LineWidth',1);
        p_tgt_px = plot(time_device, tgt_px, '--', 'Color','r');
        p_tgt_py = plot(time_device, tgt_py, '--', 'Color','b');
        
        if TRIAL.task_cond == 0
            if TRIAL.choice == 1
                p_dist_to_cue = plot(time_device, dist_to_cue_high_rew,'k', 'LineWidth',2);
            else
                p_dist_to_cue = plot(time_device, dist_to_cue_low_rew,'k', 'LineWidth',2);
            end
        else
            p_dist_to_cue = plot(time_device, dist_to_cue,'k', 'LineWidth',2);
        end
        p_dist_to_end = plot(time_device, dist_to_end,'magenta', 'LineWidth',2);
        p_dist_from_start_tgt = plot(time_device, dist_from_start_tgt,'g', 'LineWidth',0.5);
        
        ylabel('pos (deg)');
        xlabel('time (ms)');

        hAx_(2) = subplot(2,1,2);
        hold on;
        ylabel('angle diff (deg)');

        tgt_dir_vector      = [TRIAL.cue_x - TRIAL.start_x, TRIAL.cue_y - TRIAL.start_y];
        unit_tgt_dir_vector = tgt_dir_vector./norm(tgt_dir_vector);
        
        sac_dir_vector_p      = [eye_px_filt - TRIAL.start_x, eye_py_filt - TRIAL.start_y];
        unit_sac_dir_vector_p = sac_dir_vector_p./vecnorm(sac_dir_vector_p')';
        
        ang_diff_p = acosd(unit_sac_dir_vector_p*unit_tgt_dir_vector');
        p_ang_diff_p = plot(time_device,ang_diff_p,'-','Color','r','LineWidth',0.5);
        
        sac_dir_vector_v      = [eye_vx_filt, eye_vy_filt];
        unit_sac_dir_vector_v = sac_dir_vector_v./vecnorm(sac_dir_vector_v')';
        
        ang_diff_v = acosd(unit_sac_dir_vector_v*unit_tgt_dir_vector');
        p_ang_diff_v = plot(time_device,ang_diff_v,'-','Color','b','LineWidth',0.25);
        
        xline(TRIAL.time_state_cue_present,'-','cue present')
        xline(TRIAL.time_state_sac_onset,'-','sac start')

        for counter_sac = 1 : length(SACS_ALL_TRIAL.validity)
            tag_ = SACS_ALL_TRIAL.tag(counter_sac);
            if     tag_==1
                color_ = 'r';
            elseif tag_==2
                color_ = [0.6350    0.0780    0.1840]; % 'r';
            elseif tag_==3
                color_ = [0.4 0.4 0.4];
            elseif tag_==4
                color_ = 'b';
            elseif tag_==5
                color_ = [0    0.4470    0.7410]; % 'b';
            elseif tag_==6
                color_ = [0.4940    0.1840    0.5560]; % 'm';
            elseif tag_==7
                color_ = [0.4940    0.1840    0.5560]; % 'm';
            elseif tag_==8
                color_ = [0.6 0.6 0.6];
            elseif tag_==9
                color_ = [0.6 0.6 0.6];
            elseif tag_==10
                color_ = [0.6 0.6 0.6];
            elseif tag_ == 11
                color_ = [0.4660 0.6740 0.1880]; % greenish
            elseif tag_ == 13
                color_ = [0.3010 0.7450 0.9330];
            elseif ((tag_ == 14)) % for tag 14, ignore
                continue;
            else
                color_ = [0.6 0.6 0.6];
            end
            label_ = [meta_data.sac_tag_list{tag_} '_' num2str(SACS_ALL_TRIAL.count(counter_sac))];
            if ((tag_ == 11) || (tag_== 13))
                if tag_ == 11
                    legend_text_ = 'prim no corr';
                elseif tag_ == 13
                    legend_text_ = 'corr no db corr';
                end
                x_1 = xline(hAx_(1), SACS_ALL_TRIAL.time_offset(counter_sac), '-','LineWidth',2 ,'color', color_);
                legend(x_1,legend_text_);
                legend(hAx_(1),'boxoff');
                x_2 = xline(hAx_(2), SACS_ALL_TRIAL.time_offset(counter_sac), '-','LineWidth',2 ,'color', color_);
                legend(x_2,legend_text_);
                legend(hAx_(2),'boxoff');

            else
                xline(hAx_(1), SACS_ALL_TRIAL.time_onset(counter_sac), '-', label_, 'color', color_, 'interpret', 'none', 'FontSize', 12)
                xline(hAx_(2), SACS_ALL_TRIAL.time_onset(counter_sac), '-', label_, 'color', color_, 'interpret', 'none', 'FontSize', 12)
            end
            
        end
        subplot(2,1,1);
        legend([p_eye_px, p_eye_py, p_dist_to_cue,p_dist_to_end,p_dist_from_start_tgt],{'x','y','dist to cue','dist to end','dist from start'});
        legend('boxoff');

        subplot(2,1,2);
        legend([p_ang_diff_p, p_ang_diff_v],{'p','v'});
        legend('boxoff');
        
        if TRIAL.jump_cond == 1
            sgtitle_ = 'jump';
        else
            sgtitle_ = 'no jump';
        end
        if TRIAL.task_cond == 1
            sgtitle_ = sprintf('%s, forced', sgtitle_);
        else
            sgtitle_ = sprintf('%s, choice', sgtitle_);
        end
        sgtitle(sgtitle_);

        linkaxes(hAx_,'x')
        ESN_Beautify_Plot(hFig_, [16,9], 9)
%             pause;
    end
    %% Save Saccades
    SACS_ALL(counter_valid_trial) = SACS_ALL_TRIAL;
    counter_valid_trial = counter_valid_trial + 1;
    
end

%% Arrange 'SACS_ALL_DATA'
clearvars -except meta_data ...
    trials_data SACS_ALL eye_str event_data;
clearvars('SACS_ALL_DATA'); 
sac_data = struct;
field_names_SACS_ALL_DATA = fieldnames(SACS_ALL);
SACS_ALL_DATA_cell = struct2cell(SACS_ALL);
for counter_fields = 1 : 1 : length(field_names_SACS_ALL_DATA)
    SACS_ALL_DATA_field_cell = SACS_ALL_DATA_cell(counter_fields,:,:);
    SACS_ALL_DATA_field_cell = reshape(SACS_ALL_DATA_field_cell, 1, []);
    SACS_ALL_DATA_field_mat = cell2mat(SACS_ALL_DATA_field_cell);
    sac_data.(field_names_SACS_ALL_DATA{counter_fields}) = SACS_ALL_DATA_field_mat;
end

meta_data.duration = (trials_data.time_end(end) - trials_data.time_start(1));

%% ADD FIXATION VARIABLES
sac_data = MAF_Fix_Sorter(sac_data,trials_data);

%% FIX time_visual with photodiode
if ~isempty(event_data)
    % Loop over saccades and correct time_visual for tags 1,2,3,6,7
    num_sacs = numel(sac_data.tag);
    REFRENCE_TIME = event_data.align_states.BEHAVE_EB_xcorr_time_1K(:); % the initial ind will be drawn from BEHAVE_EB_xcorr_time_1K
    for counter_sac = 1 : num_sacs
        event_time = sac_data.time_visual(1,counter_sac);
        if isnan(event_time)
            % if the event_time is nan, then skip the event.
            continue;
        end
        event_ind = find(REFRENCE_TIME >= event_time, 1, 'first');
        if isempty(event_ind)
            % if the event_ind is empty, then skip the event.
            sac_data.validity(1, counter_sac) = false;
            continue;
        end
        
        tag_ = sac_data.tag(1, counter_sac);
        if (tag_ == 1) || (tag_ == 2) || (tag_ == 3) || (tag_ == 6) || (tag_ == 7)
            % 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3 % 'back_center_success' tag 6 % 'back_center_prim' tag 7
            if isfield(event_data.align_photodiode,'EPHYS_PD_aligned_ind_1K')
                event_ind_converted = event_data.align_photodiode.EPHYS_PD_aligned_ind_1K(event_ind);
                % Re-compute correct reaction time
                visual_ind_converted = event_ind_converted;
                sac_onset_time = sac_data.time_onset(1,counter_sac);
                sac_onset_ind = find(REFRENCE_TIME >= sac_onset_time, 1, 'first');
                sac_onset_ind_converted = event_data.align_states.EPHYS_EB_aligned_ind_1K(sac_onset_ind);
                reaction_time = sac_onset_ind_converted - visual_ind_converted;
            elseif isfield(event_data.align_photodiode,'MONITOR_PD_aligned_ind_1K')
                event_ind_converted = event_data.align_photodiode.MONITOR_PD_aligned_ind_1K(event_ind);
                % Re-compute correct reaction time
                visual_ind_converted = event_ind_converted; 
                sac_onset_time = sac_data.time_onset(1,counter_sac);
                sac_onset_ind = find(REFRENCE_TIME >= sac_onset_time, 1, 'first');
                reaction_time = sac_onset_ind - visual_ind_converted; % 'visual_ind_converted' is in BEHAVE time domain, so no need to convert 'sac_onset_ind'
            end
            if isempty(reaction_time)
                % if the event_ind is empty, then skip the event.
                sac_data.validity(1, counter_sac) = false;
                continue;
            end
            sac_data.reaction(counter_sac) = reaction_time;
            sac_data.time_visual(counter_sac) = sac_data.time_onset(counter_sac) - reaction_time/1e3;
        end
    end
end

end