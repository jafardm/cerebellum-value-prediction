function TRIALS_DATA = JDM_trial_organizer(data, meta_data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:    data: behavior data for  one recording session
%            meta_data: output of MAF_monkey_behavior_sac.m 
%                               includes meta_data including file address,
%                               # trials, tag list, sampling frequency, etc.
%
% Output:   TRIAL_DATA: behavior data specified for all trials of this session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ANALYZE TRIALS
data_fieldnames = fieldnames(data);
num_trials = sum(cell2mat(arrayfun(@(x) strfind(x,'trial_'),data_fieldnames))==1);

% Reference sampling frequency (downsample everything to 1K)
sampling_freq   = 1e3;
for counter_trial = 1 :  num_trials-1
    % Extract Trial Variables
    clearvars -except meta_data ...
        TRIALS_DATA_ALL TRIALS data counter_trial flag_figure sampling_freq name_time_device;

    % Get trial struct
    trial_struct = data.(sprintf('trial_%d',counter_trial));

    % Trial start/end time
    TRIAL.time_start = double(trial_struct.tgt_time_data(1));
    TRIAL.time_end   = double(trial_struct.tgt_time_data(end));

    % Trial variables
    TRIAL.start_x                    = ESN_Round(trial_struct.start_x, 0.01);
    TRIAL.start_y                    = ESN_Round(trial_struct.start_y, 0.01);
    TRIAL.cue_x                      = ESN_Round(trial_struct.cue_x, 0.01);
    TRIAL.cue_y                      = ESN_Round(trial_struct.cue_y, 0.01);
    TRIAL.end_x                      = ESN_Round(trial_struct.end_x, 0.01);
    TRIAL.end_y                      = ESN_Round(trial_struct.end_y, 0.01);
    TRIAL.iss_x                      = ESN_Round(TRIAL.end_x - TRIAL.cue_x, 0.01);
    TRIAL.iss_y                      = ESN_Round(TRIAL.end_y - TRIAL.cue_y, 0.01);
    TRIAL.reward_area                = ESN_Round(data.rew_area, 0.01);
    TRIAL.time_pursuit               = double(data.pursuit_dur);
    if isfield(trial_struct, 'pursuit_x')
        TRIAL.pursuit_x                  = ESN_Round(trial_struct.pursuit_x, 0.01);
        TRIAL.pursuit_y                  = ESN_Round(trial_struct.pursuit_y, 0.01);                              
    else
        TRIAL.pursuit_x                  = nan;
        TRIAL.pursuit_y                  = nan;
    end
    TRIAL.time_state_str_pursuit     = double(trial_struct.state_start_t_str_tgt_pursuit);
    TRIAL.time_state_str_present     = double(trial_struct.state_start_t_str_tgt_present);
    TRIAL.time_state_str_fixation    = double(trial_struct.state_start_t_str_tgt_fixation);
    TRIAL.time_state_cue_present     = double(trial_struct.state_start_t_cue_tgt_present);
    TRIAL.time_state_sac_detect_on   = double(trial_struct.state_start_t_detect_sac_start);
    TRIAL.time_state_sac_onset       = double(trial_struct.state_start_t_saccade);
    TRIAL.time_state_sac_detect_off  = double(trial_struct.state_start_t_detect_sac_end);
    TRIAL.time_state_reward          = double(trial_struct.state_start_t_deliver_rew);
    TRIAL.time_state_end_fixation    = double(trial_struct.state_start_t_end_tgt_fixation);
    TRIAL.time_state_incorrect_sac   = double(trial_struct.state_start_t_incorrect_saccade);
    TRIAL.time_state_trial_success   = double(trial_struct.state_start_t_trial_success);
    TRIAL.time_iti                   = double(data.ITI);
    TRIAL.time_punishment            = double(data.pun_time);
    TRIAL.time_fixation              = double(data.min_fix_time);

    %%%%%%%%%%%%%%%%%  REWARD PARAMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(trial_struct, 'choice'), trial_struct.choice = nan; end
    if ~isfield(trial_struct, 'jump_cond'), trial_struct.jump_cond = nan; end
    if ~isfield(trial_struct, 'high_rew_tgt_num'), trial_struct.high_rew_tgt_num = nan; end
    if ~isfield(trial_struct, 'low_rew_tgt_num'), trial_struct.low_rew_tgt_num = nan; end
    if ~isfield(trial_struct, 'tgt_cond'), trial_struct.tgt_cond= nan; end

    TRIAL.tgt_num                    = trial_struct.tgt_num;            % # tgt. num, under forced condition
    TRIAL.cue_x_high_rew             = trial_struct.cue_x_high_rew;     % # only applicable under choice condition
    TRIAL.cue_x_low_rew              = trial_struct.cue_x_low_rew;      % # only applicable under choice condition
    TRIAL.cue_y_high_rew             = trial_struct.cue_y_high_rew;     % # only applicable under choice condition
    TRIAL.cue_y_low_rew              = trial_struct.cue_y_low_rew;      % # only applicable under choice condition
    TRIAL.high_rew_tgt_num           = trial_struct.high_rew_tgt_num;   % # tgt. num. for high reward, if under high reward condition
    TRIAL.low_rew_tgt_num            = trial_struct.low_rew_tgt_num;    % # tgt. num for low reward, if under low reward condition
    TRIAL.tgt_cond                   = trial_struct.low_rew_tgt_num;    % # tgt. under forced condition

    % choice = 0, force = 1
     if strcmp(strcat(trial_struct.task_cond),'forced')   % # Trial type
       TRIAL.task_cond = 1; 
       TRIAL.choice = nan;
     elseif strcmp(strcat(trial_struct.task_cond),'choice')
        TRIAL.task_cond = 0;  
     end
     
     % condition under forced: high reward = 1, low reward = 0
      if strcmp(strcat(trial_struct.rew_cond),'h')  % # high or low reward condition, under forced condition
       TRIAL.rew_cond =  1;  
     elseif strcmp(strcat(trial_struct.rew_cond),'l')
        TRIAL.rew_cond = 0;  
      end
    % choice: high reward = 1, low reward = 0
     if strcmp(strcat(trial_struct.choice),'h')     % # animal's choice, under choice condition
       TRIAL.choice  = 1;  
     elseif strcmp(strcat(trial_struct.choice),'l')
        TRIAL.choice = 0; 
     else
         TRIAL.choice = nan;
     end

      % high reward target = 1, low reward target= 0
     if strcmp(strcat(trial_struct.tgt_cond),'h')     % # stimulus presented to monkey
       TRIAL.tgt_cond  = 1;  
     elseif strcmp(strcat(trial_struct.tgt_cond),'l')
        TRIAL.tgt_cond = 0; 
     else
         TRIAL.tgt_cond = nan;
     end

    % jump = 1, no jump = 0
     if strcmp(strcat(trial_struct.jump_cond),'y')  % # stimulus jump under reward prediction task
            TRIAL.jump_cond  = 1;  
     elseif strcmp(strcat(trial_struct.jump_cond),'n')
             TRIAL.jump_cond = 0; 
     else 
             TRIAL.jump_cond = nan;  
     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure the sampling frequency is in Hz
    if meta_data.sampling_freq == 1
        meta_data.sampling_freq = 1e3;
    elseif meta_data.sampling_freq == 2
        meta_data.sampling_freq = 2e3;
    end

    % Trial time series
    if isfield(data,'eye_rx_raw_data') % if eye_rx_raw_data exists in data struct, we're in chamber A!
        device_time_data = double(data.device_time_data) * 1e-3;
        if meta_data.sampling_freq == 1e3
            trial_time_start = trial_struct.tgt_time_data(1);
            trial_time_stop = trial_struct.tgt_time_data(end);
            trial_idx_start = find(device_time_data>=trial_time_start,1);
            trial_idx_stop = find(device_time_data>=trial_time_stop,1);
            TRIAL.time_device = device_time_data(trial_idx_start: trial_idx_stop);
            eye_r_px_raw = data.eye_rx_raw_data(trial_idx_start: trial_idx_stop);
            eye_r_py_raw = data.eye_ry_raw_data(trial_idx_start: trial_idx_stop);
            %eye_r_blink  = data.eye_r_blink_data(trial_idx_start: trial_idx_stop);
            eye_l_px_raw = data.eye_lx_raw_data(trial_idx_start: trial_idx_stop);
            eye_l_py_raw = data.eye_ly_raw_data(trial_idx_start: trial_idx_stop);
            %eye_l_blink  = data.eye_l_blink_data(trial_idx_start: trial_idx_stop);
        else
            error('not supported sampling rate for chamber A!');
        end
    elseif isfield(trial_struct,'eye_rx_raw_data') % if eye_rx_raw_data exists in each trial folder separately, we're in chamber B!
        if meta_data.sampling_freq == 2e3
            TRIAL.time_device        = trial_struct.device_time_data(1:2:end);
            eye_r_px_raw = resample(trial_struct.eye_rx_raw_data,1,2);
            eye_r_py_raw = resample(trial_struct.eye_ry_raw_data,1,2);
%             eye_r_blink  = resample(trial_struct.eye_r_blink_data,1,2);
            eye_l_px_raw = resample(trial_struct.eye_lx_raw_data,1,2);
            eye_l_py_raw = resample(trial_struct.eye_ly_raw_data,1,2);
%             eye_l_blink  = resample(trial_struct.eye_l_blink_data,1,2);
        else
            error('not supported sampling rate for chamber B!');
        end
    end

    TRIAL.time_1K = TRIAL.time_device(1) : 0.001 : TRIAL.time_device(end);
    
    invalid_r_idx = isnan(TRIAL.time_device);
    device_r_time = TRIAL.time_device;
    invalid_l_idx = isnan(TRIAL.time_device);
    device_l_time = TRIAL.time_device;
    % Remove invalid indices (mainly from blinking)
    invalid_r_idx = invalid_r_idx | isnan(eye_r_px_raw);
    invalid_r_idx = invalid_r_idx | isnan(eye_r_py_raw);
    % invalid_r_idx = invalid_r_idx | eye_r_blink;
    eye_r_px_raw(invalid_r_idx) = []; 
    eye_r_py_raw(invalid_r_idx) = [];
    device_r_time(invalid_r_idx) = [];
    
    invalid_l_idx = invalid_l_idx | isnan(eye_l_px_raw);
    invalid_l_idx = invalid_l_idx | isnan(eye_l_py_raw);
    % invalid_l_idx = invalid_l_idx | eye_l_blink;
    eye_l_px_raw(invalid_l_idx) = []; 
    eye_l_py_raw(invalid_l_idx) = [];
    device_l_time(invalid_l_idx) = [];
    
    % Check if either eye_r_px_raw, eye_r_py_raw, eye_l_px_raw or eye_l_py_raw was all NaN, something is wrong with this trial, ignore it for now!
    if all(isnan(eye_r_px_raw(:))) || all(isnan(eye_r_py_raw(:))) || all(isnan(eye_l_px_raw(:))) || all(isnan(eye_l_py_raw(:)))
        warning("There is no pupil data, Check the raw data!!!!!!")
        continue;
    end
    
    % Reconstruct eye_r data
    TRIAL.eye_r_px_raw = interp1(device_r_time, eye_r_px_raw, TRIAL.time_1K, 'linear', 'extrap');
    TRIAL.eye_r_py_raw = interp1(device_r_time, eye_r_py_raw, TRIAL.time_1K, 'linear', 'extrap');
   
    % Reconstruct eye_l data
    TRIAL.eye_l_px_raw = interp1(device_l_time, eye_l_px_raw, TRIAL.time_1K, 'linear', 'extrap');
    TRIAL.eye_l_py_raw = interp1(device_l_time, eye_l_py_raw, TRIAL.time_1K, 'linear', 'extrap');
    
    % Reconstruct tgt data
    TRIAL.tgt_px = interp1(trial_struct.tgt_time_data, trial_struct.tgt_x_data, TRIAL.time_1K, 'nearest', 'extrap');
    TRIAL.tgt_py = interp1(trial_struct.tgt_time_data, trial_struct.tgt_y_data, TRIAL.time_1K, 'nearest', 'extrap');

    % Transform to screen coordinates
    r_eye_tracked = strcmp(meta_data.eye_tracked,'right');
    l_eye_tracked = strcmp(meta_data.eye_tracked,'left');


    if isfield(data,'device_time_data') % chamber A
        if r_eye_tracked
            r_cal_matrix = trial_struct.left_cal_matrix; % use left_cal_matrix from raw_data
            if size(r_cal_matrix,1) == 1
                r_cal_matrix = reshape(r_cal_matrix,3,3);
                r_cal_matrix = r_cal_matrix';
            end
        end

    else                                % chamber B
        if r_eye_tracked
            r_cal_matrix = trial_struct.right_cal_matrix; % may change every trial due to real time correction
            if size(r_cal_matrix,1) == 1
                r_cal_matrix = reshape(r_cal_matrix,3,3);
                r_cal_matrix = r_cal_matrix';
            end
        end
        if l_eye_tracked
            l_cal_matrix = trial_struct.left_cal_matrix; % may change every trial due to real time correction
            if size(l_cal_matrix,1) == 1
                l_cal_matrix = reshape(l_cal_matrix,3,3);
                l_cal_matrix = l_cal_matrix';
            end
        end

    end
    
    % filter params
    cutoff_freq = 100.0;
    [b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
    
      % @Nazanin
    if isfield(data,'device_time_data')  % chamber A
        if r_eye_tracked
            eye_p_raw = [TRIAL.eye_l_px_raw',TRIAL.eye_l_py_raw',ones(length(TRIAL.eye_l_py_raw),1)]; % use left eye position in raw data
            eye_p     = eye_p_raw*r_cal_matrix;
            TRIAL.eye_r_px = eye_p(:,1)';
            TRIAL.eye_r_py = eye_p(:,2)';
            TRIAL.eye_r_vx = diff(TRIAL.eye_r_px)./diff(TRIAL.time_1K); 
            TRIAL.eye_r_vx = [TRIAL.eye_r_vx(1) TRIAL.eye_r_vx];
            TRIAL.eye_r_vy = diff(TRIAL.eye_r_py)./diff(TRIAL.time_1K); 
            TRIAL.eye_r_vy = [TRIAL.eye_r_vy(1) TRIAL.eye_r_vy];
            TRIAL.eye_r_vm = sqrt(TRIAL.eye_r_vx.^2 + TRIAL.eye_r_vy.^2);
            
            TRIAL.eye_r_px_filt = filtfilt(b_butter,a_butter,TRIAL.eye_r_px);
            TRIAL.eye_r_py_filt = filtfilt(b_butter,a_butter,TRIAL.eye_r_py);
            TRIAL.eye_r_vx_filt = diff(TRIAL.eye_r_px_filt)./diff(TRIAL.time_1K); 
            TRIAL.eye_r_vx_filt = [TRIAL.eye_r_vx_filt(1) TRIAL.eye_r_vx_filt];
            TRIAL.eye_r_vy_filt = diff(TRIAL.eye_r_py_filt)./diff(TRIAL.time_1K); 
            TRIAL.eye_r_vy_filt = [TRIAL.eye_r_vy_filt(1) TRIAL.eye_r_vy_filt];
            TRIAL.eye_r_vm_filt = sqrt(TRIAL.eye_r_vx_filt.^2 + TRIAL.eye_r_vy_filt.^2);
        end
    else                                % chamber B
        if r_eye_tracked
            eye_p_raw = [TRIAL.eye_r_px_raw',TRIAL.eye_r_py_raw',ones(length(TRIAL.eye_r_py_raw),1)];
            eye_p     = eye_p_raw*r_cal_matrix;
            TRIAL.eye_r_px = eye_p(:,1)';
            TRIAL.eye_r_py = eye_p(:,2)';
            TRIAL.eye_r_vx = diff(TRIAL.eye_r_px)./diff(TRIAL.time_1K); 
            TRIAL.eye_r_vx = [TRIAL.eye_r_vx(1) TRIAL.eye_r_vx];
            TRIAL.eye_r_vy = diff(TRIAL.eye_r_py)./diff(TRIAL.time_1K); 
            TRIAL.eye_r_vy = [TRIAL.eye_r_vy(1) TRIAL.eye_r_vy];
            TRIAL.eye_r_vm = sqrt(TRIAL.eye_r_vx.^2 + TRIAL.eye_r_vy.^2);
            
            TRIAL.eye_r_px_filt = filtfilt(b_butter,a_butter,TRIAL.eye_r_px);
            TRIAL.eye_r_py_filt = filtfilt(b_butter,a_butter,TRIAL.eye_r_py);
            TRIAL.eye_r_vx_filt = diff(TRIAL.eye_r_px_filt)./diff(TRIAL.time_1K); 
            TRIAL.eye_r_vx_filt = [TRIAL.eye_r_vx_filt(1) TRIAL.eye_r_vx_filt];
            TRIAL.eye_r_vy_filt = diff(TRIAL.eye_r_py_filt)./diff(TRIAL.time_1K); 
            TRIAL.eye_r_vy_filt = [TRIAL.eye_r_vy_filt(1) TRIAL.eye_r_vy_filt];
            TRIAL.eye_r_vm_filt = sqrt(TRIAL.eye_r_vx_filt.^2 + TRIAL.eye_r_vy_filt.^2);
        end
        if l_eye_tracked
            eye_p_raw = [TRIAL.eye_l_px_raw',TRIAL.eye_l_py_raw',ones(length(TRIAL.eye_l_py_raw),1)];
            eye_p     = eye_p_raw*l_cal_matrix;
            TRIAL.eye_l_px = eye_p(:,1)';
            TRIAL.eye_l_py = eye_p(:,2)';
            TRIAL.eye_l_vx = diff(TRIAL.eye_l_px)./diff(TRIAL.time_1K); 
            TRIAL.eye_l_vx = [TRIAL.eye_l_vx(1) TRIAL.eye_l_vx];
            TRIAL.eye_l_vy = diff(TRIAL.eye_l_py)./diff(TRIAL.time_1K); 
            TRIAL.eye_l_vy = [TRIAL.eye_l_vy(1) TRIAL.eye_l_vy];
            TRIAL.eye_l_vm = sqrt(TRIAL.eye_l_vx.^2 + TRIAL.eye_l_vy.^2);
            
            TRIAL.eye_l_px_filt = filtfilt(b_butter,a_butter,TRIAL.eye_l_px);
            TRIAL.eye_l_py_filt = filtfilt(b_butter,a_butter,TRIAL.eye_l_py);
            TRIAL.eye_l_vx_filt = diff(TRIAL.eye_l_px_filt)./diff(TRIAL.time_1K); 
            TRIAL.eye_l_vx_filt = [TRIAL.eye_l_vx_filt(1) TRIAL.eye_l_vx_filt];
            TRIAL.eye_l_vy_filt = diff(TRIAL.eye_l_py_filt)./diff(TRIAL.time_1K); 
            TRIAL.eye_l_vy_filt = [TRIAL.eye_l_vy_filt(1) TRIAL.eye_l_vy_filt];
            TRIAL.eye_l_vm_filt = sqrt(TRIAL.eye_l_vx_filt.^2 + TRIAL.eye_l_vy_filt.^2);
        end
    end
    
    % Trial state indices
    TRIAL.time = TRIAL.time_1K;
    TRIAL.ind_state_str_pursuit    = find(TRIAL.time>TRIAL.time_state_str_pursuit(end), 1, 'first');
    TRIAL.ind_state_str_present    = find(TRIAL.time>TRIAL.time_state_str_present(end), 1, 'first');
    TRIAL.ind_state_str_fixation   = find(TRIAL.time>TRIAL.time_state_str_fixation(end), 1, 'first');
    TRIAL.ind_state_cue_present    = find(TRIAL.time>TRIAL.time_state_cue_present(end), 1, 'first');
    TRIAL.ind_state_sac_detect_on  = find(TRIAL.time>TRIAL.time_state_sac_detect_on(end), 1, 'first');
    TRIAL.ind_state_sac_onset      = find(TRIAL.time>TRIAL.time_state_sac_onset(end), 1, 'first');
    TRIAL.ind_state_sac_detect_off = find(TRIAL.time>TRIAL.time_state_sac_detect_off(end), 1, 'first');
    TRIAL.ind_state_reward         = find(TRIAL.time>TRIAL.time_state_reward(end), 1, 'first');
    TRIAL.ind_state_end_fixation   = find(TRIAL.time>TRIAL.time_state_end_fixation(end), 1, 'first');
    TRIAL.ind_state_next_trial     = find(TRIAL.time>TRIAL.time_state_trial_success(end), 1, 'first');
    
    % Compute the start position Bias
    inds_start_fixation = TRIAL.ind_state_cue_present - round(TRIAL.time_fixation * 1000) : 1 : TRIAL.ind_state_cue_present;
    
    if r_eye_tracked
        state_start_eye_r_px_start_fixation = TRIAL.eye_r_px_filt(inds_start_fixation);
        state_start_eye_r_py_start_fixation = TRIAL.eye_r_py_filt(inds_start_fixation);
        start_x_bias = mean(state_start_eye_r_px_start_fixation) - TRIAL.start_x;
        start_y_bias = mean(state_start_eye_r_py_start_fixation) - TRIAL.start_y;
        TRIAL.start_rx_bias = start_x_bias;
        TRIAL.start_ry_bias = start_y_bias;
    end
    if l_eye_tracked
        state_start_eye_l_px_start_fixation = TRIAL.eye_l_px_filt(inds_start_fixation);
        state_start_eye_l_py_start_fixation = TRIAL.eye_l_py_filt(inds_start_fixation);
        start_x_bias = mean(state_start_eye_l_px_start_fixation) - TRIAL.start_x;
        start_y_bias = mean(state_start_eye_l_py_start_fixation) - TRIAL.start_y;
        TRIAL.start_lx_bias = start_x_bias;
        TRIAL.start_ly_bias = start_y_bias;
    end
    
    % Build TRIALS
    TRIALS(counter_trial) = TRIAL;
    
    %% DEBUGGING FSM
    flag_plot_trial = 0;
    if flag_plot_trial
        hFig_ = figure(1);
        clf(hFig_);
        disp(counter_trial);
        
        time_device = ESN_Round(TRIAL.time_1K,0.001);
        if r_eye_tracked
            eye_px_raw = TRIAL.eye_r_px;
            eye_py_raw = TRIAL.eye_r_py;
        elseif l_eye_tracked
            eye_px_raw = TRIAL.eye_l_px;
            eye_py_raw = TRIAL.eye_l_py;
        end
        eye_vx_raw = [0, movmean(diff(eye_px_raw)./diff(time_device),[1,0])];
        eye_vy_raw = [0, movmean(diff(eye_py_raw)./diff(time_device),[1,0])];
        eye_vm_raw = sqrt(eye_vx_raw.^2 + eye_vy_raw.^2);

        tgt_px = TRIAL.tgt_px;
        tgt_py = TRIAL.tgt_py;

        if TRIAL.task_cond == 0
            dist_to_cue_low_rew  = sqrt((TRIAL.cue_x_low_rew - eye_px_raw).^2 + (TRIAL.cue_y_low_rew - eye_py_raw).^2);
            dist_to_cue_high_rew = sqrt((TRIAL.cue_x_high_rew - eye_px_raw).^2 + (TRIAL.cue_y_high_rew - eye_py_raw).^2);
        else
            dist_to_cue = sqrt((TRIAL.cue_x - eye_px_raw).^2 + (TRIAL.cue_y - eye_py_raw).^2);
            dist_to_end = sqrt((TRIAL.end_x - eye_px_raw).^2 + (TRIAL.end_y - eye_py_raw).^2);
        end
        dist_from_start_tgt = sqrt((TRIAL.start_x - eye_px_raw).^2 + (TRIAL.start_y - eye_py_raw).^2);

        hAx_(1) = subplot(2,1,1);
        hold on;
        p_eye_px = plot(time_device, eye_px_raw, '-','Color','r','LineWidth',1);
        p_eye_py = plot(time_device, eye_py_raw, '-','Color','b','LineWidth',1);
        if TRIAL.task_cond == 0 % choice trial
            state_times  = [TRIAL.time_state_cue_present, TRIAL.time_state_str_pursuit, TRIAL.time_state_incorrect_sac, TRIAL.time_state_end_fixation];
            state_values = [ones(size(TRIAL.time_state_cue_present)).*1, ones(size(TRIAL.time_state_str_pursuit)).*2, ones(size(TRIAL.time_state_incorrect_sac)).*3, ones(size(TRIAL.time_state_end_fixation)).*4];
            [sorted_state_times, sorted_state_times_idx] = sort(state_times,'ascend');
            sorted_state_values = state_values(sorted_state_times_idx);
            % Find period of low and high reward target presentation
            num_cue_state = length(TRIAL.time_state_cue_present);
            idx_time_tgt  = zeros(size(time_device));
            for counter_sac_state = 1 : num_cue_state
                cue_present_time = TRIAL.time_state_cue_present(counter_sac_state);
                sorted_cue_state_idx = find(sorted_state_times >= cue_present_time,1,'first');
                if counter_sac_state == num_cue_state
                    idx_time_tgt = idx_time_tgt | ((time_device >= ESN_Round(cue_present_time,0.001)) & (time_device <= ESN_Round(TRIAL.time_state_end_fixation,0.001)));
                else
                    idx_time_tgt = idx_time_tgt | ((time_device >= ESN_Round(cue_present_time,0.001)) & (time_device <= ESN_Round(sorted_state_times(sorted_cue_state_idx+1),0.001)));
                end
            end
            tgt_px_low_rew  = tgt_px;
            tgt_py_low_rew  = tgt_py;
            tgt_px_high_rew = tgt_px;
            tgt_py_high_rew = tgt_py;
            tgt_px_low_rew(idx_time_tgt)  = TRIAL.cue_x_low_rew;
            tgt_py_low_rew(idx_time_tgt)  = TRIAL.cue_y_low_rew;
            tgt_px_high_rew(idx_time_tgt) = TRIAL.cue_x_high_rew;
            tgt_py_high_rew(idx_time_tgt) = TRIAL.cue_y_high_rew;
            p_tgt_px_low_rew  = plot(time_device, tgt_px_low_rew, '--', 'Color','r');
            p_tgt_py_low_rew = plot(time_device, tgt_py_low_rew, '--', 'Color','b');
            p_tgt_px_high_rew  = plot(time_device, tgt_px_high_rew, '--', 'Color','r');
            p_tgt_py_high_rew = plot(time_device, tgt_py_high_rew, '--', 'Color','b');
            p_dist_to_cue_low_rew  = plot(time_device, dist_to_cue_low_rew,'k', 'LineWidth',2);
            p_dist_to_cue_high_rew = plot(time_device, dist_to_cue_high_rew,'k', 'LineWidth',2);
        else
            p_tgt_px = plot(time_device, tgt_px, '--', 'Color','r');
            p_tgt_py = plot(time_device, tgt_py, '--', 'Color','b');
            p_dist_to_cue = plot(time_device, dist_to_cue,'k', 'LineWidth',2);
            if TRIAL.jump_cond == 1
                p_dist_to_end = plot(time_device, dist_to_end,'magenta', 'LineWidth',2);
            end
        end
        p_dist_from_start_tgt = plot(time_device, dist_from_start_tgt,'g', 'LineWidth',0.5);
        p_rew_dist = yline(TRIAL.reward_area/2,'-','rew dist','Color',[0.4940 0.1840 0.5560]);
        
        xline(TRIAL.time_state_str_pursuit,'-','pursuit');
        xline(TRIAL.time_state_cue_present,'-','cue present');
%         xline(TRIAL.time_state_sac_detect_on,'-','detect sac');
        xline(TRIAL.time_state_sac_onset,'-','sac');
        if ~isempty(TRIAL.time_state_incorrect_sac)
            xline(TRIAL.time_state_incorrect_sac,'-','incorrect sac');
        end

        ylabel('pos (deg)');
        xlabel('time (ms)');
        
        if TRIAL.task_cond == 0
            legend([p_eye_px,p_eye_py,p_dist_to_cue_low_rew,p_dist_to_cue_high_rew,p_dist_from_start_tgt],{'x','y','dist to cue low rew','dist to cue high rew','dist to start tgt'});
        else
            if TRIAL.jump_cond == 1
                legend([p_eye_px,p_eye_py,p_dist_to_cue,p_dist_to_end,p_dist_from_start_tgt],{'x','y','dist to cue','dist to end','dist to start tgt'})
            else
                legend([p_eye_px,p_eye_py,p_dist_to_cue,p_dist_from_start_tgt],{'x','y','dist to cue','dist to start tgt'})
            end
        end
        legend('boxoff');

        hAx_(2) = subplot(2,1,2);
        hold on;
        ylabel('angle diff (deg)');
   
        if TRIAL.task_cond == 0
            tgt_low_dir_vector      = [TRIAL.cue_x_low_rew - TRIAL.start_x, TRIAL.cue_y_low_rew - TRIAL.start_y];
            unit_tgt_low_rew_dir_vector = reshape(tgt_low_dir_vector./norm(tgt_low_dir_vector),2,[]);
            tgt_high_dir_vector      = [TRIAL.cue_x_high_rew - TRIAL.start_x, TRIAL.cue_y_high_rew - TRIAL.start_y];
            unit_tgt_high_rew_dir_vector = reshape(tgt_high_dir_vector./norm(tgt_high_dir_vector),2,[]);
        else
            tgt_dir_vector      = [TRIAL.cue_x - TRIAL.start_x, TRIAL.cue_y - TRIAL.start_y];
            unit_tgt_dir_vector = reshape(tgt_dir_vector./norm(tgt_dir_vector),2,[]);
        end
        sac_dir_vector_v      = [eye_vx_raw; eye_vy_raw];
        unit_sac_dir_vector_v = sac_dir_vector_v./vecnorm(sac_dir_vector_v);
        unit_sac_dir_vector_v(isnan(unit_sac_dir_vector_v)) = 0;

        if TRIAL.task_cond == 0
            ang_diff_v_low_rew = acosd(unit_tgt_low_rew_dir_vector'*unit_sac_dir_vector_v);
            ang_diff_v_high_rew = acosd(unit_tgt_high_rew_dir_vector'*unit_sac_dir_vector_v);
        else
            ang_diff_v = acosd(unit_tgt_dir_vector'*unit_sac_dir_vector_v);
        end

        num_sac_state = length(TRIAL.time_state_sac_onset);
        idx_time_sac  = zeros(size(time_device));
        for counter_sac_state = 1 : num_sac_state
            if counter_sac_state == num_sac_state
                idx_time_sac = idx_time_sac | ((time_device >= ESN_Round(TRIAL.time_state_sac_onset(counter_sac_state),0.001)) & (time_device <= ESN_Round(TRIAL.time_state_sac_detect_off(end),0.001)));
            else
                if ~isempty(TRIAL.time_state_incorrect_sac) && (length(TRIAL.time_state_incorrect_sac) >= counter_sac_state)
                    idx_time_sac = idx_time_sac | ((time_device >= ESN_Round(TRIAL.time_state_sac_onset(counter_sac_state),0.001)) & (time_device <= ESN_Round(TRIAL.time_state_incorrect_sac(counter_sac_state),0.001)));
                end
            end
        end
        
        if TRIAL.task_cond == 0
            p_ang_diff_v_low_rew = plot(time_device(idx_time_sac),ang_diff_v_low_rew(idx_time_sac),'o','Color','b','MarkerSize',8);
            p_ang_diff_v_high_rew = plot(time_device(idx_time_sac),ang_diff_v_high_rew(idx_time_sac),'o','Color','r','MarkerSize',8);
        else
            p_ang_diff_v = plot(time_device(idx_time_sac),ang_diff_v(idx_time_sac),'o','Color','b','MarkerSize',8);
        end
        p_ang_threshold = yline(45,'-','ang threshold','Color',[0.4940 0.1840 0.5560]);
        
        if TRIAL.jump_cond == 1
            sgtitle_ = 'jump';
        else
            sgtitle_ = 'no jump';
        end
        sgtitle(sgtitle_);

        linkaxes(hAx_,'x')
        ESN_Beautify_Plot(hFig_, [16,9], 9)
    end
end
%% TRIAL_DATA
clearvars -except meta_data TRIALS_DATA_ALL TRIALS_DATA TRIALS flag_figure;
clearvars('TRIALS_DATA'); TRIALS_DATA = struct;
field_names_TRIALS = fieldnames(TRIALS);
for counter_fields = 1 : 1 : length(field_names_TRIALS)
    for counter_trials = 1 : 1 : length(TRIALS)
        variable_TRIALS_ = TRIALS(counter_trials).(field_names_TRIALS{counter_fields});
        % handling an error which the variable_TRIALS_ was []
        if isempty(variable_TRIALS_)
            variable_TRIALS_ = nan;
        end
        variable_TRIALS_ = variable_TRIALS_(:);
        if max(size(variable_TRIALS_)) > 1
            variable_TRIALS_ = mat2cell(variable_TRIALS_, size(variable_TRIALS_,1), size(variable_TRIALS_,2));
        end
        % the field does not exist in TRIALS_DATA
        if ~isfield(TRIALS_DATA, field_names_TRIALS{counter_fields})
            TRIALS_DATA.(field_names_TRIALS{counter_fields}) = [];
        end
        variable_TRIALS_DATA_ = TRIALS_DATA.(field_names_TRIALS{counter_fields});
        % variable_TRIALS_ is cell array
        if iscell(variable_TRIALS_)
            % variable_TRIALS_DATA_ is cell array
            % variables are compatible (both are cell), add new data
            if iscell(variable_TRIALS_DATA_)
                variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
            end
            % variable_TRIALS_DATA_ is matrix array
            % convert variable_TRIALS_DATA_ to cell, and add new data
            if isnumeric(variable_TRIALS_DATA_)
                variable_TRIALS_DATA_ = num2cell(variable_TRIALS_DATA_);
                variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
            end   
        end
        % variable_TRIALS_ is matrix array
        if isnumeric(variable_TRIALS_)
            % variable_TRIALS_DATA_ is cell array
            % convert variable_TRIALS_ to cell, and add new data
            if iscell(variable_TRIALS_DATA_)
                variable_TRIALS_ = mat2cell(variable_TRIALS_, size(variable_TRIALS_,1), size(variable_TRIALS_,2));
                variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
            end
            % variable_TRIALS_DATA_ is matrix array
            % variables are compatible (both are matrix), add new data
            if isnumeric(variable_TRIALS_DATA_)
                variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
            end
        end
        TRIALS_DATA.(field_names_TRIALS{counter_fields}) = variable_TRIALS_DATA_;
        
    end
end

% Build TRIALS_DATA_ALL
TRIALS_DATA_ALL = TRIALS_DATA;

% Arrange 'TRIALS_DATA'
clearvars -except meta_data TRIALS_DATA_ALL TRIALS_DATA flag_figure;
clearvars('TRIALS_DATA'); TRIALS_DATA = struct;
field_names_TRIALS_DATA_ALL = fieldnames(TRIALS_DATA_ALL);
for counter_fields = 1 : 1 : length(field_names_TRIALS_DATA_ALL)
    for counter_files = 1 : 1 : length(TRIALS_DATA_ALL)
        variable_TRIALS_DATA_ALL_ = TRIALS_DATA_ALL(counter_files).(field_names_TRIALS_DATA_ALL{counter_fields});
        % the field does not exist in TRIALS_DATA
        if ~isfield(TRIALS_DATA, field_names_TRIALS_DATA_ALL{counter_fields})
            TRIALS_DATA.(field_names_TRIALS_DATA_ALL{counter_fields}) = [];
        end
        variable_TRIALS_DATA_ = TRIALS_DATA.(field_names_TRIALS_DATA_ALL{counter_fields});
        variable_TRIALS_DATA_ = horzcat(variable_TRIALS_DATA_, variable_TRIALS_DATA_ALL_);
        TRIALS_DATA.(field_names_TRIALS_DATA_ALL{counter_fields}) = variable_TRIALS_DATA_;
    end
end

end