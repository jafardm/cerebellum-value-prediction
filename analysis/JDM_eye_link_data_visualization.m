
%% SACCADE PARAMETER
% Eye pos. filter parameter
sampling_freq = 1000.0;
cutoff_freq = 100.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
% Primary saccade parameter
params_prim.MinPeakHeight       = 75.0; % deg/s
params_prim.MinPeakProminence   = 50; % data points
params_prim.rough_threshold     = 50.0; % deg/s
params_prim.fine_threshold      = 20.0; % deg/s
params_prim.sampling_freq       = 1000.0; % Hz
params_prim.cutoff_freq         = 50.0; % Hz
params_prim.window_half_length  = 4; % data points
params_prim.prominence_or_first = 'prominent'; % which peak to select, 'prominent' or 'first'
% Corrective saccade parameter
params_corr.MinPeakHeight       = 50.0; % deg/s
params_corr.MinPeakProminence   = 50; % data points
params_corr.rough_threshold     = 40.0; % deg/s
params_corr.fine_threshold      = 20.0; % deg/s
params_corr.sampling_freq       = 1000.0; % Hz
params_corr.cutoff_freq         = 75.0; % Hz
params_corr.window_half_length  = 4; % data points
params_corr.prominence_or_first = 'prominent'; % which peak to select, 'prominent' or 'first'
threshold_pos = 3.5; % in deg., to determine validity of saccade depending on position of eye wrt. target of interest
%% LOAD DATA
path_data_monkey_data = 'Y:\\Ephys\data_132F\';
month                 = '2024-03';
session               = '2024-03-01';
rec_dir = dir2(fullfile(path_data_monkey_data,month,session));
rec_dir = rec_dir([rec_dir.isdir],:);

for idx_d = 1:size(rec_dir,1)
    data_files = dir2(fullfile(rec_dir(idx_d).folder,rec_dir(idx_d).name));
    
    analyzed_dir = dir2(fullfile(data_files(1).folder,data_files(1).name, 'behavior_data\eye\','*ANALYZED*'));
    raw_dir      = dir2(fullfile(data_files(2).folder,data_files(2).name,'*.mat*'));

    sacdata = load(fullfile(analyzed_dir.folder,analyzed_dir.name));
    data = load(fullfile(raw_dir.folder,raw_dir.name));
    data = data.data;
    % eyelink timestamps
    device_time_data = double(data.device_time_data) * 1e-3;

    data_fieldnames = fieldnames(data);
    num_trials = sum(arrayfun(@(str) contains(str,'trial'),data_fieldnames));

    for counter_trial = 9:num_trials
    
    disp(counter_trial)
    trial_data  = data.(sprintf('trial_%d',counter_trial));
  
    % The eye x-y coordination for current trial
    start_t = trial_data.tgt_time_data(1);
    stop_t = trial_data.tgt_time_data(end);
    start_idx = find(device_time_data>=start_t,1);
    stop_idx = find(device_time_data>=stop_t,1);
    eye_x = data.eye_lx_raw_data(start_idx: stop_idx);
    eye_y = data.eye_ly_raw_data(start_idx: stop_idx);
    
    % calibration matrix for current trial
    caL_matrix = trial_data.left_cal_matrix;  
    
       if size(caL_matrix,1) == 1
            if length(caL_matrix) == 18
                caL_matrix = caL_matrix(1:9);
            end
            caL_matrix = transpose(reshape(caL_matrix,3,3));
       end

    time_trial = device_time_data(start_idx:stop_idx)';
    
    eye_p_raw = [eye_x',eye_y',ones(length(eye_x),1)];
    eye_p     = eye_p_raw*caL_matrix;
    
    % Transform to screen coordinates 
    eye_px = eye_p(:,1);
    eye_py = eye_p(:,2);
    
    invalid_idx = isnan(time_trial);
    % Remove invalid indices (mainly from blinking)
    invalid_idx = invalid_idx | isnan(eye_px);
    invalid_idx = invalid_idx | isnan(eye_py);
  
    eye_px(invalid_idx) = []; 
    eye_py(invalid_idx) = [];
    time_trial(invalid_idx) = [];
%     eye_px = interp1(time_target, x_target, double(time_trial), 'previous');
%     eye_py = interp1(time_target, y_target, double(time_trial), 'previous');
    
    
    % Filter eye data
    eye_px_filt = filtfilt(b_butter,a_butter,eye_px);
    eye_py_filt = filtfilt(b_butter,a_butter,eye_py);
    
    % Compute eye speed
    eye_vx_filt = diff(eye_px_filt)./diff(time_trial); 
    eye_vx_filt = [eye_vx_filt(1); eye_vx_filt]; % make it the same length as position data
    eye_vy_filt = diff(eye_py_filt)./diff(time_trial); 
    eye_vy_filt = [eye_vy_filt(1); eye_vy_filt];
    eye_vm_filt = sqrt(eye_vx_filt.^2+eye_vy_filt.^2);
    
    % Get real time tgt data and resample @ 1000 Hz 
    time_real   = trial_data.tgt_time_data; % time registered during FSM loop
    tgt_px_real = trial_data.tgt_x_data;
    tgt_py_real = trial_data.tgt_y_data;
    
    tgt_px      = zeros(size(time_trial));
    tgt_py      = zeros(size(time_trial));
    
    pre_idx_1K = 1;
    for counter_t = 1 : length(time_real)
        t = time_real(counter_t);
        
        pst_idx_1k = find(time_trial >= t, 1, 'first');
        
        tgt_px(pre_idx_1K:pst_idx_1k) = tgt_px_real(counter_t);
        tgt_py(pre_idx_1K:pst_idx_1k) = tgt_py_real(counter_t);
        
        pre_idx_1K = pst_idx_1k;
    end
    
    % Save to trial struct.
    trial_data.time_1K        = time_trial;
    trial_data.eye_px_filt    = eye_px_filt;
    trial_data.eye_py_filt    = eye_py_filt;
    trial_data.eye_vx_filt    = eye_vx_filt;
    trial_data.eye_vy_filt    = eye_vy_filt;
    trial_data.eye_vm_filt    = eye_vm_filt;
    trial_data.start_x        = data.horz_offset; % to be fixed; save directly from experiment
    trial_data.start_y        = data.vert_offset;
    trial_data.tgt_px_filt    = tgt_px;
    trial_data.tgt_py_filt    = tgt_py;

    tags_ = 1;
    idx_sac = find(sacdata.sac_data.trial_num == counter_trial & ismember(sacdata.sac_data.tag,tags_));

    tags_times_onset = [(sacdata.sac_data.time_onset(idx_sac) - sacdata.trials_data.time_start(counter_trial))]*1000;
    tags_times_offset = [(sacdata.sac_data.time_offset(idx_sac) - sacdata.trials_data.time_start(counter_trial))]*1000;
    
    % Find primary saccade
    prim_sac = find_prim_sac(trial_data, threshold_pos, params_prim);
    
    % Plot trial
    figure;
    clf;
    sgtitle(sprintf('trial %d',counter_trial));
    ax_1 = subplot(4,1,1);
    ax = gca;
    hold on;
    title('processed');
    ylabel('pos (deg)');
    xlabel('time (ms)');
    
    time_plot = trial_data.time_1K.*1000;
    time_offset = trial_data.time_1K(1)*1000;
    time_plot = time_plot - time_offset;
    % eye_blink = trial_data.eye_blink_filt;
    xlim([time_plot(1),time_plot(end)]);
    p_eye_px = plot(time_plot, trial_data.eye_px_filt, '-','Color','r','LineWidth',1);
    p_eye_py = plot(time_plot, trial_data.eye_py_filt, '-','Color','b','LineWidth',1);
    p_tgt_px = plot(time_plot, trial_data.tgt_px_filt, '--', 'Color','r');
    p_tgt_py = plot(time_plot, trial_data.tgt_py_filt, '--', 'Color','b');
    
    dist_to_tgt = sqrt((trial_data.tgt_px_filt - trial_data.eye_px_filt).^2 + (trial_data.tgt_py_filt - trial_data.eye_py_filt).^2);
    dist_from_start_tgt = sqrt((trial_data.start_x - trial_data.eye_px_filt).^2 + (trial_data.start_y - trial_data.eye_py_filt).^2);
    
    p_dist_to_tgt = plot(time_plot, dist_to_tgt,'k', 'LineWidth',2);
    p_dist_from_start_tgt = plot(time_plot, dist_from_start_tgt,'g', 'LineWidth',0.5);
    
    plot_state_time(trial_data.state_start_t_cue_tgt_present.*1000,time_offset,[],'k');
    plot_state_time(trial_data.state_start_t_incorrect_saccade.*1000,time_offset,'incorrect sac','r');
    plot_state_time(trial_data.state_start_t_saccade.*1000,time_offset,'sac','k');
    
    if prim_sac.validity
        prim_sac_color = 'g';
    else
        prim_sac_color = 'r';
    end
    xline(time_plot(prim_sac.ind_start),prim_sac_color, 'LineWidth',1.5);
    xline(time_plot(prim_sac.ind_finish),prim_sac_color,'LineWidth',1.5);
    
    legend([p_eye_px, p_eye_py,p_tgt_px,p_tgt_py, p_dist_to_tgt,p_dist_from_start_tgt],{'x','y','dist to tgt','dist from start'});
    legend('boxoff');
    %%%%%%
    ax_2 = subplot(4,1,2);
    hold on;
    xlim([time_plot(1),time_plot(end)]);
    ylim([min(trial_data.eye_vm_filt),max(trial_data.eye_vm_filt)]);
    ylabel('vel (deg/s)');
    
    p_eye_vm = plot(time_plot, trial_data.eye_vm_filt, '-','Color','r','LineWidth',1.5);
    
    sac_detect_threshold = data.sac_detect_threshold;
    sac_on_off_threshold = data.sac_on_off_threshold;
    yline(sac_detect_threshold,'--','sac detect');
    yline(sac_on_off_threshold,'--','sac on off');
    
    eye_px_real = trial_data.eye_x_data;
    eye_py_real = trial_data.eye_y_data;
    blink_idx   = trial_data.eye_x_data >= 1000;
    eye_px_real(blink_idx) = nan;
    eye_py_real(blink_idx) = nan;
    
    time_plot_real = trial_data.tgt_time_data.*1000;
    time_plot_real = time_plot_real - trial_data.time_1K(1)*1000;
    num_vel_samp = 3;
    
    eye_vx_real  = zeros(length(time_plot_real),1);
    eye_vy_real  = zeros(length(time_plot_real),1);
    for counter_t = num_vel_samp : length(time_plot_real)
        idx = counter_t - num_vel_samp + 1 : counter_t;
        
        eye_vx_real(counter_t) = mean(diff(eye_px_real(idx))./diff(trial_data.tgt_time_data(idx)));
        eye_vy_real(counter_t) = mean(diff(eye_py_real(idx))./diff(trial_data.tgt_time_data(idx)));
    end
    eye_vm_real   = sqrt(eye_vx_real.^2 + eye_vy_real.^2);
    p_eye_vm_real = plot(time_plot_real, eye_vm_real, '-','Color','b','LineWidth',0.5);
    
    lambda     = 0.7;
    eye_vx_new = nan(length(time_plot_real),1);
    eye_vy_new = nan(length(time_plot_real),1);
    eye_vx_new(2) = (eye_px_real(2) - eye_px_real(1))/(trial_data.tgt_time_data(2) - trial_data.tgt_time_data(1));
    eye_vy_new(2) = (eye_py_real(2) - eye_py_real(1))/(trial_data.tgt_time_data(2) - trial_data.tgt_time_data(1));
    for counter_t = 3 : length(time_plot_real)
        eye_vx_new(counter_t) = lambda*eye_vx_new(counter_t-1) + (1-lambda)*(eye_px_real(counter_t) - eye_px_real(counter_t-1))/(trial_data.tgt_time_data(counter_t) - trial_data.tgt_time_data(counter_t-1));
        eye_vy_new(counter_t) = lambda*eye_vy_new(counter_t-1) + (1-lambda)*(eye_py_real(counter_t) - eye_py_real(counter_t-1))/(trial_data.tgt_time_data(counter_t) - trial_data.tgt_time_data(counter_t-1));
    end
    eye_vm_new   = sqrt(eye_vx_new.^2 + eye_vy_new.^2);
    p_eye_vm_new = plot(time_plot_real, eye_vm_new, '-','Color','g','LineWidth',0.5);
    
    plot_state_time(trial_data.state_start_t_cue_tgt_present.*1000,time_offset,[],'k');
    plot_state_time(trial_data.state_start_t_incorrect_saccade.*1000,time_offset,'incorrect sac','r');
    plot_state_time(trial_data.state_start_t_saccade.*1000,time_offset,'sac','k');
    
    
    legend([p_eye_vm, p_eye_vm_real,p_eye_vm_new],{'post','real','new'});
    legend('boxoff');
    %%%%%%
    ax_3 = subplot(4,1,3);
    hold on;
    title('raw');
    ylabel('pos (deg)');
    xlim([time_plot(1),time_plot(end)]);
    ylim(ax.YLim);
    
    plot(time_plot_real, eye_px_real, '-','Color','r','LineWidth',1);
    plot(time_plot_real, eye_py_real, '-','Color','b','LineWidth',1);
    plot(time_plot_real, trial_data.tgt_x_data, '--', 'Color','r');
    plot(time_plot_real, trial_data.tgt_y_data, '--', 'Color','b');
    
    dist_to_tgt = sqrt((trial_data.tgt_x_data - eye_px_real).^2 + (trial_data.tgt_y_data - eye_py_real).^2);
    plot(time_plot_real, dist_to_tgt,'k', 'LineWidth',2);
    
    plot_state_time(trial_data.state_start_t_cue_tgt_present.*1000,time_offset,[],'k');
    plot_state_time(trial_data.state_start_t_incorrect_saccade.*1000,time_offset,'incorrect sac','r');
    plot_state_time(trial_data.state_start_t_saccade.*1000,time_offset,'sac','k');
    %%%%%%
    ax_4 = subplot(4,1,4);
    hold on;
    ylabel('angle diff (deg)');
    
    tgt_dir_vector      = [trial_data.cue_x - trial_data.start_x; trial_data.cue_y - trial_data.start_y];
    unit_tgt_dir_vector = tgt_dir_vector./norm(tgt_dir_vector);
    
    sac_dir_vector_p      = [eye_px_real - trial_data.start_x; eye_py_real - trial_data.start_y];
    unit_sac_dir_vector_p = sac_dir_vector_p./vecnorm(sac_dir_vector_p')';
    
    ang_diff_p = acosd(unit_sac_dir_vector_p'.*unit_tgt_dir_vector');
    p_ang_diff_p = plot(time_plot_real,ang_diff_p,'-','Color','r','LineWidth',1);
    
    sac_dir_vector_v      = [eye_vx_real, eye_vy_real];
    unit_sac_dir_vector_v = sac_dir_vector_v./vecnorm(sac_dir_vector_v);
    
    ang_diff_v = acosd(unit_sac_dir_vector_v*unit_tgt_dir_vector);
    p_ang_diff_v = plot(time_plot_real,ang_diff_v,'-','Color','b','LineWidth',1);
    
    plot_state_time(trial_data.state_start_t_cue_tgt_present.*1000,time_offset,[],'k');
    plot_state_time(trial_data.state_start_t_incorrect_saccade.*1000,time_offset,'incorrect sac','r');
    plot_state_time(trial_data.state_start_t_saccade.*1000,time_offset,'sac','k');
    if isempty(tags_times_onset)
        title('%s','There is no tag 1')
    else
        xline(tags_times_onset,'Color','g','LineWidth',1.5)
        xline(tags_times_offset,'Color','g','LineWidth',1.5)
    end
    % legend([p_ang_diff_p, p_ang_diff_v],{'p','v'});
    legend('boxoff');
    
    linkaxes([ax_1,ax_2,ax_3,ax_4],'x');
    close all
    end
end

%% FIND PRIMARY SACCADE
function prim_sac = find_prim_sac(trial_data, threshold_pos, params_prim)
    time_start = trial_data.state_start_t_cue_tgt_present;
    time_end   = trial_data.state_start_t_deliver_rew;
    idx_start  = find(trial_data.time_1K >= time_start(end),1,'first');
    idx_end    = find(trial_data.time_1K >= time_end(end),1,'first');
    
    idx_search_start = idx_start;
    idx_search_end   = idx_end + 50; % increase search range
    
    sac = ESN_Sac_Finder(trial_data.eye_vm_filt,idx_search_start, idx_search_end, params_prim);
    
    % Validation
    prim_sac.validity   = sac.validity;
    prim_sac.inds      = sac.inds;
    prim_sac.ind_start  = sac.ind_start;
    prim_sac.ind_vmax   = sac.ind_vmax;
    prim_sac.ind_finish = sac.ind_finish;
    prim_sac.eye_px_start = trial_data.eye_px_filt(prim_sac.ind_start);
    prim_sac.eye_py_start = trial_data.eye_py_filt(prim_sac.ind_start);
    prim_sac.eye_px_finish = trial_data.eye_px_filt(prim_sac.ind_finish);
    prim_sac.eye_py_finish = trial_data.eye_py_filt(prim_sac.ind_finish);   
    prim_sac.eye_amp_x = prim_sac.eye_px_finish - prim_sac.eye_px_start;
    prim_sac.eye_amp_y = prim_sac.eye_py_finish - prim_sac.eye_py_start;
    prim_sac.eye_amp_m = sqrt(prim_sac.eye_amp_x.^2 + prim_sac.eye_amp_y.^2);
    prim_sac.eye_vm_max = trial_data.eye_vm_filt(prim_sac.ind_vmax);

    % Distance btwn eye start position and start target
    prim_sac.diff_start = sqrt((prim_sac.eye_px_start-trial_data.start_x).^2 +(prim_sac.eye_py_start-trial_data.start_y).^2);
    % Dist. btwn. eye finish pos. and cue tgt.
    prim_sac.diff_finish = sqrt((prim_sac.eye_px_finish-trial_data.cue_x).^2 +(prim_sac.eye_py_finish-trial_data.cue_y).^2);
    % If the saccade doesn't meet req., invalidiate
    if ((prim_sac.diff_start >= threshold_pos) || (prim_sac.diff_finish >= threshold_pos))
        prim_sac.validity = 0;
    end
end
%% FIND CORR SACCADE
function corr_sac = find_corr_sac(trial_data, threshold_pos, params_corr)
    time_start = trial_data.state_start_t_detect_sac_end;
    time_end   = trial_data.state_start_t_trial_success;
    idx_start  = find(trial_data.time_2K >= time_start(end),1,'first');
    idx_end    = find(trial_data.time_2K >= time_end(end),1,'first');
    
    idx_search_start = idx_start + 100;
    idx_search_end   = idx_end + 50; % increase search range
    
    sac = ESN_Sac_Finder(trial_data.eye_vm_filt,idx_search_start, idx_search_end, params_corr);
    
    % Validation
    corr_sac.validity   = sac.validity;
    corr_sac.inds      = sac.inds;
    corr_sac.ind_start  = sac.ind_start;
    corr_sac.ind_vmax   = sac.ind_vmax;
    corr_sac.ind_finish = sac.ind_finish;
    corr_sac.eye_px_start = trial_data.eye_px_filt(corr_sac.ind_start);
    corr_sac.eye_py_start = trial_data.eye_py_filt(corr_sac.ind_start);
    corr_sac.eye_px_finish = trial_data.eye_px_filt(corr_sac.ind_finish);
    corr_sac.eye_py_finish = trial_data.eye_py_filt(corr_sac.ind_finish);   
    corr_sac.eye_amp_x = corr_sac.eye_px_finish - corr_sac.eye_px_start;
    corr_sac.eye_amp_y = corr_sac.eye_py_finish - corr_sac.eye_py_start;
    corr_sac.eye_amp_m = sqrt(corr_sac.eye_amp_x.^2 + corr_sac.eye_amp_y.^2);
    corr_sac.eye_vm_max = trial_data.eye_vm_filt(corr_sac.ind_vmax);

    % Distance btwn eye start position and start target
    corr_sac.diff_start = sqrt((corr_sac.eye_px_start-trial_data.cue_x).^2 +(corr_sac.eye_py_start-trial_data.cue_y).^2);
    % Dist. btwn. eye finish pos. and cue tgt.
    corr_sac.diff_finish = sqrt((corr_sac.eye_px_finish-trial_data.end_x).^2 +(corr_sac.eye_py_finish-trial_data.end_y).^2);
    % If the saccade doesn't meet req., invalidiate
    if ((corr_sac.diff_start >= threshold_pos) || (corr_sac.diff_finish >= threshold_pos))
        corr_sac.validity = 0;
    end
end
%% PLOT STATE
function plot_state_time(state_time, time_offset,state_name, color)
    state_time = state_time - time_offset;
    for i = 1 : length(state_time)
        xline(state_time(i),'-',state_name,'Color',color,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
    end
end