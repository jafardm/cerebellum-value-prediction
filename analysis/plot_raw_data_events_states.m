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
path_data_monkey_data = 'Z:\Ephys\data_65F\';
month                 = '2025-01';
session               = '2025-01-01';
rec_dir = dir2(fullfile(path_data_monkey_data,month,session));
rec_dir = rec_dir([rec_dir.isdir],:);

for idx_d = 1:size(rec_dir,1)
    data_files = dir2(fullfile(rec_dir(idx_d).folder,rec_dir(idx_d).name));
    
    analyzed_dir = dir2(fullfile(data_files(1).folder,data_files(3).name, 'behavior_data\eye\','*ANALYZED*'));
    raw_dir      = dir2(fullfile(data_files(2).folder,data_files(4).name,'*.mat*'));

    sacdata = load(fullfile(analyzed_dir.folder,analyzed_dir.name));
    data = load(fullfile(raw_dir.folder,raw_dir.name));
    data = data.data;
    % eyelink timestamps
    device_time_data = double(data.device_time_data) * 1e-3;

    data_fieldnames = fieldnames(data);
    num_trials = sum(arrayfun(@(str) contains(str,'trial'),data_fieldnames));

    for counter_trial = 1:num_trials
    
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
    trial_data.eye_px_filt    = eye_px_filt  - data.horz_offset;
    trial_data.eye_py_filt    = eye_py_filt  - data.vert_offset;
    trial_data.eye_vx_filt    = eye_vx_filt;
    trial_data.eye_vy_filt    = eye_vy_filt;
    trial_data.eye_vm_filt    = eye_vm_filt;
    trial_data.start_x        = data.horz_offset; % to be fixed; save directly from experiment
    trial_data.start_y        = data.vert_offset;
    trial_data.tgt_px_filt    = tgt_px - data.horz_offset;
    trial_data.tgt_py_filt    = tgt_py - data.vert_offset;

    % Plot trial
    fig = figure('Units','normalized','Position',[0,0,1,1]);
    clf;
    sgtitle(sprintf('trial %d',counter_trial));
    ax = gca;
    hold on;
    ylabel('pos (deg)');
    xlabel('time (ms)');
    
    time_plot = trial_data.time_1K.*1000;
    time_offset = (trial_data.time_1K(1))*1000;
    time_plot = time_plot - time_offset;
    % eye_blink = trial_data.eye_blink_filt;
    xlim([time_plot(1),time_plot(end)]);
    p_eye_px = plot(time_plot, trial_data.eye_px_filt, '-','Color','r','LineWidth',2);
    p_eye_py = plot(time_plot, trial_data.eye_py_filt, '-','Color','b','LineWidth',2);
    p_tgt_px = plot(time_plot, trial_data.tgt_px_filt, '--', 'Color','r','LineWidth',2);
    p_tgt_py = plot(time_plot, trial_data.tgt_py_filt, '--', 'Color','b','LineWidth',2);
   


    JDM_plot_state_time(trial_data.state_start_t_str_tgt_fixation.*1000,time_offset,'Fixation','k');
    JDM_plot_state_time(trial_data.state_start_t_cue_tgt_present.*1000,time_offset,'Primary target','k');
    JDM_plot_state_time(trial_data.state_start_t_saccade.*1000,time_offset,'Saccade','k');
    JDM_plot_state_time(trial_data.state_start_t_end_tgt_fixation.*1000,time_offset,'End target','k'); 
    JDM_plot_state_time(trial_data.state_start_t_trial_success.*1000,time_offset,'Trial succes','k');
    
    legend([p_eye_px, p_eye_py,p_tgt_px,p_tgt_py],{'x','y','tgt x','tgt y'});
    legend('boxoff');
    ESN_Beautify_Plot(fig,[20,8]);
    pause;
     
   close(gcf);
    end
end
%%
function JDM_plot_state_time(state_time, time_offset,state_name, color)
    state_time = state_time - time_offset;
    for i = 1 : length(state_time)
        xline(state_time(i), '-', state_name, ...
    'Color', color, ...
    'LabelHorizontalAlignment', 'center', ...
    'LabelVerticalAlignment', 'middle', ...
    'LineWidth', 2); % Add the LineWidth property

    end
end

