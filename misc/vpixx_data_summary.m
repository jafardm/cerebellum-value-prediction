%% SACCADE PARAMETER
% Eye pos. filter parameter
sampling_freq = 2000;
cutoff_freq = 100;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
% Primary saccade parameter
params_prim.MinPeakHeight       = 75.0; % deg/s
params_prim.MinPeakProminence   = 50; % data points
params_prim.rough_threshold     = 50.0; % deg/s
params_prim.fine_threshold      = 20.0; % deg/s
params_prim.sampling_freq       = 2000.0; % Hz
params_prim.cutoff_freq         = 50.0; % Hz
params_prim.window_half_length  = 4; % data points
params_prim.prominence_or_first = 'prominent'; % which peak to select, 'prominent' or 'first'
% Corrective saccade parameter
params_corr.MinPeakHeight       = 50.0; % deg/s
params_corr.MinPeakProminence   = 50; % data points
params_corr.rough_threshold     = 40.0; % deg/s
params_corr.fine_threshold      = 20.0; % deg/s
params_corr.sampling_freq       = 2000.0; % Hz
params_corr.cutoff_freq         = 75.0; % Hz
params_corr.window_half_length  = 4; % data points
params_corr.prominence_or_first = 'prominent'; % which peak to select, 'prominent' or 'first'
threshold_pos = 3.5; % in deg., to determine validity of saccade depending on position of eye wrt. target of interest
%% LOAD DATA
loaded_data = load('Z:\Ephys\data_132F\2023-02\2023-02-02\2023-02-02_11-56-33\raw_data\random_corrective_saccades_115643.mat');
data = loaded_data.data;
%% PLOT TRIAL (SIMPLE SACCADE)
trial_idx = 2;
which_eye = 'l'; % 'l' - left; 'r' - right

figure;

trial_data = data.(sprintf('trial_%d',trial_idx));

cal_matrix = trial_data.cal_matrix; % may change every trial due to real time correction
cal_matrix = reshape(cal_matrix,3,3);
cal_matrix = cal_matrix';
vpixx_time = trial_data.vpixx_time_data;
time_2K    = vpixx_time(1):0.0005:vpixx_time(end);

invalid_idx = isnan(vpixx_time);

eye_px_raw = trial_data.(sprintf('eye_%sx_raw_data',which_eye));
eye_py_raw = trial_data.(sprintf('eye_%sy_raw_data',which_eye));
eye_blink  = trial_data.(sprintf('eye_%s_blink_data',which_eye));

% Remove invalid indices (mainly from blinking)
invalid_idx = invalid_idx | isnan(eye_px_raw);
invalid_idx = invalid_idx | isnan(eye_py_raw);
invalid_idx = invalid_idx | eye_blink;
eye_px_raw(invalid_idx) = []; 
eye_py_raw(invalid_idx) = [];
vpixx_time(invalid_idx) = [];

% Resampling at 2000 Hz
eye_px_raw = interp1(vpixx_time, eye_px_raw, time_2K,'linear','extrap');
eye_py_raw = interp1(vpixx_time, eye_py_raw, time_2K, 'linear','extrap');
eye_p_raw  = [eye_px_raw', eye_py_raw', ones(length(eye_px_raw),1)];

% Transform to screen coordinates 
eye_p  = eye_p_raw*cal_matrix;
eye_px = eye_p(:,1);
eye_py = eye_p(:,2);

% Filter eye data
eye_px_filt = filtfilt(b_butter,a_butter,eye_px);
eye_py_filt = filtfilt(b_butter,a_butter,eye_py);

% Get real time tgt data and resample @ 2000 Hz 
time_real   = trial_data.tgt_time_data; % time registered during FSM loop
tgt_px_real = trial_data.tgt_x_data;
tgt_py_real = trial_data.tgt_y_data;
tgt_px      = zeros(size(time_2K));
tgt_py      = zeros(size(time_2K));

pre_idx_2K = 1;
for counter_t = 1 : length(time_real)
    t = time_real(counter_t);
    
    post_idx_2K = find(time_2K >= t, 1, 'first');
    
    tgt_px(pre_idx_2K:post_idx_2K) = tgt_px_real(counter_t);
    tgt_py(pre_idx_2K:post_idx_2K) = tgt_py_real(counter_t);
    
    pre_idx_2K = post_idx_2K;
end

% Compute eye speed
eye_vx_filt = diff(eye_px_filt)./diff(time_2K'); 
eye_vx_filt = [eye_vx_filt(1); eye_vx_filt]; % make it the same length as position data
eye_vy_filt = diff(eye_py_filt)./diff(time_2K'); 
eye_vy_filt = [eye_vy_filt(1); eye_vy_filt];
eye_vm_filt = sqrt(eye_vx_filt.^2+eye_vy_filt.^2);

% Save to trial struct.
trial_data.time_2K        = time_2K;
trial_data.eye_px_filt    = eye_px_filt;
trial_data.eye_py_filt    = eye_py_filt;
trial_data.eye_vx_filt    = eye_vx_filt;
trial_data.eye_vy_filt    = eye_vy_filt;
trial_data.eye_vm_filt    = eye_vm_filt;
trial_data.start_x        = data.horz_offset; % to be fixed; save directly from experiment
trial_data.start_y        = data.vert_offset;
trial_data.tgt_px_filt    = tgt_px;
trial_data.tgt_py_filt    = tgt_py;
trial_data.eye_blink_filt = eye_blink;
% Reshape data
data_length      = length(time_2K);
data_length_real = length(time_real);
trial_data_field = fieldnames(trial_data);
for counter_field = 1 : length(trial_data_field)
    field_name = trial_data_field{counter_field};
    field_data = trial_data.(field_name); 
    if ((size(field_data,1) == data_length) || (size(field_data,2) == data_length))
        trial_data.(field_name) = reshape(field_data, data_length,1);     
    elseif ((size(field_data,1) == data_length_real) || (size(field_data,2) == data_length_real))
        trial_data.(field_name) = reshape(field_data, data_length_real,1);      
    end
end


% Find primary saccade
prim_sac = find_prim_sac(trial_data, threshold_pos, params_prim);

% Plot trial
sgtitle(sprintf('trial %d',trial_idx));
ax_1 = subplot(4,1,1);
ax = gca;
hold on;
title('processed');
ylabel('pos (deg)');
xlabel('time (ms)');

time_plot = trial_data.time_2K.*1000;
time_offset = trial_data.time_2K(1)*1000;
time_plot = time_plot - time_offset;
eye_blink = trial_data.eye_blink_filt;
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

legend([p_eye_px, p_eye_py, p_dist_to_tgt,p_dist_from_start_tgt],{'x','y','dist to tgt','dist from start'});
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
time_plot_real = time_plot_real - trial_data.time_2K(1)*1000;
num_vel_samp = 3;
eye_vx_real  = nan(length(time_plot_real),1);
eye_vy_real  = nan(length(time_plot_real),1);
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
unit_sac_dir_vector_v = sac_dir_vector_v./vecnorm(sac_dir_vector_v')';

ang_diff_v = acosd(unit_sac_dir_vector_v*unit_tgt_dir_vector);
p_ang_diff_v = plot(time_plot_real,ang_diff_v,'-','Color','b','LineWidth',1);

plot_state_time(trial_data.state_start_t_cue_tgt_present.*1000,time_offset,[],'k');
plot_state_time(trial_data.state_start_t_incorrect_saccade.*1000,time_offset,'incorrect sac','r');
plot_state_time(trial_data.state_start_t_saccade.*1000,time_offset,'sac','k');

legend([p_ang_diff_p, p_ang_diff_v],{'p','v'});
legend('boxoff');

linkaxes([ax_1,ax_2,ax_3,ax_4],'x');
%% PLOT TRIAL (CORR SACCADE) 
trial_idx = 3;
which_eye = 'l'; % 'l' - left; 'r' - right

figure;

trial_data = data.(sprintf('trial_%d',trial_idx));

cal_matrix = trial_data.left_cal_matrix; % may change every trial due to real time correction
if size(cal_matrix,1) == 1
    cal_matrix = reshape(cal_matrix,3,3);
    cal_matrix = cal_matrix';
end
vpixx_time = trial_data.device_time_data;
time_2K    = vpixx_time(1):0.0005:vpixx_time(end);

invalid_idx = isnan(vpixx_time);

eye_px_raw = trial_data.(sprintf('eye_%sx_raw_data',which_eye));
eye_py_raw = trial_data.(sprintf('eye_%sy_raw_data',which_eye));
eye_blink  = trial_data.(sprintf('eye_%s_blink_data',which_eye));

% Remove invalid indices (mainly from blinking)
invalid_idx = invalid_idx | isnan(eye_px_raw);
invalid_idx = invalid_idx | isnan(eye_py_raw);
invalid_idx = invalid_idx | eye_blink;
eye_px_raw(invalid_idx) = []; 
eye_py_raw(invalid_idx) = [];
vpixx_time(invalid_idx) = [];

% Resampling at 2000 Hz
eye_px_raw = interp1(vpixx_time, eye_px_raw, time_2K,'linear','extrap');
eye_py_raw = interp1(vpixx_time, eye_py_raw, time_2K, 'linear','extrap');
eye_p_raw  = [eye_px_raw', eye_py_raw', ones(length(eye_px_raw),1)];

% Transform to screen coordinates 
eye_p  = eye_p_raw*cal_matrix;
eye_px = eye_p(:,1);
eye_py = eye_p(:,2);

% Filter eye data
eye_px_filt = filtfilt(b_butter,a_butter,eye_px);
eye_py_filt = filtfilt(b_butter,a_butter,eye_py);

% Get real time tgt data and resample @ 2000 Hz 
time_real   = trial_data.tgt_time_data; % time registered during FSM loop
tgt_px_real = trial_data.tgt_x_data;
tgt_py_real = trial_data.tgt_y_data;
tgt_px      = zeros(size(time_2K));
tgt_py      = zeros(size(time_2K));
pre_idx_2K = 1;
for counter_t = 1 : length(time_real)
    t = time_real(counter_t);
    
    post_idx_2K = find(time_2K >= t, 1, 'first');
    
    tgt_px(pre_idx_2K:post_idx_2K) = tgt_px_real(counter_t);
    tgt_py(pre_idx_2K:post_idx_2K) = tgt_py_real(counter_t);
    
    pre_idx_2K = post_idx_2K;
end

% Compute eye speed
eye_vx_filt = diff(eye_px_filt)./diff(time_2K'); 
eye_vx_filt = [eye_vx_filt(1); eye_vx_filt]; % make it the same length as position data
eye_vy_filt = diff(eye_py_filt)./diff(time_2K'); 
eye_vy_filt = [eye_vy_filt(1); eye_vy_filt];
eye_vm_filt = sqrt(eye_vx_filt.^2+eye_vy_filt.^2);

% Save to trial struct.
trial_data.time_2K        = time_2K;
trial_data.eye_px_filt    = eye_px_filt;
trial_data.eye_py_filt    = eye_py_filt;
trial_data.eye_vx_filt    = eye_vx_filt;
trial_data.eye_vy_filt    = eye_vy_filt;
trial_data.eye_vm_filt    = eye_vm_filt;
trial_data.start_x        = data.horz_offset; % to be fixed; save directly from experiment
trial_data.start_y        = data.vert_offset;
trial_data.tgt_px_filt    = tgt_px;
trial_data.tgt_py_filt    = tgt_py;
trial_data.eye_blink_filt = eye_blink;
% Reshape data
data_length      = length(time_2K);
data_length_real = length(time_real);
trial_data_field = fieldnames(trial_data);
for counter_field = 1 : length(trial_data_field)
    field_name = trial_data_field{counter_field};
    field_data = trial_data.(field_name); 
    if ((size(field_data,1) == data_length) || (size(field_data,2) == data_length))
        trial_data.(field_name) = reshape(field_data, data_length,1);     
    elseif ((size(field_data,1) == data_length_real) || (size(field_data,2) == data_length_real))
        trial_data.(field_name) = reshape(field_data, data_length_real,1);      
    end
end


% Find saccades
prim_sac = find_prim_sac(trial_data, threshold_pos, params_prim);
corr_sac = find_corr_sac(trial_data, threshold_pos, params_corr);

% Plot trial
sgtitle(sprintf('trial %d',trial_idx));
ax_1 = subplot(4,1,1);
ax = gca;
hold on;
title('processed');
ylabel('pos (deg)');
xlabel('time (ms)');

time_plot = trial_data.time_2K.*1000;
time_offset = trial_data.time_2K(1)*1000;
time_plot = time_plot - time_offset;
eye_blink = trial_data.eye_blink_filt;
xlim([time_plot(1),time_plot(end)]);
p_eye_px = plot(time_plot, trial_data.eye_px_filt, '-','Color','r','LineWidth',1);
p_eye_py = plot(time_plot, trial_data.eye_py_filt, '-','Color','b','LineWidth',1);
p_tgt_px = plot(time_plot, trial_data.tgt_px_filt, '--', 'Color','r');
p_tgt_py = plot(time_plot, trial_data.tgt_py_filt, '--', 'Color','b');

dist_to_cue = sqrt((trial_data.cue_x - trial_data.eye_px_filt).^2 + (trial_data.cue_y - trial_data.eye_py_filt).^2);
dist_to_end = sqrt((trial_data.end_x - trial_data.eye_px_filt).^2 + (trial_data.end_y - trial_data.eye_py_filt).^2);
dist_from_start_tgt = sqrt((trial_data.start_x - trial_data.eye_px_filt).^2 + (trial_data.start_y - trial_data.eye_py_filt).^2);

p_dist_to_cue = plot(time_plot, dist_to_cue,'k', 'LineWidth',2);
p_dist_to_end = plot(time_plot, dist_to_end,'magenta', 'LineWidth',2);
p_dist_from_start_tgt = plot(time_plot, dist_from_start_tgt,'g', 'LineWidth',0.5);

plot_state_time(trial_data.state_start_t_deliver_rew.*1000,time_offset,'delive rew','k');
plot_state_time(trial_data.state_start_t_incorrect_saccade.*1000,time_offset,'incorrect sac','r');
plot_state_time(trial_data.state_start_t_saccade.*1000,time_offset,'sac','k');
% plot_state_time(trial_data.state_start_t_detect_sac_end.*1000,time_offset,'end','k');
% plot_state_time(trial_data.state_start_t_trial_success.*1000,time_offset,'success','k');

if prim_sac.validity
    prim_sac_color = 'g';
else
    prim_sac_color = 'r';
end
xline(time_plot(prim_sac.ind_start),prim_sac_color, 'LineWidth',1.5);
xline(time_plot(prim_sac.ind_finish),prim_sac_color,'LineWidth',1.5);
if corr_sac.validity
    corr_sac_color = 'g';
else
    corr_sac_color = 'r';
end
xline(time_plot(corr_sac.ind_start),corr_sac_color, 'LineWidth',1.5);
xline(time_plot(corr_sac.ind_finish),corr_sac_color,'LineWidth',1.5);

legend([p_eye_px, p_eye_py, p_dist_to_cue,p_dist_to_end,p_dist_from_start_tgt],{'x','y','dist to cue','dist to end','dist from start'});
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
time_plot_real = time_plot_real - trial_data.time_2K(1)*1000;
num_vel_samp = 3;
eye_vx_real  = nan(length(time_plot_real),1);
eye_vy_real  = nan(length(time_plot_real),1);
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
dist_to_cue = sqrt((trial_data.cue_x - eye_px_real).^2 + (trial_data.cue_y - eye_py_real).^2);
plot(time_plot_real, dist_to_tgt,'k', 'LineWidth',2);
plot(time_plot_real, dist_to_cue,'r', 'LineWidth',2);

plot_state_time(trial_data.state_start_t_cue_tgt_present.*1000,time_offset,[],'k');
plot_state_time(trial_data.state_start_t_incorrect_saccade.*1000,time_offset,'incorrect sac','r');
plot_state_time(trial_data.state_start_t_saccade.*1000,time_offset,'sac','k');
%%%%%%
ax_4 = subplot(4,1,4);
hold on;
ylabel('angle diff (deg)');

tgt_dir_vector      = [trial_data.cue_x - trial_data.start_x, trial_data.cue_y - trial_data.start_y];
unit_tgt_dir_vector = tgt_dir_vector./norm(tgt_dir_vector);

sac_dir_vector_p      = [eye_px_real - trial_data.start_x, eye_py_real - trial_data.start_y];
unit_sac_dir_vector_p = sac_dir_vector_p./vecnorm(sac_dir_vector_p')';

ang_diff_p = acosd(unit_sac_dir_vector_p*unit_tgt_dir_vector');
p_ang_diff_p = plot(time_plot_real,ang_diff_p,'-','Color','r','LineWidth',1);

sac_dir_vector_v      = [eye_vx_real, eye_vy_real];
unit_sac_dir_vector_v = sac_dir_vector_v./vecnorm(sac_dir_vector_v')';

ang_diff_v = acosd(unit_sac_dir_vector_v*unit_tgt_dir_vector');
p_ang_diff_v = plot(time_plot_real,ang_diff_v,'-','Color','b','LineWidth',1);

plot_state_time(trial_data.state_start_t_cue_tgt_present.*1000,time_offset,[],'k');
plot_state_time(trial_data.state_start_t_incorrect_saccade.*1000,time_offset,'incorrect sac','r');
plot_state_time(trial_data.state_start_t_saccade.*1000,time_offset,'sac','k');

legend([p_ang_diff_p, p_ang_diff_v],{'p','v'});
legend('boxoff');
linkaxes([ax_1,ax_2,ax_3,ax_4],'x')
%% PD DELAY, DATAPIXX SIDE
trial_idx = 3;
pd_ch = 1;
d_ch  = 1;
tgt_on_change  = -1;
tgt_off_change = 1;

trial_data = data.(sprintf('trial_%d',trial_idx));
vpixx_time = trial_data.vpixx_time_data;

% Tgt-on delay (more important)
pd_on_change_idc   = find(diff(trial_data.din_data) == tgt_on_change);
dout_on_change_idc = find(diff(trial_data.dout_data) == tgt_on_change);
fprintf('\n## Tgt on ##\n')
fprintf('Num. dout: %d, num. pd: %d; num. missing signals: %d\n', length(dout_on_change_idc), length(pd_on_change_idc),length(dout_on_change_idc)-length(pd_on_change_idc));
max_delay = 0.100; % if more than this, invalid (s) 
counter_invalid_delay = 0;
counter_delay         = 0;
on_screen_delays      = [];
for counter_change = 1 : length(dout_on_change_idc)
    dout_on_change_idx = dout_on_change_idc(counter_change);
    [~,min_pd_on_change_idx] = min(abs(pd_on_change_idc - dout_on_change_idx));
    pd_on_change_idx         = pd_on_change_idc(min_pd_on_change_idx);
    screen_delay             = vpixx_time(pd_on_change_idx) - vpixx_time(dout_on_change_idx);
    % Screen delay can't be negative nor too long
    if (screen_delay < 0) || (screen_delay > max_delay) || isempty(pd_on_change_idx)
        counter_invalid_delay                =  counter_invalid_delay + 1;
    else
        counter_delay = counter_delay + 1;
        on_screen_delays(counter_delay) = vpixx_time(pd_on_change_idx) - vpixx_time(dout_on_change_idx);
    end
end
fprintf('# invalid delays: %d\n',counter_invalid_delay);
mean_delay_msec = mean(on_screen_delays).*1000;
std_delay_msec  = std(on_screen_delays).*1000;
fprintf('Delay (ms), mean: %1.2f; std: %1.2f\n',mean_delay_msec,std_delay_msec);
% Tgt-off delay (less important)   
pd_off_change_idc   = find(diff(trial_data.din_data) == tgt_off_change);
dout_off_change_idc = find(diff(trial_data.dout_data) == tgt_off_change);
fprintf('## Tgt off ##\n')
fprintf('# dout: %d, # pd: %d; # missing signals: %d\n', length(dout_off_change_idc), length(pd_off_change_idc),length(dout_off_change_idc)-length(pd_off_change_idc));
max_delay = 0.100; % if more than this, invalid (s) 
counter_invalid_delay = 0;
counter_delay         = 0;
off_screen_delays     = [];
for counter_change = 1 : length(dout_off_change_idc)
    dout_on_change_idx = dout_off_change_idc(counter_change);
    [~,min_pd_on_change_idx] = min(abs(pd_off_change_idc - dout_on_change_idx));
    pd_off_change_idx        = pd_off_change_idc(min_pd_on_change_idx);
    screen_delay             = vpixx_time(pd_off_change_idx) - vpixx_time(dout_on_change_idx);
    % Screen delay can't be negative
    if screen_delay < 0
        counter_invalid_delay                =  counter_invalid_delay + 1;
    elseif screen_delay > max_delay
        counter_invalid_delay = counter_invalid_delay + 1;
    else
        counter_delay = counter_delay + 1;
        off_screen_delays(counter_delay) = vpixx_time(pd_off_change_idx) - vpixx_time(dout_on_change_idx);
    end
end
fprintf('# invalid delays: %d\n',counter_invalid_delay);
mean_delay_msec = mean(off_screen_delays).*1000;
std_delay_msec  = std(off_screen_delays).*1000;
fprintf('Delay (ms), mean: %1.2f; std: %1.2f\n',mean_delay_msec,std_delay_msec)
%% STATE AND DIO/PD CHANGES
trial_idx = 5;
pd_ch = 1;
d_ch  = 1;
tgt_on_change  = -1;
tgt_off_change = 1;

figure;
% State
sgtitle(sprintf('trial %d',trial_idx));
ax_1 = subplot(2,1,1);
hold on;
xlabel('time (ms)');

time_offset = trial_data.time_2K(1)*1000;

plot_state_time(trial_data.state_start_t_str_tgt_pursuit.*1000,time_offset,'pursuit',[0,0,0]);
plot_state_time(trial_data.state_start_t_str_tgt_present.*1000,time_offset,'str present',[0 0.4470 0.7410]);
% plot_state_time(trial_data.state_start_t_str_tgt_fixation.*1000,time_offset,'str fixation',[0.8500 0.3250 0.0980]);
plot_state_time(trial_data.state_start_t_cue_tgt_present.*1000,time_offset,'cue',[0.9290 0.6940 0.1250]);
plot_state_time(trial_data.state_start_t_saccade.*1000,time_offset,'sac',[0.4940 0.1840 0.5560]);
plot_state_time(trial_data.state_start_t_deliver_rew.*1000,time_offset,'rew',[0.4660 0.6740 0.1880]);
plot_state_time(trial_data.state_start_t_trial_success.*1000,time_offset,'success',[0.3010 0.7450 0.9330]);
plot_state_time(trial_data.state_start_t_incorrect_saccade.*1000,time_offset,'incorrect',[0.6350 0.0780 0.1840]);

% DIO & PD
pd_on_change_idc    = find(diff(trial_data.din_data) == tgt_on_change);
dout_on_change_idc  = find(diff(trial_data.dout_data) == tgt_on_change);
pd_off_change_idc   = find(diff(trial_data.din_data) == tgt_off_change);
dout_off_change_idc = find(diff(trial_data.dout_data) == tgt_off_change);

pd_on_change_times    = trial_data.vpixx_time_data(pd_on_change_idc);
dout_on_change_times  = trial_data.vpixx_time_data(dout_on_change_idc);
pd_off_change_times   = trial_data.vpixx_time_data(pd_off_change_idc);
dout_off_change_times = trial_data.vpixx_time_data(dout_off_change_idc);

ax_2 = subplot(2,1,2);
hold on;

plot_state_time(dout_on_change_times.*1000,time_offset,'d on',[0,0,0]);
plot_state_time(pd_on_change_times.*1000,time_offset,'pd on',[0.6350 0.0780 0.1840]);
plot_state_time(dout_off_change_times.*1000,time_offset,'d off',[0,0,0]);
plot_state_time(pd_off_change_times.*1000,time_offset,'pd off',[0.6350 0.0780 0.1840]);

linkaxes([ax_1,ax_2],'x')

%% (OLD) SCREEN DELAY
cmd_t_data = data.trial_1.state_start_t_str_tgt_pursuit;  
din_data          = data.trial_1.din_data;
diff_din_data     = diff(din_data);

din_data_diff_idx = find((diff_din_data) == -1);

vpixx_time_data   = data.trial_1.vpixx_time_data;
din_data_diff_t   =  vpixx_time_data(din_data_diff_idx);

clearvars screen_delay
% for counter_stim = 1 : length(cmd_t_data)-1
%     cmd_t   = cmd_t_data(counter_stim);
%     cmd_idx = find(vpixx_time_data >= cmd_t, 1);
%     
%     screen_idx = din_data_diff_idx(find(din_data_diff_idx >= cmd_idx,1));
%     screen_t   = vpixx_time_data(screen_idx); 
%     screen_delay(counter_stim) = screen_t - cmd_t;
% end

% Sometimes PD signal not detected so fewer in number than number of stimulus presentations,
% so better to compute screens based on PD signal, ignoring those missing cases
for counter_pd = 2 : length(din_data_diff_t)
    pd_t = din_data_diff_t(counter_pd);
    
    % Find the latest stimulus presentation time
    prior_stim_t = cmd_t_data(cmd_t_data < pd_t);
    late_stim_t = prior_stim_t(end);
    screen_delay(counter_pd) = pd_t - late_stim_t;
end

for counter_stim = 1 : length(cmd_t_data)
    xline(cmd_t_data(counter_stim),'-r');
    hold on;
end
for counter_pd = 1 : length(din_data_diff_t)
    xline(din_data_diff_t(counter_pd),'-b');
    hold on;
end   

ind_cmd = zeros(size(vpixx_time_data));
for counter_stim = 1 : length(cmd_t_data)
    cmd_t = cmd_t_data(counter_stim);
    ind_cmd(find(vpixx_time_data>cmd_t,1,'first')) = 1;
end
hold on;
plot(ind_cmd,'r');
%% (OLD) HISTOGRAM OF CUE PRESENT. TIMES
% Currently bug in the code; maybe the times should come bf/after target.draw(), not just target.flip()
cue_present_time = [];
for t = 1:num_trial
    trial_field_name = sprintf('trial_%d',t);
    trial_data = data.trial_data.(trial_field_name);
    num_cue_present = length(trial_data.time_cue_present_before); % within a trial, cue could be presented multiple times if animal doesn't go thru with the trial
    for c_idx = 1:num_cue_present
        cue_present_time = [cue_present_time, trial_data.time_cue_present_after(c_idx) - trial_data.time_cue_present_before(c_idx)];
    end
end
% TIME MARKER INTERVALS
figure
time_marker_diff = [];
tgt_x_diff = [];
tgt_y_diff = [];
for t = 1:num_trial
    trial_field_name = sprintf('trial_%d',t);
    trial_data = data.trial_data.(trial_field_name);
    time_marker_diff = [time_marker_diff, diff(trial_data.time_data)];
    tgt_x_diff = [tgt_x_diff, diff(trial_data.tgt_x_data)];
    tgt_y_diff = [tgt_y_diff, diff(trial_data.tgt_y_data)];
end

time_marker_freq = 1./time_marker_diff; % ideally, there shouldn't be any frequencies that fall around 60 Hz, indicating the screen isn't updating as fast as it should.
Histogram
hist_time_marker_freq_axes = subplot(3,2,[1,2]);
hist_time_marker_freq = histogram(hist_time_marker_freq_axes,time_marker_freq);
hist_time_marker_freq.NumBins = 200;
hist_time_marker_freq.Normalization = 'probability';
ylabel('probability'); xlabel('Frequency of finite state machine loop (Hz)'); 
% Time marker freq. vs. change in distance of target
tgt_diff = sqrt(tgt_x_diff.^2 + tgt_y_diff.^2);
figure; hold on;
plot(tgt_diff,'.');
yyaxis right
plot(time_marker_freq,'or');
%% (OLD) PROCESS TRIAL DATA
clearvars -except data cal_matrix num_trial
% Eye pos. filter parameter
sampling_freq = 2000;
cutoff_freq = 100;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
% Primary saccade parameter
params_prim.MinPeakHeight       = 100.0; % deg/s
params_prim.MinPeakProminence   = 40; % data points
params_prim.rough_threshold     = 20.0; % deg/s
params_prim.fine_threshold      = 10.0; % deg/s
params_prim.sampling_freq       = 2000.0; % Hz
params_prim.cutoff_freq         = 50.0; % Hz
params_prim.window_half_length  = 4; % data points
params_prim.prominence_or_first = 'prominent'; % which peak to select, 'prominent' or 'first'
threshold_pos = 3.5; % in deg., to determine validity of saccade depending on position of eye wrt. target of interest
for t = 1:num_trial
    trial_field_name = sprintf('trial_%d',t);
    trial_data = data.trial_data.(trial_field_name);
    vpixx_time = trial_data.vpixx_time; 
    time_2K = vpixx_time(1):0.0005:vpixx_time(end);
    invalid_idc = isnan(vpixx_time);
    % Eye position
    eye_px_raw = trial_data.left_horizontal_eye_raw;   invalid_idc = invalid_idc | isnan(eye_px_raw);
    eye_py_raw = trial_data.left_vertical_eye_raw;   invalid_idc = invalid_idc | isnan(eye_py_raw);
    % Removing invalid indices 
    eye_px_raw(invalid_idc) = []; 
    eye_py_raw(invalid_idc) = [];
    vpixx_time(invalid_idc) = [];
    % Resampling at 2000 Hz
    tgt_x_data = interp1(trial_data.time_data,trial_data.tgt_x_data,time_2K,'linear','extrap'); 
    tgt_y_data = interp1(trial_data.time_data,trial_data.tgt_y_data,time_2K,'linear','extrap');
    eye_px_raw = interp1(vpixx_time, eye_px_raw, time_2K,'linear','extrap');
    eye_py_raw = interp1(vpixx_time, eye_py_raw, time_2K, 'linear','extrap');
    eye_p_raw = [eye_px_raw',eye_py_raw',ones(length(eye_px_raw),1)];
    % Transform to screen coordinates and filter
    eye_p = eye_p_raw*cal_matrix;
    eye_px = eye_p(:,1);
    eye_px_filt = filtfilt(b_butter,a_butter,eye_px);
    eye_py = eye_p(:,2);
    eye_py_filt = filtfilt(b_butter,a_butter,eye_py);
    % Eye speed
    vpixx_time_diff = diff(vpixx_time);
    eye_vx_filt = diff(eye_px_filt)./diff(time_2K'); eye_vx_filt = [eye_vx_filt(1); eye_vx_filt];
    eye_vy_filt = diff(eye_py_filt)./diff(time_2K'); eye_vy_filt = [eye_vy_filt(1); eye_vy_filt];
    eye_vm_filt = sqrt(eye_vx_filt.^2+eye_vy_filt.^2);
    % Primary saccade
    time_cue_present_after = trial_data.time_cue_present_after;
    idx_cue_present_after = find(time_2K >= time_cue_present_after(end),1,'first');
    time_saccade_start = trial_data.time_saccade_start;
    time_saccade_end = trial_data.time_saccade_end;
    vpixx_idx_search_start = idx_cue_present_after;
    vpixx_idx_search_end = find(time_2K >= time_saccade_end(end),1,'first');
    vpixx_idx_search_end = vpixx_idx_search_end + 50;
    output_ = ESN_Sac_Finder(eye_vm_filt,vpixx_idx_search_start, vpixx_idx_search_end, params_prim);
    SAC_PRIM_DATA.validity   = output_.validity;
    SAC_PRIM_DATA.inds      = output_.inds;
    SAC_PRIM_DATA.ind_start  = output_.ind_start;
    SAC_PRIM_DATA.ind_vmax   = output_.ind_vmax;
    SAC_PRIM_DATA.ind_finish = output_.ind_finish;
    SAC_PRIM_DATA.eye_px_start = eye_px_filt(SAC_PRIM_DATA.ind_start);
    SAC_PRIM_DATA.eye_py_start = eye_py_filt(SAC_PRIM_DATA.ind_start);
    SAC_PRIM_DATA.eye_px_finish = eye_px_filt(SAC_PRIM_DATA.ind_finish);
    SAC_PRIM_DATA.eye_py_finish = eye_py_filt(SAC_PRIM_DATA.ind_finish);   
    SAC_PRIM_DATA.eye_amp_x = SAC_PRIM_DATA.eye_px_finish - SAC_PRIM_DATA.eye_px_start;
    SAC_PRIM_DATA.eye_amp_y = SAC_PRIM_DATA.eye_py_finish - SAC_PRIM_DATA.eye_py_start;
    SAC_PRIM_DATA.eye_amp_m = sqrt(SAC_PRIM_DATA.eye_amp_x.^2 + SAC_PRIM_DATA.eye_amp_y.^2);
    SAC_PRIM_DATA.eye_vm_max = eye_vm_filt(SAC_PRIM_DATA.ind_vmax);
    reaction_2K = SAC_PRIM_DATA.ind_start - idx_cue_present_after;
    SAC_PRIM_DATA.reaction = reaction_2K./2; % in ms
    % Distance btwn eye start position and start target
    diff_start = sqrt((SAC_PRIM_DATA.eye_px_start-trial_data.start_x).^2 +(SAC_PRIM_DATA.eye_py_start-trial_data.start_y).^2);
    % Dist. btwn. eye finish pos. and cue tgt.
    diff_finish = sqrt((SAC_PRIM_DATA.eye_px_finish-trial_data.cue_x).^2 +(SAC_PRIM_DATA.eye_py_finish-trial_data.cue_y).^2);
    % If the saccade doesn't meet req., invalidiate
    if ((diff_start >= threshold_pos) || (diff_finish >= threshold_pos))
        SAC_PRIM_DATA.validity = 0;
    end
    % Save data
    trial_data.time_2K = time_2K;
    trial_data.tgt_x_data = tgt_x_data;
    trial_data.tgt_y_data = tgt_y_data;
    trial_data.eye_px_filt = eye_px_filt;
    trial_data.eye_py_filt = eye_py_filt;
    trial_data.eye_vx_filt = eye_vx_filt;
    trial_data.eye_vy_filt = eye_vy_filt;
    trial_data.eye_vm_filt = eye_vm_filt;
    trial_data.SAC_PRIM_DATA = SAC_PRIM_DATA;
    
    data.trial_data.(trial_field_name) = trial_data; 
end
%% (OLD) PLOT PRIM. SAC. SCATTER PLOT, REACTION TIME, SAC. MAX. SPEED
sac_px = {}; sac_py = {};
reaction = []; vm_max = [];
valid_sac_count = 0;
num_trial = length(fieldnames(data.trial_data));

for t = 1:num_trial
    trial_field_name = sprintf('trial_%d',t);
    trial_data = data.trial_data.(trial_field_name);
    eye_px_filt = trial_data.eye_px_filt;
    eye_py_filt = trial_data.eye_py_filt;
    % Plot only valid primSac.
    if trial_data.SAC_PRIM_DATA.validity == 1
        valid_sac_count = valid_sac_count + 1;
        % Cue
        cue_p(valid_sac_count,:) = [trial_data.cue_x, trial_data.cue_y];
        % Eye
        primSac_end_idx = trial_data.SAC_PRIM_DATA.ind_finish;
        eye_p_end(valid_sac_count,:) = [eye_px_filt(primSac_end_idx), eye_py_filt(primSac_end_idx)];
        primSac_idx = trial_data.SAC_PRIM_DATA.inds;
        sac_px{valid_sac_count} = trial_data.eye_px_filt(trial_data.SAC_PRIM_DATA.ind_start:trial_data.SAC_PRIM_DATA.ind_finish);
        sac_py{valid_sac_count} = trial_data.eye_py_filt(trial_data.SAC_PRIM_DATA.ind_start:trial_data.SAC_PRIM_DATA.ind_finish);
        reaction = [reaction, trial_data.SAC_PRIM_DATA.reaction];
        vm_max = [vm_max, trial_data.SAC_PRIM_DATA.eye_vm_max];
    end
end
% PrimSac. scatter plot
prim_sac_scatter_fig = figure;
prim_sac_scatter_axes = subplot(3,2,[1:4]); hold(prim_sac_scatter_axes,'on');

num_target = 8;
target_colors = [0 0.4470 0.7410;
                 0.8500 0.3250 0.0980;
                 0.9290 0.6940 0.1250;
                 0.4940 0.1840 0.5560;
                 0.4660 0.6740 0.1880;
                 0.3010 0.7450 0.9330;
                 0.6350 0.0780 0.1840;
                 0 0 0];
unique_cue_p = unique(cue_p,'rows');
for t = 1:valid_sac_count
    cue_p_trial = cue_p(t,:);
    eye_p_end_trial = eye_p_end(t,:);
    
    unique_cue_idx = find(unique_cue_p == cue_p_trial,1,'first');
    cue_plot = plot(prim_sac_scatter_axes,cue_p_trial(1),cue_p_trial(2),'*','Color',target_colors(unique_cue_idx,:),'MarkerSize',10);
    eye_end_plot = plot(prim_sac_scatter_axes,eye_p_end_trial(1),eye_p_end_trial(2),'o','Color',target_colors(unique_cue_idx,:),'MarkerSize',10);
    eye_plot = plot(prim_sac_scatter_axes,sac_px{t},sac_py{t},'k');
end
prim_sac_scatter_axes.XLim = [trial_data.start_x-6,trial_data.start_x+6]; prim_sac_scatter_axes.XLabel.String = 'Pos. (deg)';
prim_sac_scatter_axes.YLim = [trial_data.start_y-6,trial_data.start_y+6]; prim_sac_scatter_axes.YLabel.String = 'Pos. (deg)';
prim_sac_scatter_axes.Title.String = {'Prim. sac. trajectory & corresponding target',...
                                        sprintf('Num. trial: %d. Num valid sac.: %d.',num_trial, valid_sac_count)};
% Max. sac. speed
prim_sac_vm_max_axes = subplot(3,2,5);
histogram(prim_sac_vm_max_axes, vm_max);
prim_sac_vm_max_axes.XLabel.String = 'Speed (deg/s)'; prim_sac_vm_max_axes.YLabel.String = '#';
prim_sac_vm_max_axes.Title.String = {'Max. sac. speed',sprintf('Target dist.:  deg')};
% Reaction time
prim_sac_reaction_axes = subplot(3,2,6);
histogram(prim_sac_reaction_axes, reaction);
prim_sac_reaction_axes.XLabel.String = 'Time (ms)'; prim_sac_reaction_axes.YLabel.String = '#';
prim_sac_reaction_axes.Title.String = 'Reaction time (not actual)';
% %% PLOT TRIAL
% trial_fig = figure;
% trial_axes = axes(trial_fig); hold(trial_axes,'on');
% good_trial_list = [1,5]; 
% bad_trial_list = [2,3];
% for ax_idx = 1:2
%     trial_ax(ax_idx) = subplot(3,2,4+ax_idx);
%     hold(trial_ax(ax_idx),'on');
% end
% trial_ax(1) = subplot(3,2,5); trial_ax(2) = subplot(3,2,6);
% ax_idx = 0;
% for t = bad_trial_list
%     ax_idx = ax_idx + 1;
%     trial_field_name = sprintf('trial_%d',t);
%     trial_data = data.trial_data.(trial_field_name);
%     vpixx_idx_cue_present = find(trial_data.time_2K >= trial_data.time_cue_present_before(end),1,'first');
%     vpixx_idx_saccade_end = find(trial_data.time_2K >= trial_data.time_saccade_end(end),1,'first');
%     plot_idc = vpixx_idx_cue_present-1000:vpixx_idx_saccade_end+400;
%     tgt_x_plot = plot(trial_ax(ax_idx),trial_data.time_2K(plot_idc),trial_data.tgt_x_data(plot_idc),'--r','LineWidth',2);
%     tgt_y_plot = plot(trial_ax(ax_idx),trial_data.time_2K(plot_idc),trial_data.tgt_y_data(plot_idc),'--b','LineWidth',2);
%     eye_x_plot = plot(trial_ax(ax_idx),trial_data.time_2K(plot_idc),trial_data.eye_px_filt(plot_idc),'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
%     eye_y_plot = plot(trial_ax(ax_idx),trial_data.time_2K(plot_idc),trial_data.eye_py_filt(plot_idc),'Color',[0 0.4470 0.7410],'LineWidth',2);
%     xline(trial_ax(ax_idx),trial_data.time_2K(trial_data.SACS_PRIM_DATA.ind_start),'--','primSac. onset','LabelHorizontalAlignment','center');
%     xline(trial_ax(ax_idx),trial_data.time_2K(trial_data.SACS_PRIM_DATA.ind_finish),'--','primSac. offset','LabelHorizontalAlignment','center');
%     trial_ax(ax_idx).YLim = [min([trial_data.tgt_y_data(plot_idc),trial_data.tgt_x_data(plot_idc)])-4,max([trial_data.tgt_y_data(plot_idc),trial_data.tgt_x_data(plot_idc)])+10];
%     legend([tgt_x_plot,tgt_y_plot,eye_x_plot,eye_y_plot],{'tgt_x','tgt_y','eye_x','eye_y'},'Box','off');
%     xlabel(trial_ax(ax_idx),'Time (s)'); ylabel(trial_ax(ax_idx),'Pos. (deg)');
%     pause
%     cla(trial_axes);
% end

%% (OLD) REACTION TIME OF CORR SAC
num_trial = sum(arrayfun(@(z) contains(z,'trial'), fieldnames(data)));
for counter_trial = 1 : num_trial
    disp(counter_trial);
    trial_data = data.(sprintf('trial_%d',counter_trial));

    cal_matrix = trial_data.cal_matrix; % may change every trial due to real time correction
    if size(cal_matrix,1) == 1
        cal_matrix = reshape(cal_matrix,3,3);
        cal_matrix = cal_matrix';
    end
    vpixx_time = trial_data.vpixx_time_data;
    time_2K    = vpixx_time(1):0.0005:vpixx_time(end);
    
    invalid_idx = isnan(vpixx_time);

    eye_px_raw = trial_data.(sprintf('eye_%sx_raw_data',which_eye));
    eye_py_raw = trial_data.(sprintf('eye_%sy_raw_data',which_eye));
    eye_blink  = trial_data.(sprintf('eye_%s_blink_data',which_eye));

    % Remove invalid indices (mainly from blinking)
    invalid_idx = invalid_idx | isnan(eye_px_raw);
    invalid_idx = invalid_idx | isnan(eye_py_raw);
    invalid_idx = invalid_idx | eye_blink;
    eye_px_raw(invalid_idx) = []; 
    eye_py_raw(invalid_idx) = [];
    vpixx_time(invalid_idx) = [];

    % Resampling at 2000 Hz
    eye_px_raw = interp1(vpixx_time, eye_px_raw, time_2K,'linear','extrap');
    eye_py_raw = interp1(vpixx_time, eye_py_raw, time_2K, 'linear','extrap');
    eye_p_raw  = [eye_px_raw', eye_py_raw', ones(length(eye_px_raw),1)];

    % Transform to screen coordinates 
    eye_p  = eye_p_raw*cal_matrix;
    eye_px = eye_p(:,1);
    eye_py = eye_p(:,2);

    % Filter eye data
    eye_px_filt = filtfilt(b_butter,a_butter,eye_px);
    eye_py_filt = filtfilt(b_butter,a_butter,eye_py);

    % Get real time tgt data and resample @ 2000 Hz 
    time_real   = trial_data.tgt_time_data; % time registered during FSM loop
    tgt_px_real = trial_data.tgt_x_data;
    tgt_py_real = trial_data.tgt_y_data;
    tgt_px      = zeros(size(time_2K));
    tgt_py      = zeros(size(time_2K));
    for counter_t = 1 : length(time_real)
        t = time_real(counter_t);

        post_idx_2K = find(time_2K >= t, 1, 'first');

        tgt_px(pre_idx_2K:post_idx_2K) = tgt_px_real(counter_t);
        tgt_py(pre_idx_2K:post_idx_2K) = tgt_py_real(counter_t);

        pre_idx_2K = post_idx_2K;
    end
    
    % Compute eye speed
    eye_vx_filt = diff(eye_px_filt)./diff(time_2K'); 
    eye_vx_filt = [eye_vx_filt(1); eye_vx_filt]; % make it the same length as position data
    eye_vy_filt = diff(eye_py_filt)./diff(time_2K'); 
    eye_vy_filt = [eye_vy_filt(1); eye_vy_filt];
    eye_vm_filt = sqrt(eye_vx_filt.^2+eye_vy_filt.^2);
    
    % Save to trial struct.
    trial_data.time_2K        = time_2K;
    trial_data.eye_px_filt    = eye_px_filt;
    trial_data.eye_py_filt    = eye_py_filt;
    trial_data.eye_vx_filt    = eye_vx_filt;
    trial_data.eye_vy_filt    = eye_vy_filt;
    trial_data.eye_vm_filt    = eye_vm_filt;
    trial_data.start_x        = data.horz_offset; % to be fixed; save directly from experiment
    trial_data.start_y        = data.vert_offset;
    trial_data.tgt_px_filt    = tgt_px;
    trial_data.tgt_py_filt    = tgt_py;
    trial_data.eye_blink_filt = eye_blink;

    % Find saccades
    prim_sac = find_prim_sac(trial_data, threshold_pos, params_prim);
    corr_sac = find_corr_sac(trial_data, threshold_pos, params_corr);
    
    reaction_time(counter_trial) = corr_sac.ind_start - prim_sac.ind_finish;

end
%% (OLD) PLOT SACCADE MAIN SEQUENCE
flag_plot_trial=0;
[SACS_ALL_DATA] = JSP_Sac_Sorter_Vpixx(data,flag_plot_trial);

ramon_data = load('C:\Users\jays3\OneDrive - Johns Hopkins\Shadmehr Lab\Neurophysiology\Data\population_data\2020-12\2020-12-14\2020-12-14_12-49-15\analyzed_data\201214_124915_ANALYZED_ALL_SACS');

mirza_DATA = load('C:\Users\jays3\OneDrive - Johns Hopkins\Shadmehr Lab\Neurophysiology\Data\dwell_time\Mirza\2021-02-01\2021-02-01_13-11-53\210201_131153_random_target_backstep_dual_pump_131934_ANALYZED_corrSac_sorted');
[mirza_data.SACS_ALL_DATA, TRIALS_DATA, EXPERIMENT_PARAMS] = ESN_Sac_Sorter(mirza_DATA.TRIALS_DATA, mirza_DATA.EXPERIMENT_PARAMS);
figure;
tag = 1; % other_irrelv
sac_idx = ramon_data.SACS_ALL_DATA.tag == tag;
scatter(ramon_data.SACS_ALL_DATA.eye_r_amp_m(sac_idx), ramon_data.SACS_ALL_DATA.eye_r_vm_max(sac_idx),150,'.'); hold;
sac_idx = mirza_data.SACS_ALL_DATA.tag == tag;
scatter(mirza_data.SACS_ALL_DATA.eye_r_amp_m(sac_idx), mirza_data.SACS_ALL_DATA.eye_r_vm_max(sac_idx),150,'.');
scatter(SACS_ALL_DATA.eye_amp_m, SACS_ALL_DATA.eye_vm_max,50,[0.4660 0.6740 0.1880],'.'); 
legend('Ramon','Mirza', 'Ada');
figure;
scatter(ramon_data.SACS_ALL_DATA.eye_r_amp_m, ramon_data.SACS_ALL_DATA.eye_r_vm_max,150,'.'); hold;
scatter(mirza_data.SACS_ALL_DATA.eye_r_amp_m, mirza_data.SACS_ALL_DATA.eye_r_vm_max,150,'.');
scatter(SACS_ALL_DATA.eye_amp_m, SACS_ALL_DATA.eye_vm_max,50,[0.4660 0.6740 0.1880],'.'); 
legend('Ramon','Mirza', 'Ada');
figure;
sac_idx = ramon_data.SACS_ALL_DATA.tag == tag;
sac_idx = sac_idx & (ramon_data.SACS_ALL_DATA.eye_r_amp_m < 5);
histogram(ramon_data.SACS_ALL_DATA.eye_r_vm_max(sac_idx))
title('Ramon less than 5 deg; task irrelevant');
figure;
sac_idx = SACS_ALL_DATA.eye_amp_m < 5;
histogram(SACS_ALL_DATA.eye_vm_max(sac_idx))
title('Ada less than 5 deg');
%% FIND PRIMARY SACCADE
function prim_sac = find_prim_sac(trial_data, threshold_pos, params_prim)
    time_start = trial_data.state_start_t_cue_tgt_present;
    time_end   = trial_data.state_start_t_deliver_rew;
    idx_start  = find(trial_data.time_2K >= time_start(end),1,'first');
    idx_end    = find(trial_data.time_2K >= time_end(end),1,'first');
    
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