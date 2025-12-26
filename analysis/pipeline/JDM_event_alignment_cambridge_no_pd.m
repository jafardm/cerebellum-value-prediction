function JDM_event_alignment_cambridge_no_pd(path_to_raw, path_to_save)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cambridge probe, chamber A 
% Creates behavior and EPHYS events, then aligns them in time such that
% eventually, the data we have from monkey's behavior is well
% aligned in time to the EPHYS signals.
%
% Inputs:    path_to_raw: path to raw data. 
%                         E.g. ['Y:\', filesep, 'EPHYS', filesep, 'data_132F', filesep, ..., 'raw_data']
%
%            path_to_save: path to save EPHYS event data
%
%
% Output:   Two figures. Figure1 shows the state_combined alignment (that uses a 
%           combination of random signal and photodiode command to align behavior and EPHYS).
%           Figure2 shows the Photodiode_combined alignment (that aligns
%           Behavior's photodiode command to EPHYS's actual photodiode signal).
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load EPHYS EVENT DATA
file_name = 'all_channels.events';
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(path_to_raw(end), filesep);path_to_raw = [path_to_raw filesep];end
[ch_data, ch_time, ch_info] = load_open_ephys_data([path_to_raw file_name]);
EPHYS.CH_EVE.ch_data = ch_data;
EPHYS.CH_EVE.ch_time = ch_time;
EPHYS.CH_EVE.ch_info = ch_info;
EPHYS.file_name_CH_EVE = file_name;
EPHYS.file_path_CH_EVE = path_to_raw;
fprintf(' --> Completed. \n')
EPHYS.path_to_save = path_to_save;

EPHYS.debug_figures = false;

%% load a .continuous file
file_name_ = dir(fullfile(path_to_raw,'*_1.continuous'));
file_name = file_name_(1).name;
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(path_to_raw(end), filesep);path_to_raw = [path_to_raw filesep];end
[~, ch_time, ~] = load_open_ephys_data([path_to_raw file_name]);
EPHYS.time_30K   = double(ch_time);
fprintf(' --> Completed. \n')

%% load BEHAVE DATA
file_dir = dir([path_to_raw '*.mat']);
fprintf(['Loading ', file_dir.name, ' ... ']);
data = load(fullfile(file_dir.folder,file_dir.name),'data');
BEHAVE = data.data;
fprintf(' --> Completed. \n')

%% Build BEHAVE Alignment EVENTS
clearvars -except EPHYS BEHAVE name_time_device device_freq cutoff_freq
fprintf(['Building BEHAVE Alignment events', ' ... ']);

device_freq = 1e3;
name_time_device = sprintf('time_%dK',round(device_freq/1000));

time_reference = ESN_Round(double(BEHAVE.device_time_data)./1000, 0.001);
length_time    = length(time_reference);
time_regular   = time_reference(1) : 1e-3 : time_reference(end);
BEHAVE.Alignment.(name_time_device)=time_reference;

% Convert the states of EyeLink's digital input from decimal to binary
din_state_data = dec2bin(BEHAVE.din_state_data);



% Random signal
rnd_signal_ch = 7;
rnd_signal = str2num(din_state_data(:,end-rnd_signal_ch));
time_random_signal_rise     = ESN_Round( double(BEHAVE.din_time_data(rnd_signal == 1))./1000 , 0.001);
time_random_signal_fall     = ESN_Round( double(BEHAVE.din_time_data(rnd_signal == 0))./1000 , 0.001);

% PD command signal
dout_pd_signal_ch = 5;
dout_pd_signal = str2num(din_state_data(:,end-dout_pd_signal_ch));
time_dout_photodiode_rise     = ESN_Round( double(BEHAVE.din_time_data(dout_pd_signal == 1))./1000 , 0.001);
time_dout_photodiode_fall     = ESN_Round( double(BEHAVE.din_time_data(dout_pd_signal == 0))./1000 , 0.001);

% Actual PD signal (Replaced with PD command)
%pd_signal_ch = 0;
%pd_signal = str2num(din_state_data(:,end-pd_signal_ch));
time_photodiode_rise     = time_dout_photodiode_rise;
time_photodiode_fall     = time_dout_photodiode_fall;

variable_list = {...
    '_random_signal_rise' ,'_random_signal_fall', ...
    '_dout_photodiode_rise','_dout_photodiode_fall',...
    '_photodiode_rise','_photodiode_fall'};


for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'BEHAVE.Alignment.time' variable_name ' = ' 'time' variable_name ';']);
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'time_temp_' ' = ' 'time' variable_name ';']);
    time_temp_(end+1) = max([time_regular(end), time_temp_(end)])+1;
    event_temp_       = false(length_time, 1);
    counter_temp_     = find(time_temp_ >= time_regular(1), 1, 'first');
    eval([ 'time'    variable_name ' = ' 'time_temp_'    ';']);
    eval([ 'event'   variable_name ' = ' 'event_temp_'   ';']);
    eval([ 'counter' variable_name ' = ' 'counter_temp_' ';']);
end
for counter_time_point = 1 : length_time
    time_point_     = time_regular(counter_time_point);
    
    % random_signal_rise
    if time_point_ >= time_random_signal_rise(  counter_random_signal_rise)
        event_random_signal_rise(    counter_time_point) = true;
        counter_random_signal_rise   = counter_random_signal_rise   + 1;
    end
    % random_signal_fall
    if time_point_ >= time_random_signal_fall(  counter_random_signal_fall)
        event_random_signal_fall(    counter_time_point) = true;
        counter_random_signal_fall   = counter_random_signal_fall   + 1;
    end
    % dout_photodiode_rise
    if time_point_ >= time_dout_photodiode_rise(  counter_dout_photodiode_rise)
        event_dout_photodiode_rise(    counter_time_point) = true;
        counter_dout_photodiode_rise   = counter_dout_photodiode_rise   + 1;
    end
    % dout_photodiode_fall
    if time_point_ >= time_dout_photodiode_fall(  counter_dout_photodiode_fall)
        event_dout_photodiode_fall(    counter_time_point) = true;
        counter_dout_photodiode_fall   = counter_dout_photodiode_fall   + 1;
    end
    % photodiode_rise
    if time_point_ >= time_photodiode_rise(     counter_photodiode_rise)
        event_photodiode_rise(       counter_time_point) = true;
        counter_photodiode_rise      = counter_photodiode_rise      + 1;
    end
    % photodiode_fall
    if time_point_ >= time_photodiode_fall(     counter_photodiode_fall)
        event_photodiode_fall(       counter_time_point) = true;
        counter_photodiode_fall      = counter_photodiode_fall      + 1;
    end
end

event_random_signal       = false(length_time, 1);
event_dout_photodiode     = false(length_time, 1);
event_photodiode          = false(length_time, 1);
flag_random_signal        = false;
flag_dout_photodiode      = false;
flag_photodiode           = false;
for counter_time_point = 1 : length_time
    flag_random_signal   = flag_random_signal  ||   event_random_signal_rise( counter_time_point);
    flag_random_signal   = flag_random_signal  && (~event_random_signal_fall( counter_time_point));
    flag_dout_photodiode = flag_dout_photodiode          ||   event_dout_photodiode_rise(         counter_time_point);
    flag_dout_photodiode = flag_dout_photodiode          && (~event_dout_photodiode_fall(         counter_time_point));
    flag_photodiode      = flag_photodiode          ||   event_photodiode_rise(         counter_time_point);
    flag_photodiode      = flag_photodiode          && (~event_photodiode_fall(         counter_time_point));
    event_random_signal( counter_time_point)  = flag_random_signal;
    event_dout_photodiode(counter_time_point) = flag_dout_photodiode;
    event_photodiode(    counter_time_point)  = flag_photodiode;   
end

% Reconstruct data
BEHAVE.random_signal = interp1(time_reference,double(event_random_signal),time_regular,'nearest','extrap');
BEHAVE.dout_pd_signal = interp1(time_reference,double(event_dout_photodiode),time_regular,'nearest','extrap');
BEHAVE.pd_signal = interp1(time_reference,double(event_photodiode),time_regular,'nearest','extrap');
BEHAVE.(name_time_device) = time_regular;

BEHAVE.Alignment.event_state_combined      = double(BEHAVE.dout_pd_signal) .*1 + double(BEHAVE.random_signal) .*2;
BEHAVE.Alignment.event_photodiode_combined = double(BEHAVE.dout_pd_signal) .*1;
BEHAVE.Alignment.(name_time_device)        = BEHAVE.(name_time_device);
fprintf(' --> Completed. \n')

%% Build EPHYS Alignment events
clearvars -except EPHYS BEHAVE name_time_device device_freq cutoff_freq
fprintf(['Building EPHYS Alignment events', ' ... ']);
EPHYS.CH_EVE.data = [EPHYS.CH_EVE.ch_time(:) EPHYS.CH_EVE.ch_data(:) EPHYS.CH_EVE.ch_info.eventId(:)];

time_reference      = ESN_Round(EPHYS.time_30K(1),0.001) : 1/device_freq : ESN_Round(EPHYS.time_30K(end),0.001);
length_time         = length(time_reference);

time_random_signal_rise = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 1) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_random_signal_fall = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 1) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
time_dout_photodiode_rise = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 0) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_dout_photodiode_fall = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 0) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
time_photodiode_rise     = time_dout_photodiode_rise + 8; % add 8ms delay to PD command
time_photodiode_fall     = time_dout_photodiode_fall + 8;

variable_list = {...
    '_random_signal_rise' ,'_random_signal_fall', ...
    '_dout_photodiode_rise','_dout_photodiode_fall',...
    '_photodiode_rise','_photodiode_fall'};

EPHYS.Alignment.(name_time_device) = time_reference;
for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'EPHYS.Alignment.time' variable_name ' = ' 'time' variable_name ';']);
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'time_temp_' ' = ' 'time' variable_name ';']);
    time_temp_(end+1) = max([time_reference(end), time_temp_(end)])+1;
    event_temp_       = false(length_time, 1);
    counter_temp_     = find(time_temp_ >= time_reference(1), 1, 'first');
    eval([ 'time'    variable_name ' = ' 'time_temp_'    ';']);
    eval([ 'event'   variable_name ' = ' 'event_temp_'   ';']);
    eval([ 'counter' variable_name ' = ' 'counter_temp_' ';']);
end
for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point);
    
    % random_signal_rise
    if time_point_ >= time_random_signal_rise(  counter_random_signal_rise)
        event_random_signal_rise(    counter_time_point) = true;
        counter_random_signal_rise   = counter_random_signal_rise   + 1;
    end
    % random_signal_fall
    if time_point_ >= time_random_signal_fall(  counter_random_signal_fall)
        event_random_signal_fall(    counter_time_point) = true;
        counter_random_signal_fall   = counter_random_signal_fall   + 1;
    end
    % dout_photodiode_rise
    if time_point_ >= time_dout_photodiode_rise(  counter_dout_photodiode_rise)
        event_dout_photodiode_rise(    counter_time_point) = true;
        counter_dout_photodiode_rise   = counter_dout_photodiode_rise   + 1;
    end
    % dout_photodiode_fall
    if time_point_ >= time_dout_photodiode_fall(  counter_dout_photodiode_fall)
        event_dout_photodiode_fall(    counter_time_point) = true;
        counter_dout_photodiode_fall   = counter_dout_photodiode_fall   + 1;
    end
    % photodiode_rise
    if time_point_ >= time_photodiode_rise(     counter_photodiode_rise)
        event_photodiode_rise(       counter_time_point) = true;
        counter_photodiode_rise      = counter_photodiode_rise      + 1;
    end
    % photodiode_fall
    if time_point_ >= time_photodiode_fall(     counter_photodiode_fall)
        event_photodiode_fall(       counter_time_point) = true;
        counter_photodiode_fall      = counter_photodiode_fall      + 1;
    end
end

event_random_signal       = false(length_time, 1);
event_dout_photodiode     = false(length_time, 1);
event_photodiode          = false(length_time, 1);
flag_random_signal        = false;
flag_dout_photodiode      = false;
flag_photodiode           = false;
for counter_time_point = 1 : length_time
    flag_random_signal   = flag_random_signal  ||   event_random_signal_rise( counter_time_point);
    flag_random_signal   = flag_random_signal  && (~event_random_signal_fall( counter_time_point));
    flag_dout_photodiode = flag_dout_photodiode          ||   event_dout_photodiode_rise(         counter_time_point);
    flag_dout_photodiode = flag_dout_photodiode          && (~event_dout_photodiode_fall(         counter_time_point));
    flag_photodiode      = flag_photodiode          ||   event_photodiode_rise(         counter_time_point);
    flag_photodiode      = flag_photodiode          && (~event_photodiode_fall(         counter_time_point));
    event_random_signal( counter_time_point)  = flag_random_signal;
    event_dout_photodiode(counter_time_point) = flag_dout_photodiode;
    event_photodiode(    counter_time_point)  = flag_photodiode;   
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'EPHYS.Alignment.event' variable_name ' = ' 'event' variable_name ';']);
end

event_state_combined = ...
    double(event_dout_photodiode)  .* 1 + double(event_random_signal) .*2;
event_photodiode_combined = ...
    double(event_photodiode)  .* 1;

EPHYS.Alignment.event_state_combined      = event_state_combined;
EPHYS.Alignment.event_photodiode_combined = event_photodiode_combined;

fprintf(' --> Completed. \n');
%% ALIGN EPHYS and BEHAVE state_combined through xcorr and dtw
clearvars -except EPHYS BEHAVE name_time_device device_freq cutoff_freq
fprintf(['Aligning EPHYS and BEHAVE state_combined', ' ... ']);
EPHYS_time_xK              = EPHYS.Alignment.(name_time_device);
EPHYS_time_30K             = EPHYS.time_30K;
BEHAVE_time_xK             = BEHAVE.Alignment.(name_time_device);
EPHYS_state_combined       = EPHYS.Alignment.event_state_combined;
BEHAVE_state_combined      = BEHAVE.Alignment.event_state_combined;

% state_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(EPHYS_state_combined, BEHAVE_state_combined); % cross-correlate signals with each other
dur_diff = length(BEHAVE_state_combined) - length(EPHYS_state_combined);
% the difference between starting ephys and behavior can not be more than 2
% min
xcorr_value(abs(xcorr_lag-dur_diff)>=12e4 & abs(xcorr_lag+dur_diff)>=12e4) = nan;
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);

if  sample_diff > 0 % EPHYS signal later than BEHAVE
    EPHYS_EB_xcorr_state_combined_xK  = EPHYS_state_combined(  abs(sample_diff):end);
    EPHYS_EB_xcorr_time_xK            = EPHYS_time_xK(            abs(sample_diff):end);
    BEHAVE_EB_xcorr_state_combined_xK = BEHAVE_state_combined;
    BEHAVE_EB_xcorr_time_xK           = BEHAVE_time_xK;
elseif sample_diff < 0 % EPHYS signal earlier than BEHAVE
    EPHYS_EB_xcorr_state_combined_xK  = EPHYS_state_combined;
    EPHYS_EB_xcorr_time_xK            = EPHYS_time_xK;
    BEHAVE_EB_xcorr_state_combined_xK = BEHAVE_state_combined( abs(sample_diff):end);
    BEHAVE_EB_xcorr_time_xK           = BEHAVE_time_xK(           abs(sample_diff):end);
end
% state_combined: make the vectors the same size
if length(BEHAVE_EB_xcorr_state_combined_xK) ~= length(EPHYS_EB_xcorr_state_combined_xK)
    min_length = min([ length(BEHAVE_EB_xcorr_state_combined_xK),  length(EPHYS_EB_xcorr_state_combined_xK)]);
    EPHYS_EB_xcorr_state_combined_xK  = EPHYS_EB_xcorr_state_combined_xK(  1:min_length);
    EPHYS_EB_xcorr_time_xK            = EPHYS_EB_xcorr_time_xK(            1:min_length);
    BEHAVE_EB_xcorr_state_combined_xK = BEHAVE_EB_xcorr_state_combined_xK( 1:min_length);
    BEHAVE_EB_xcorr_time_xK           = BEHAVE_EB_xcorr_time_xK(           1:min_length);
end

% %% state_combined: fdind the Dynamic Time Warp (DTW) between 2 time series
% low pass filter the signals to generate a sinusoid around rises and falls
% Depending on the frequency of change in random signal, apply this filter or not
% 'cutoff_freq' of 25 was of old
cutoff_freq = nan;
if isnan(cutoff_freq)
    EPHYS_EB_xcorr_state_combined_xK_filt = EPHYS_EB_xcorr_state_combined_xK;
    BEHAVE_EB_xcorr_state_combined_xK_filt = BEHAVE_EB_xcorr_state_combined_xK;
else
    [b_butter,a_butter] = butter(3,(cutoff_freq/(device_freq/2)), 'low');
    EPHYS_EB_xcorr_state_combined_xK_filt  = filtfilt(b_butter,a_butter,EPHYS_EB_xcorr_state_combined_xK);
    BEHAVE_EB_xcorr_state_combined_xK_filt = filtfilt(b_butter,a_butter,BEHAVE_EB_xcorr_state_combined_xK);
end

% break the dtw analayises to smaller chunks, dtw does not work with large vectors
ind_edge_width = ceil(length(EPHYS_EB_xcorr_time_xK ) / 20e3);
ind_edges = round(linspace(1, length(EPHYS_EB_xcorr_time_xK ), ind_edge_width));
ind_edges(1) = 0;
% init and loop over chunks
EPHYS_EB_inds_DTW  = cell((length(ind_edges)-1), 1);
BEHAVE_EB_inds_DTW = cell((length(ind_edges)-1), 1);
parfor counter_chunk = 1 : (length(ind_edges)-1)
    inds_chunk = ( (ind_edges(counter_chunk)+1) : 1 : (ind_edges(counter_chunk+1)) )';
    EPHYS_EB_state_combined_chunk  = EPHYS_EB_xcorr_state_combined_xK_filt(inds_chunk);
    BEHAVE_EB_state_combined_chunk = BEHAVE_EB_xcorr_state_combined_xK_filt(inds_chunk);
%     s1 = subplot(2,1,1);
%     plot(EPHYS_EB_state_combined_chunk);
%     s2 = subplot(2,1,2);
%     plot(BEHAVE_EB_state_combined_chunk);
%     linkaxes([s1,s2]);
    [~,ix,iy] = dtw(EPHYS_EB_state_combined_chunk,BEHAVE_EB_state_combined_chunk, 50, 'absolute');  % allow upto 20ms warp
    EPHYS_EB_inds_DTW{counter_chunk}  = ix(:) + inds_chunk(1) - 1;
    BEHAVE_EB_inds_DTW{counter_chunk} = iy(:) + inds_chunk(1) - 1;
end
EPHYS_EB_inds_DTW  = cell2mat(EPHYS_EB_inds_DTW);
BEHAVE_EB_inds_DTW = cell2mat(BEHAVE_EB_inds_DTW);
EPHYS_EB_inds      = ( 1 : 1 : length(EPHYS_EB_xcorr_time_xK ) )';
BEHAVE_EB_inds     = ( 1 : 1 : length(BEHAVE_EB_xcorr_time_xK) )';
% dtw works by replicating the inds to match the two signals, here we
% reverse the replicated inds to generate two matched signals but with the
% size of original signals.
EB_ind_convert_from_EPHYS_to_BEHAVE = nan(size(EPHYS_EB_inds));
EB_ind_convert_from_BEHAVE_to_EPHYS = nan(size(BEHAVE_EB_inds));
for counter_ind = 1 : 1 : length(EPHYS_EB_inds_DTW)
    ind_EPHYS_EB_DTW  = EPHYS_EB_inds_DTW(counter_ind);
    ind_BEHAVE_EB_DTW = BEHAVE_EB_inds_DTW(counter_ind);
    % 'EB_ind_convert_from_EPHYS_to_BEHAVE' - indices are BEHAVE; values are corresponding order in EPHYS indices
    % e.g., 'EB_ind_convert_from_EPHYS_to_BEHAVE(2)' - the order of EPHYS index corresponding to 2nd index in BEHAVE
    EB_ind_convert_from_EPHYS_to_BEHAVE(ind_BEHAVE_EB_DTW) = ind_EPHYS_EB_DTW;
    EB_ind_convert_from_BEHAVE_to_EPHYS(ind_EPHYS_EB_DTW)  = ind_BEHAVE_EB_DTW;
end
% Find indices in 30 KHz time points (Ephys), corresponding to times in x KHz time points (Ephys)
time_reference      = EPHYS_time_30K(:);
length_time         = length(time_reference);
time_EPHYS_EB_xcorr_xK = EPHYS_EB_xcorr_time_xK;
time_EPHYS_EB_xcorr_xK(end+1) = max([time_reference(end), time_EPHYS_EB_xcorr_xK(end)])+1;
event_EPHYS_EB_xcorr_30K       = nan(length(EPHYS_EB_xcorr_time_xK), 1);
counter_EPHYS_EB_xcorr     = find(time_EPHYS_EB_xcorr_xK >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point); % 30 KHz
    if time_point_ >= time_EPHYS_EB_xcorr_xK(  counter_EPHYS_EB_xcorr)
        event_EPHYS_EB_xcorr_30K(    counter_EPHYS_EB_xcorr) = counter_time_point; % values here should theoretically be separated by 30 indices (30 KHz/ 1 KHz)
        counter_EPHYS_EB_xcorr   = counter_EPHYS_EB_xcorr   + 1;
    end
end
% Find indices in x KHz time points (Ephys), corresponding to times in x KHz time points (Ephys)
time_reference      = EPHYS_time_xK(:);
length_time         = length(time_reference);
time_EPHYS_EB_xcorr_xK = EPHYS_EB_xcorr_time_xK;
time_EPHYS_EB_xcorr_xK(end+1) = max([time_reference(end), time_EPHYS_EB_xcorr_xK(end)])+1;
event_EPHYS_EB_xcorr_xK       = nan(length(EPHYS_EB_xcorr_time_xK), 1);
counter_EPHYS_EB_xcorr     = find(time_EPHYS_EB_xcorr_xK >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_EPHYS_EB_xcorr_xK(  counter_EPHYS_EB_xcorr)
        event_EPHYS_EB_xcorr_xK(    counter_EPHYS_EB_xcorr) = counter_time_point;
        counter_EPHYS_EB_xcorr   = counter_EPHYS_EB_xcorr   + 1;
    end
end
% Find indices in x KHz time points (BEHAVE), corresponding to times in x KHz time points (BEHAVE)
time_reference      = BEHAVE_time_xK(:);
length_time         = length(time_reference);
time_BEHAVE_EB_xcorr_xK = BEHAVE_EB_xcorr_time_xK;
time_BEHAVE_EB_xcorr_xK(end+1) = max([time_reference(end), time_BEHAVE_EB_xcorr_xK(end)])+1;
event_BEHAVE_EB_xcorr_xK       = nan(length(BEHAVE_EB_xcorr_time_xK), 1);
counter_BEHAVE_EB_xcorr     = find(time_BEHAVE_EB_xcorr_xK >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_BEHAVE_EB_xcorr_xK(  counter_BEHAVE_EB_xcorr)
        event_BEHAVE_EB_xcorr_xK(    counter_BEHAVE_EB_xcorr) = counter_time_point;
        counter_BEHAVE_EB_xcorr   = counter_BEHAVE_EB_xcorr   + 1;
    end
end

EPHYS_EB_xcorr_ind_30K   = event_EPHYS_EB_xcorr_30K;
EPHYS_EB_xcorr_ind_xK    = event_EPHYS_EB_xcorr_xK;
BEHAVE_EB_xcorr_ind_xK   = event_BEHAVE_EB_xcorr_xK;
EPHYS_EB_aligned_ind_30K = EPHYS_EB_xcorr_ind_30K(EB_ind_convert_from_EPHYS_to_BEHAVE); % indices are BEHAVE indices; values are EPHYS indices 
EPHYS_EB_aligned_ind_xK  = EPHYS_EB_xcorr_ind_xK( EB_ind_convert_from_EPHYS_to_BEHAVE); 
BEHAVE_EB_aligned_ind_xK = BEHAVE_EB_xcorr_ind_xK(EB_ind_convert_from_BEHAVE_to_EPHYS); % indices are EPHYS indices; values are BEHAVE indices 

EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_30K            = EPHYS_EB_aligned_ind_30K;
EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K             = EPHYS_EB_aligned_ind_xK;
EPHYS.CH_EVE.align_states.BEHAVE_EB_aligned_ind_1K            = BEHAVE_EB_aligned_ind_xK;
EPHYS.CH_EVE.align_states.EB_ind_convert_from_BEHAVE_to_EPHYS = EB_ind_convert_from_BEHAVE_to_EPHYS;
EPHYS.CH_EVE.align_states.EB_ind_convert_from_EPHYS_to_BEHAVE = EB_ind_convert_from_EPHYS_to_BEHAVE;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_time_1K              = EPHYS_EB_xcorr_time_xK;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_time_1K             = BEHAVE_EB_xcorr_time_xK;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_ind_30K              = EPHYS_EB_xcorr_ind_30K;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_ind_1K               = EPHYS_EB_xcorr_ind_xK;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K              = BEHAVE_EB_xcorr_ind_xK;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_state_combined_1K    = EPHYS_EB_xcorr_state_combined_xK;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_state_combined_1K   = BEHAVE_EB_xcorr_state_combined_xK;
% EPHYS.CH_EVE.align_states.state_description            = BEHAVE.Alignment.state_description;
fprintf(' --> Completed. \n');

if EPHYS.debug_figures
clf(figure(1));ax_(1)=subplot(2,1,1);plot(EPHYS_state_combined);ax_(2)=subplot(2,1,2);plot(BEHAVE_state_combined);linkaxes(ax_,'x');
clf(figure(2));subplot(2,1,1);plot(EPHYS_EB_xcorr_state_combined_xK);subplot(2,1,2);plot(BEHAVE_EB_xcorr_state_combined_xK);
clf(figure(3));ax_(1)=subplot(2,1,1);plot(EPHYS_EB_xcorr_state_combined_xK);ax_(2)=subplot(2,1,2);plot(BEHAVE_EB_xcorr_state_combined_xK(EB_ind_convert_from_BEHAVE_to_EPHYS));linkaxes(ax_,'x');
end
% Mandatory checking
figure;
ax_(1)=subplot(2,1,1);
plot(EPHYS_EB_xcorr_state_combined_xK);
hold on;
plot(BEHAVE_EB_xcorr_state_combined_xK(EB_ind_convert_from_BEHAVE_to_EPHYS));
title('EPHYS: digital events registered, BEHAVE: digital events registered') 
ax_(2)=subplot(2,1,2);
plot(EPHYS_EB_xcorr_state_combined_xK' - ...
    BEHAVE_EB_xcorr_state_combined_xK...
    (EB_ind_convert_from_BEHAVE_to_EPHYS));
title({'difference','converted to most closely matched EPHYS time domain',...
       '2 patterns should match; otherwise something incorrect','few missing events are acceptable'})
linkaxes(ax_,'xy');
xlabel('sample');
ylim([min(EPHYS_state_combined)-1.5,max(EPHYS_state_combined)+0.5])
% xlim([5000,10000]);
sgtitle('Digital alignment: press any key to proceed');
pause
close 

disp(['length EPHYS: ' num2str( EPHYS.Alignment.(name_time_device)(end)-EPHYS.Alignment.(name_time_device)(1) )])
disp(['length BEHAVE: ' num2str( BEHAVE.Alignment.(name_time_device)(end)-BEHAVE.Alignment.(name_time_device)(1) )])
disp(['xcorr diff: ' num2str( sample_diff )])

%% ALIGN EPHYS and BEHAVE photodiode_combined through xcorr and dtw
clearvars -except EPHYS BEHAVE name_time_device device_freq cutoff_freq
fprintf(['Aligning EPHYS and BEHAVE photodiode_combined', ' ... ']);
EPHYS_time_xK              = EPHYS.Alignment.(name_time_device);
EPHYS_time_30K             = EPHYS.time_30K;
BEHAVE_time_xK             = BEHAVE.Alignment.(name_time_device);
EPHYS_photodiode_combined  = EPHYS.Alignment.event_photodiode_combined;
BEHAVE_photodiode_combined = BEHAVE.Alignment.event_photodiode_combined;

% photodiode_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(EPHYS_photodiode_combined-1/2, BEHAVE_photodiode_combined-1/2); % cross-correlate signals with each other
% dur_diff = length(BEHAVE_photodiode_combined) - length(EPHYS_photodiode_combined);
% % the difference between starting ephys and behavior can not be more than 2
% % min
% xcorr_value(abs(xcorr_lag-dur_diff)>=12e4 & abs(xcorr_lag+dur_diff)>=12e4) = nan;
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);
if sample_diff > 0
    EPHYS_PD_xcorr_photodiode_combined_xK  = EPHYS_photodiode_combined(  abs(sample_diff):end);
    EPHYS_PD_xcorr_time_xK                 = EPHYS_time_xK(                 abs(sample_diff):end);
    BEHAVE_PD_xcorr_photodiode_combined_xK = BEHAVE_photodiode_combined;
    BEHAVE_PD_xcorr_time_xK                = BEHAVE_time_xK;
elseif sample_diff < 0
    EPHYS_PD_xcorr_photodiode_combined_xK  = EPHYS_photodiode_combined;
    EPHYS_PD_xcorr_time_xK                 = EPHYS_time_xK;
    BEHAVE_PD_xcorr_photodiode_combined_xK = BEHAVE_photodiode_combined( abs(sample_diff):end);
    BEHAVE_PD_xcorr_time_xK                = BEHAVE_time_xK(                abs(sample_diff):end);
end
% photodiode_combined: make the vectors the same size
if length(BEHAVE_PD_xcorr_photodiode_combined_xK) ~= length(EPHYS_PD_xcorr_photodiode_combined_xK)
    min_length = min([ length(BEHAVE_PD_xcorr_photodiode_combined_xK),  length(EPHYS_PD_xcorr_photodiode_combined_xK)]);
    EPHYS_PD_xcorr_photodiode_combined_xK  = EPHYS_PD_xcorr_photodiode_combined_xK(  1:min_length);
    EPHYS_PD_xcorr_time_xK                 = EPHYS_PD_xcorr_time_xK(                 1:min_length);
    BEHAVE_PD_xcorr_photodiode_combined_xK = BEHAVE_PD_xcorr_photodiode_combined_xK( 1:min_length);
    BEHAVE_PD_xcorr_time_xK                = BEHAVE_PD_xcorr_time_xK(                1:min_length);
end

% photodiode_combined: find the Dynamic Time Warp (DTW) between 2 time series
% low pass filter the signals to generate a sinusoid around rises and falls
% Depending on the frequency of change in random signal, apply this filter or not
% 'cutoff_freq' of 25 was of old
cutoff_freq = nan;
if isnan(cutoff_freq)
    EPHYS_PD_xcorr_photodiode_combined_xK_filt  = EPHYS_PD_xcorr_photodiode_combined_xK;
    BEHAVE_PD_xcorr_photodiode_combined_xK_filt = BEHAVE_PD_xcorr_photodiode_combined_xK;
else
    [b_butter,a_butter] = butter(3,(cutoff_freq/(device_freq/2)), 'low');
    EPHYS_PD_xcorr_photodiode_combined_xK_filt  = filtfilt(b_butter,a_butter,EPHYS_PD_xcorr_photodiode_combined_xK);
    BEHAVE_PD_xcorr_photodiode_combined_xK_filt = filtfilt(b_butter,a_butter,BEHAVE_PD_xcorr_photodiode_combined_xK);
end

% break the dtw analayises to smaller chunks, dtw does not work with large vectors
ind_edge_width = ceil(length(EPHYS_PD_xcorr_time_xK ) / 5000);
ind_edges = round(linspace(1, length(EPHYS_PD_xcorr_time_xK ), ind_edge_width));
ind_edges(1) = 0;
% init and loop over chunks
EPHYS_PD_inds_DTW  = cell((length(ind_edges)-1), 1);
BEHAVE_PD_inds_DTW = cell((length(ind_edges)-1), 1);
parfor counter_chunk = 1 : (length(ind_edges)-1)
    inds_chunk = ( (ind_edges(counter_chunk)+1) : 1 : (ind_edges(counter_chunk+1)) )';
    EPHYS_PD_photodiode_combined_chunk  = EPHYS_PD_xcorr_photodiode_combined_xK_filt(inds_chunk);
    BEHAVE_PD_photodiode_combined_chunk = BEHAVE_PD_xcorr_photodiode_combined_xK_filt(inds_chunk);
    [~,ix,iy] = dtw(EPHYS_PD_photodiode_combined_chunk,BEHAVE_PD_photodiode_combined_chunk, 50, 'absolute'); % allow upto 50ms warp
    EPHYS_PD_inds_DTW{counter_chunk}  = ix(:) + inds_chunk(1) - 1;
    BEHAVE_PD_inds_DTW{counter_chunk} = iy(:) + inds_chunk(1) - 1;
end
EPHYS_PD_inds_DTW  = cell2mat(EPHYS_PD_inds_DTW);
BEHAVE_PD_inds_DTW = cell2mat(BEHAVE_PD_inds_DTW);
EPHYS_PD_inds      = ( 1 : 1 : length(EPHYS_PD_xcorr_time_xK ) )';
BEHAVE_PD_inds     = ( 1 : 1 : length(BEHAVE_PD_xcorr_time_xK) )';
% dtw works by replicating the inds to match the two signals, here we
% reverse the replicated inds to generate two matched signals but with the
% size of original signals.
PD_ind_convert_from_EPHYS_to_BEHAVE = nan(size(EPHYS_PD_inds));
PD_ind_convert_from_BEHAVE_to_EPHYS = nan(size(BEHAVE_PD_inds));
for counter_ind = 1 : 1 : length(EPHYS_PD_inds_DTW)
    ind_EPHYS_PD_DTW  = EPHYS_PD_inds_DTW(counter_ind);
    ind_BEHAVE_PD_DTW = BEHAVE_PD_inds_DTW(counter_ind);
    PD_ind_convert_from_EPHYS_to_BEHAVE(ind_BEHAVE_PD_DTW) = ind_EPHYS_PD_DTW;
    PD_ind_convert_from_BEHAVE_to_EPHYS(ind_EPHYS_PD_DTW)  = ind_BEHAVE_PD_DTW;
end

time_reference      = EPHYS_time_30K(:);
length_time         = length(time_reference);
time_EPHYS_PD_xcorr_xK = EPHYS_PD_xcorr_time_xK;
time_EPHYS_PD_xcorr_xK(end+1) = max([time_reference(end), time_EPHYS_PD_xcorr_xK(end)])+1;
event_EPHYS_PD_xcorr_30K       = nan(length(EPHYS_PD_xcorr_time_xK), 1);
counter_EPHYS_PD_xcorr     = find(time_EPHYS_PD_xcorr_xK >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_EPHYS_PD_xcorr_xK(  counter_EPHYS_PD_xcorr)
        event_EPHYS_PD_xcorr_30K(    counter_EPHYS_PD_xcorr) = counter_time_point;
        counter_EPHYS_PD_xcorr   = counter_EPHYS_PD_xcorr   + 1;
    end
end

time_reference      = EPHYS_time_xK(:);
length_time         = length(time_reference);
time_EPHYS_PD_xcorr_xK = EPHYS_PD_xcorr_time_xK;
time_EPHYS_PD_xcorr_xK(end+1) = max([time_reference(end), time_EPHYS_PD_xcorr_xK(end)])+1;
event_EPHYS_PD_xcorr_xK       = nan(length(EPHYS_PD_xcorr_time_xK), 1);
counter_EPHYS_PD_xcorr     = find(time_EPHYS_PD_xcorr_xK >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_EPHYS_PD_xcorr_xK(  counter_EPHYS_PD_xcorr)
        event_EPHYS_PD_xcorr_xK(    counter_EPHYS_PD_xcorr) = counter_time_point;
        counter_EPHYS_PD_xcorr   = counter_EPHYS_PD_xcorr   + 1;
    end
end

time_reference      = BEHAVE_time_xK(:);
length_time         = length(time_reference);
time_BEHAVE_PD_xcorr_xK = BEHAVE_PD_xcorr_time_xK;
time_BEHAVE_PD_xcorr_xK(end+1) = max([time_reference(end), time_BEHAVE_PD_xcorr_xK(end)])+1;
event_BEHAVE_PD_xcorr_xK       = nan(length(BEHAVE_PD_xcorr_time_xK), 1);
counter_BEHAVE_PD_xcorr     = find(time_BEHAVE_PD_xcorr_xK >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_BEHAVE_PD_xcorr_xK(  counter_BEHAVE_PD_xcorr)
        event_BEHAVE_PD_xcorr_xK(    counter_BEHAVE_PD_xcorr) = counter_time_point;
        counter_BEHAVE_PD_xcorr   = counter_BEHAVE_PD_xcorr   + 1;
    end
end

EPHYS_PD_xcorr_ind_30K   = event_EPHYS_PD_xcorr_30K;
EPHYS_PD_xcorr_ind_xK    = event_EPHYS_PD_xcorr_xK;
BEHAVE_PD_xcorr_ind_xK   = event_BEHAVE_PD_xcorr_xK;
EPHYS_PD_aligned_ind_30K = EPHYS_PD_xcorr_ind_30K(PD_ind_convert_from_EPHYS_to_BEHAVE);
EPHYS_PD_aligned_ind_xK  = EPHYS_PD_xcorr_ind_xK( PD_ind_convert_from_EPHYS_to_BEHAVE);
BEHAVE_PD_aligned_ind_xK = BEHAVE_PD_xcorr_ind_xK(PD_ind_convert_from_BEHAVE_to_EPHYS);

EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_30K               = EPHYS_PD_aligned_ind_30K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_1K                = EPHYS_PD_aligned_ind_xK;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_aligned_ind_1K               = BEHAVE_PD_aligned_ind_xK;
EPHYS.CH_EVE.align_photodiode.PD_ind_convert_from_BEHAVE_to_EPHYS    = PD_ind_convert_from_BEHAVE_to_EPHYS;
EPHYS.CH_EVE.align_photodiode.PD_ind_convert_from_EPHYS_to_BEHAVE    = PD_ind_convert_from_EPHYS_to_BEHAVE;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_time_1K                 = EPHYS_PD_xcorr_time_xK;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_time_1K                = BEHAVE_PD_xcorr_time_xK;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_ind_30K                 = EPHYS_PD_xcorr_ind_30K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_ind_1K                  = EPHYS_PD_xcorr_ind_xK;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_ind_1K                 = BEHAVE_PD_xcorr_ind_xK;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_photodiode_combined_1K  = EPHYS_PD_xcorr_photodiode_combined_xK;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_photodiode_combined_1K = BEHAVE_PD_xcorr_photodiode_combined_xK;
fprintf(' --> Completed. \n');

if EPHYS.debug_figures
clf(figure(1));ax_(1)=subplot(2,1,1);plot(EPHYS_photodiode_combined);ax_(2)=subplot(2,1,2);plot(BEHAVE_photodiode_combined);linkaxes(ax_,'x');
clf(figure(2));subplot(2,1,1);plot(EPHYS_PD_xcorr_photodiode_combined_xK);subplot(2,1,2);plot(BEHAVE_PD_xcorr_photodiode_combined_xK);
clf(figure(3));ax_(1)=suhbplot(2,1,1);plot(EPHYS_PD_xcorr_photodiode_combined_xK);ax_(2)=subplot(2,1,2);plot(BEHAVE_PD_xcorr_photodiode_combined_xK(PD_ind_convert_from_BEHAVE_to_EPHYS));linkaxes(ax_,'x');
end
% Mandatory checking
figure;
sgtitle('Photodiode alignment: press any key to proceed');
ax_(1)=subplot(2,1,1);
plot(EPHYS_PD_xcorr_photodiode_combined_xK);
hold on;
plot(BEHAVE_PD_xcorr_photodiode_combined_xK(PD_ind_convert_from_BEHAVE_to_EPHYS));
title('EPHYS: digital events registered, BEHAVE: digital events registered') 
ax_(2)=subplot(2,1,2);
plot(EPHYS_PD_xcorr_photodiode_combined_xK' - BEHAVE_PD_xcorr_photodiode_combined_xK(PD_ind_convert_from_BEHAVE_to_EPHYS));
title({'Difference','converted to most closely matched EPHYS time domain',...
       '2 patterns should match; otherwise something incorrect','few missing events are acceptable'})
linkaxes(ax_,'xy');
xlabel('sample');
ylim([min(EPHYS_photodiode_combined)-1.5,max(EPHYS_photodiode_combined)+0.5])
% xlim([5000,100000]);
pause
close 
%% Save EPHYS EVENT DATA
clearvars -except EPHYS BEHAVE name_time_device device_freq
EPHYS_time_1K    = EPHYS.Alignment.(name_time_device);
EPHYS_time_30K   = EPHYS.time_30K;
BEHAVE_time_1K   = BEHAVE.Alignment.(name_time_device);
align_photodiode = EPHYS.CH_EVE.align_photodiode;
align_states     = EPHYS.CH_EVE.align_states;
sampling_rate    = device_freq;

path_to_save = EPHYS.path_to_save;
file_name = EPHYS.file_name_CH_EVE;
[~, file_name, ~] = fileparts(file_name);
file_name = [file_name '_EVE1_aligned.mat'];
clearvars EPHYS BEHAVE
fprintf([file_name ': Saving EPHYS Event Data ...'])
save([path_to_save file_name], 'BEHAVE_time_1K','EPHYS_time_30K',...
    'EPHYS_time_1K','align_photodiode','align_states', 'file_name',...
    'path_to_save','sampling_rate','-v7.3');
fprintf(' --> Completed. \n')

end
