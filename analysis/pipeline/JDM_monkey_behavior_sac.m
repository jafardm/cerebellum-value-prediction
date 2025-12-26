function [sac_data, trials_data, meta_data] = JDM_monkey_behavior_sac...
    (mat_file_address, flag_old_app, recal_matrix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs saccade sorting and generates saccade map for all tags
%
% Inputs:    mat_file_address: behavior data file address
%            flag_old_app: to specify if we are using old or new trial organizer app
%            flag_figure: figure generation flag

% Outputs:   meta_data: includes parameters of the experiment including file address,
%                               # trials, tag list, sampling frequency, etc.
%            Sac sorter figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAF_params_funcs;

if isempty(recal_matrix)
    recal_matrix = eye(3);
    meta_data.flag_recal = 0;
else
    meta_data.flag_recal = 1;
end

% Get file_name and file_path
[path_to_raw_data,file_name,~] = fileparts(mat_file_address);
if ~strcmp(path_to_raw_data(end), filesep)
    path_to_raw_data = [path_to_raw_data filesep];
end

meta_data.path_to_rec   = path_to_raw_data(1:end-9);
meta_data.raw_data_file = file_name;

path_to_eye = [meta_data.path_to_rec, 'analyzed_data', filesep, ...
            'behavior_data', filesep, 'eye', filesep];

% Load Data
fprintf('Loading ...\n')
load(mat_file_address, 'data');
fprintf('%s: Loading complete.\n', file_name);
data_fieldnames = fieldnames(data);
meta_data.num_trials  = sum(arrayfun(@(str) contains(str,'trial'),...
    data_fieldnames));

% load Sac Finder Results
uneye_file = [path_to_eye file_name '_UNEYE.mat'];
converted_file = [path_to_eye file_name '_CONVERTED.mat'];
if not(exist(uneye_file, 'file') && exist(converted_file, 'file'))
    error('UNEYE data and converted data undetected, using original algorithm instead')
end

% load EPHYS EVENT DATA
file_name_ = dir([path_to_eye '*_EVE1_aligned.mat']);
if isempty(file_name_)
    event_data = [];
else
    file_name = file_name_(1).name;
    event_data = load([path_to_eye file_name]);
end

UNEYE_DATA = load(uneye_file,'Sac');
CONVERTED_DATA = load(converted_file,'data');
CONVERTED_DATA = CONVERTED_DATA.data;
% CONVERTED_DATA.data.eye_tracked = 'left';

meta_data.chamber = CONVERTED_DATA.chamber;
meta_data.eye_tracked = CONVERTED_DATA.eye_tracked;
meta_data.flag_old_app = flag_old_app;

% Get right_eye_tracked, and left_eye_tracked
if flag_old_app
    time_device = data.eyelink_time;
else
    % Record sampling frequency
    if strcmp(meta_data.chamber,'B')
        time_device  = data.trial_1.device_time_data;
    else
        time_device  = data.device_time_data;
    end
end

meta_data.sampling_freq = round(1/median(diff(time_device))); % in Hz

% Check if we are using the old behavior software or the new one
if flag_old_app == 1
    trials_data = MAF_trial_organizer_old_app_v2(data);
else
    trials_data = JDM_trial_organizer(data, meta_data);
end

% recalibrate and choose the tracked eye
trials_data_new = trials_data;
if strcmp(meta_data.eye_tracked,'left')
    eye_side = 'l';
else
    eye_side = 'r';
end
variable_list = { ...
    ['eye_',eye_side,'_px'], ['eye_',eye_side,'_py']; ...
    ['eye_',eye_side,'_px_filt'], ['eye_',eye_side,'_py_filt']; ...
    ['eye_',eye_side,'_vx'],['eye_',eye_side,'_vy'];...
    ['eye_',eye_side,'_vx_filt'],['eye_',eye_side,'_vy_filt']};
variable_list_new = {'eye_px', 'eye_py'; ...
    'eye_px_filt', 'eye_py_filt'; ...
    'eye_vx','eye_vy';...
    'eye_vx_filt','eye_vy_filt'};
trials_data_new = rmfield(trials_data_new,variable_list);
trials_data_new = rmfield(trials_data_new,{['eye_',eye_side,'_vm'],...
    ['eye_',eye_side,'_vm_filt']});

num_trial = length(trials_data.start_x);
for t_idx = 1:num_trial
    for v_idx = 1:size(variable_list,1)
        x_to_cal = trials_data.(variable_list{v_idx,1}){t_idx};
        y_to_cal = trials_data.(variable_list{v_idx,2}){t_idx};
        data_to_cal = [x_to_cal, y_to_cal, ones(size(x_to_cal))];
        data_cal = data_to_cal * recal_matrix;
        x_cal = data_cal(:,1);
        y_cal = data_cal(:,2);
        trials_data_new.(variable_list_new{v_idx,1}){t_idx} = x_cal;
        trials_data_new.(variable_list_new{v_idx,2}){t_idx} = y_cal;
        if contains(variable_list{v_idx,1},['eye_',eye_side,'_v'])
            if contains(variable_list{v_idx,1} ,'filt')
                variable_name = [variable_list_new{v_idx,1}(1:5),'m_filt'];
            else
                variable_name = [variable_list_new{v_idx,1}(1:5),'m'];
            end
            trials_data_new.(variable_name){t_idx} = sqrt(x_cal.^2 + y_cal.^2);
        end
    end
end

trials_data = trials_data_new;

% Sac Sorting
[sac_data, meta_data] = JDM_Sac_Sorter(trials_data, meta_data,...
    UNEYE_DATA, CONVERTED_DATA, event_data);