function [sac_data, meta_data] = JDM_Sac_Sorter_v2(trials_data, meta_data,UNEYE_DATA,CONVERTED_DATA,event_data)

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
counter_valid_trial = 1; 
for counter_trial = 1 : num_trials


% INITIALIZATION
clearvars -except meta_data trials_data counter_trial num_trials SACS_ALL ...
                    counter_valid_trial trial_num sampling_freq ...
                    UNEYE_DATA CONVERTED_DATA event_data

% EXTRACT TRIAL DATA

TRIAL = extract_trial_data(trials_data, counter_trial);

% EXTRACT SACCADES
clearvars -except meta_data trials_data counter_trial num_trials TRIAL SACS_ALL counter_valid_trial trial_num time_device...
                      sampling_freq r_eye_tracked l_eye_tracked UNEYE_DATA CONVERTED_DATA event_data
% Load saccade validity info
all_sac_validity_temp = UNEYE_DATA.Sac.(sprintf('trial_%d',counter_trial)).all_sac_validity;

% Skip this trial if there are no valid saccades
if isempty(all_sac_validity_temp)
    fprintf("No saccade in Trial %d\n", counter_trial);
    continue
end

% Get time device reference
time_device = trials_data.time_1K{1, counter_trial};

% Extract saccade timing indices
[all_sac_ind_onset, all_sac_ind_vmax, all_sac_ind_offset, all_sac_validity] = ...
    extract_saccade_indices(UNEYE_DATA, CONVERTED_DATA, counter_trial, time_device);

% CREATE SACCADES STRUCTURE
SACS_ALL_TRIAL = create_saccade_trial_struct( ...
    all_sac_validity, counter_trial, all_sac_ind_onset, ...
    all_sac_ind_vmax, all_sac_ind_offset, time_device, TRIAL, trials_data);

%% SET TAGGING THRESHOLDS
threshold_pos = 2.2;              % deg, position threshold (default 1.5, backstep 1.7)
threshold_ang = 45.0;             % deg, direction threshold
threshold_micro_saccade = 1.0;    % deg, micro saccade limit
flag_plot_trial = 0;

% TAG PRIMARY SACCADES (tag 1, 2, 3)
% 1: Primary Success — reaches target (low or high)
% 2: Primary Attempt — correct direction but fails to reach
% 3: Primary Failure — wrong direction or no attempt
SACS_ALL_TRIAL = tag_primary_saccades(SACS_ALL_TRIAL, TRIAL, ...
    threshold_pos, threshold_ang, threshold_micro_saccade);

% TAG CORRECTIVE SACCADES (tag 4, 5)
% Only applies to "jump" and "forced" trial conditions
% Scans between prim_success offset and next str_present (or trial end)
SACS_ALL_TRIAL = tag_corrective_saccades(SACS_ALL_TRIAL, TRIAL, ...
    threshold_pos, threshold_ang);

% TAG BACK-TO-CENTER SACCADES (tag 6, 7, 8)
% Scans before cue_present and after previous str_present (if any)
% 6: First return-to-center
% 7: Return to center after prim/corr
% 8: Unresolved center-related saccades (to be revisited)
SACS_ALL_TRIAL = back_to_center_saccades(TRIAL, SACS_ALL_TRIAL, threshold_pos);

% TAG IRRELEVANT SACCADES (tag 9, 10)
% 9: Target-Irrelevant — ends near target, but untagged
% 10: Other Irrelevant — untagged, ends away from target
SACS_ALL_TRIAL = tag_irrelevant_saccades(SACS_ALL_TRIAL, threshold_pos);


% TAG PRIMARY SACCADES WITH NO FOLLOW-UP CORRECTION (special case of tag 1)
% Marks prim_success not followed by corrective sac,
% but visual info at offset implies error (distance > threshold)
SACS_ALL_TRIAL = tag_primary_no_correction_saccades(SACS_ALL_TRIAL, ...
    TRIAL, threshold_pos);


% Fills missing visual saccade metrics and sets reaction times
% for 'target_irrelev' (tag==9) and 'other_irrelev' (tag==10) saccades.      
SACS_ALL_TRIAL = fill_nan_visual_and_set_irrelev_reaction(SACS_ALL_TRIAL, TRIAL);

% TAG 12
% Re-tag the irrelevant sac. possibly as db_corr_success; most are target_irrelevant sac. but some were found to be classified as 
% irrelevant sac. bc. it comes late in the trial, so the tgt. has jumped back to start by the time db_corr sac. lands at cue tgt. 
% @ the previous corr_success offset, tgt. should have jumped to a different location from where corr_success was aiming;
% dwell time could be long enough that tgt. doesn't jump back to cue bf. the first corr_success ends, in which case,
% there is no "corrective" saccade to be made
SACS_ALL_TRIAL = tag_post_correction_saccades(SACS_ALL_TRIAL, threshold_pos);

% TAG 13
   
% Tag the corr. sac. that is not followed by another corr. sac.
% Sac. is corr_success, but visual info. is what occurs at its offset
% Onset time is either the end of the trial or onset of next saccade
% Offset time is the same as that of corr_success: offset of sac.
% Visual time is the offset time
% @ its offset, tgt. should have jumped to a different location from where corr_success was aiming;
% dwell time could be long enough that tgt. doesn't jump back to cue bf. the first corr_success ends, in which case,
% there is no "corrective" saccade to be made

% Finds all corrective-success saccades (tag 4).
% For each, checks if:
% No valid next corrective saccade (i.e., not tag 12), and
% The corrective saccade failed to land near the target (>1° error).
% If both are true, creates a new saccade entry with tag 13.
% Appends all relevant fields properly, including updated visual angle and amplitude.
SACS_ALL_TRIAL = tag_unreached_corrective_saccades(SACS_ALL_TRIAL, TRIAL);

% TAG 6
    % Re-tag the 1st saccade as back_center_success if it landed at start
    % Handles case when there is 'back_center_success' but eye makes irrelevant saccades before the saccade twrd. target that initiates cue present.
    idx_sac = 1;
    SACS_ALL_TRIAL = tag_back_to_center_success(SACS_ALL_TRIAL, TRIAL, idx_sac, threshold_pos);

%  Fill up the count values
    SACS_ALL_TRIAL = tag_saccade_counting(SACS_ALL_TRIAL);

% Tag 14
    % Make another tag for 'other_irrelev' where the visual ang. is 
    % based on direction of center target from the offset pos. of previous saccade
 SACS_ALL_TRIAL = tag_virtual_visual_event_after_irrelevant_saccade(SACS_ALL_TRIAL);

  % Tag 15
  % Make another tag for 'back_center_irrelev' where the visual ang. is 
  % based on direction of center target from the offset pos. of previous saccade
  % Visual time is the same as that of tag 8
  SACS_ALL_TRIAL = tag_back_to_center_virtual_event(SACS_ALL_TRIAL, TRIAL);
    
% Debugging

plot_trial_debug(flag_plot_trial, counter_trial, TRIAL, SACS_ALL_TRIAL, time_device, meta_data);

 % Append trials
    SACS_ALL(counter_valid_trial) = SACS_ALL_TRIAL;
    counter_valid_trial = counter_valid_trial + 1;
    
end

% Arrange 'SACS_ALL_DATA'
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

% ADD FIXATION VARIABLES
sac_data = MAF_Fix_Sorter(sac_data,trials_data);

% FIX time_visual with photodiode
sac_data = fix_time_visual_with_photodiode(sac_data, event_data);

end

%% ================================  SUBFUNCTIONS ===========================

%% extract_trial_data
function TRIAL = extract_trial_data(trials_data, counter_trial)
% EXTRACT_TRIAL_DATA Extracts trial information from a trials_data struct for a specific trial index.

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
    TRIAL.stim_flag         = trials_data.stim_flag(  1,counter_trial);
    TRIAL.state_stim_t_start         = trials_data.state_stim_t_start(  1,counter_trial);

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
end

%% GET saccade inices from UNEYE
function [all_sac_ind_onset, all_sac_ind_vmax, all_sac_ind_offset, all_sac_validity] = ...
    extract_saccade_indices(UNEYE_DATA, CONVERTED_DATA, counter_trial, time_device)

    all_sac_ind_onset_temp  = UNEYE_DATA.Sac.(sprintf('trial_%d',counter_trial)).all_sac_ind_onset;
    all_sac_ind_vmax_temp   = UNEYE_DATA.Sac.(sprintf('trial_%d',counter_trial)).all_sac_ind_vmax;
    all_sac_ind_offset_temp = UNEYE_DATA.Sac.(sprintf('trial_%d',counter_trial)).all_sac_ind_offset;
    all_sac_validity_temp   = UNEYE_DATA.Sac.(sprintf('trial_%d',counter_trial)).all_sac_validity;

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

end
%% Cratea SACS_ALL data
function SACS_ALL_TRIAL = create_saccade_trial_struct(...
    all_sac_validity, counter_trial, all_sac_ind_onset, all_sac_ind_vmax, all_sac_ind_offset, ...
    time_device, TRIAL, trials_data)

    % Create the base structure
    SACS_ALL_TRIAL = struct;
    
    % Basic saccade information
    SACS_ALL_TRIAL.validity    = reshape(all_sac_validity, 1, []);
    SACS_ALL_TRIAL.trial_num   = reshape(ones(size(all_sac_validity)) * counter_trial, 1, []);
    SACS_ALL_TRIAL.tag         = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.count       = nan(size(all_sac_validity));
    SACS_ALL_TRIAL.flag_last_cue = false(size(all_sac_validity));
    
    % Timing information
    SACS_ALL_TRIAL.time_onset  = reshape(time_device(all_sac_ind_onset), 1, []);
    SACS_ALL_TRIAL.time_vmax   = reshape(time_device(all_sac_ind_vmax), 1, []);
    SACS_ALL_TRIAL.time_offset = reshape(time_device(all_sac_ind_offset), 1, []);
    
    % Position information
    SACS_ALL_TRIAL.start_x = ones(size(all_sac_ind_onset)) * TRIAL.start_x;
    SACS_ALL_TRIAL.start_y = ones(size(all_sac_ind_onset)) * TRIAL.start_y;
    SACS_ALL_TRIAL.cue_x = ones(size(all_sac_ind_onset)) * TRIAL.cue_x;
    SACS_ALL_TRIAL.cue_y = ones(size(all_sac_ind_onset)) * TRIAL.cue_y;
    SACS_ALL_TRIAL.end_x = ones(size(all_sac_ind_onset)) * TRIAL.end_x;
    SACS_ALL_TRIAL.end_y = ones(size(all_sac_ind_onset)) * TRIAL.end_y;
    
    % Stimulus information
    SACS_ALL_TRIAL.stim_flag = zeros(size(SACS_ALL_TRIAL.time_onset));
    SACS_ALL_TRIAL.state_stim_t_start = nan(size(all_sac_ind_onset));
    
    % Process stimulation times if they exist
    if ~iscell(trials_data.state_stim_t_start)
        trials_data.state_stim_t_start = num2cell(trials_data.state_stim_t_start, 1);
    end
    
    stim_times = trials_data.state_stim_t_start{1, counter_trial};
    
    if ~isnan(stim_times)
        for s = 1:length(stim_times)
            stim_time = stim_times(s);
            for k = 1:length(SACS_ALL_TRIAL.time_onset)
                if SACS_ALL_TRIAL.validity(k) == 1 && ...
                   stim_time >= SACS_ALL_TRIAL.time_onset(k) && ...
                   stim_time <= SACS_ALL_TRIAL.time_offset(k)
                   
                    SACS_ALL_TRIAL.stim_flag(k) = 1;
                    
                    if isnan(SACS_ALL_TRIAL.state_stim_t_start(k))
                        SACS_ALL_TRIAL.state_stim_t_start(k) = stim_time;
                    end
                end
            end
        end
    end
    
    % Reward task parameters
    SACS_ALL_TRIAL.tgt_num          = ones(size(all_sac_ind_onset)) * TRIAL.tgt_num;          
    SACS_ALL_TRIAL.cue_x_high_rew   = ones(size(all_sac_ind_onset)) * TRIAL.cue_x_high_rew;     
    SACS_ALL_TRIAL.cue_x_low_rew    = ones(size(all_sac_ind_onset)) * TRIAL.cue_x_low_rew;       
    SACS_ALL_TRIAL.cue_y_high_rew   = ones(size(all_sac_ind_onset)) * TRIAL.cue_y_high_rew;     
    SACS_ALL_TRIAL.cue_y_low_rew    = ones(size(all_sac_ind_onset)) * TRIAL.cue_y_low_rew;       
    SACS_ALL_TRIAL.high_rew_tgt_num = ones(size(all_sac_ind_onset)) * TRIAL.high_rew_tgt_num;  
    SACS_ALL_TRIAL.low_rew_tgt_num  = ones(size(all_sac_ind_onset)) * TRIAL.low_rew_tgt_num;   
    SACS_ALL_TRIAL.task_cond        = ones(size(all_sac_ind_onset)) * TRIAL.task_cond;
    SACS_ALL_TRIAL.choice           = ones(size(all_sac_ind_onset)) * TRIAL.choice;
    SACS_ALL_TRIAL.rew_cond         = ones(size(all_sac_ind_onset)) * TRIAL.rew_cond;
    SACS_ALL_TRIAL.jump_cond        = ones(size(all_sac_ind_onset)) * TRIAL.jump_cond;
    SACS_ALL_TRIAL.tgt_cond         = ones(size(all_sac_ind_onset)) * TRIAL.tgt_cond;
    
    % Additional metrics (initialized as NaN)
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
    
    % Duration calculation
    SACS_ALL_TRIAL.duration = reshape((SACS_ALL_TRIAL.time_offset - SACS_ALL_TRIAL.time_onset) * 1000.0, 1, []);
    
    % Eye movement metrics
    SACS_ALL_TRIAL.eye_vm_max     = reshape(TRIAL.eye_vm_filt(all_sac_ind_vmax), 1, []);
    SACS_ALL_TRIAL.eye_px_onset   = reshape(TRIAL.eye_px_filt(all_sac_ind_onset), 1, []);
    SACS_ALL_TRIAL.eye_px_offset  = reshape(TRIAL.eye_px_filt(all_sac_ind_offset), 1, []);
    SACS_ALL_TRIAL.eye_py_onset   = reshape(TRIAL.eye_py_filt(all_sac_ind_onset), 1, []);
    SACS_ALL_TRIAL.eye_py_offset  = reshape(TRIAL.eye_py_filt(all_sac_ind_offset), 1, []);
    SACS_ALL_TRIAL.eye_amp_x      = reshape(SACS_ALL_TRIAL.eye_px_offset - SACS_ALL_TRIAL.eye_px_onset, 1, []);
    SACS_ALL_TRIAL.eye_amp_y      = reshape(SACS_ALL_TRIAL.eye_py_offset - SACS_ALL_TRIAL.eye_py_onset, 1, []);
    SACS_ALL_TRIAL.eye_amp_m      = reshape(sqrt(SACS_ALL_TRIAL.eye_amp_x.^2 + SACS_ALL_TRIAL.eye_amp_y.^2), 1, []);
    SACS_ALL_TRIAL.eye_ang        = reshape(atan2d(SACS_ALL_TRIAL.eye_amp_y, SACS_ALL_TRIAL.eye_amp_x), 1, []);
    
    % Target positions
    SACS_ALL_TRIAL.tgt_px_onset   = reshape(TRIAL.tgt_px(all_sac_ind_onset), 1, []);
    SACS_ALL_TRIAL.tgt_px_offset  = reshape(TRIAL.tgt_px(all_sac_ind_offset), 1, []);
    SACS_ALL_TRIAL.tgt_py_onset   = reshape(TRIAL.tgt_py(all_sac_ind_onset), 1, []);
    SACS_ALL_TRIAL.tgt_py_offset  = reshape(TRIAL.tgt_py(all_sac_ind_offset), 1, []);
    
    % Backstep experiment parameters
    if isfield(TRIAL, 'time_state_dwell_on')
        SACS_ALL_TRIAL.time_state_dwell_on = ones(size(all_sac_ind_onset)) * TRIAL.time_state_dwell_on(end);
        SACS_ALL_TRIAL.time_dwell = ones(size(all_sac_ind_onset)) * TRIAL.time_dwell; 
        
        if isfield(TRIAL, 'time_state_dwell_off')
            SACS_ALL_TRIAL.time_state_dwell_off = ones(size(all_sac_ind_onset)) * TRIAL.time_state_dwell_off(end);
        else
            SACS_ALL_TRIAL.time_state_dwell_off = ones(size(all_sac_ind_onset)) * nan;
        end
    end
end
%% TAGS 1,2,3 algorithm
function SACS_ALL_TRIAL = tag_primary_saccades(SACS_ALL_TRIAL, TRIAL, threshold_pos, threshold_ang, threshold_micro_saccade)
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
end
%% TAGS 4 and 5
function SACS_ALL_TRIAL = tag_corrective_saccades(SACS_ALL_TRIAL, TRIAL, threshold_pos, threshold_ang)
    % Tag saccades as 'corr_success' or 'corr_fail' following a 'prim_success'

    idx_prim_success = find(SACS_ALL_TRIAL.tag == 1); % find primary success saccades

    if ~isempty(idx_prim_success) && (TRIAL.jump_cond == 1) && (TRIAL.task_cond == 1)
        for counter_prim = 1:length(idx_prim_success)
            idx_prim = idx_prim_success(counter_prim);
            time_last_cue_pres = TRIAL.time_state_cue_present(end);
            time_start_search = SACS_ALL_TRIAL.time_offset(idx_prim);

            if time_start_search < time_last_cue_pres
                time_finish_search = TRIAL.time_state_str_fixation( ...
                    find(TRIAL.time_state_str_fixation > time_start_search, 1, 'first'));
                flag_last_cue = false;
            else
                time_finish_search = TRIAL.time_end;
                flag_last_cue = true;
            end

            idx_sac = find((SACS_ALL_TRIAL.time_onset > time_start_search) & ...
                           (SACS_ALL_TRIAL.time_onset < time_finish_search), 1, 'first');

            if ~isempty(idx_sac)
                % Assign basic timing and visual properties
                SACS_ALL_TRIAL.time_visual(idx_sac)      = time_start_search;
                SACS_ALL_TRIAL.time_auditory(idx_sac)    = time_start_search;
                SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = TRIAL.cue_x;
                SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = TRIAL.cue_y;
                SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.end_x;
                SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.end_y;
                SACS_ALL_TRIAL.reaction(idx_sac)         = ...
                    (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;

                % Vector amplitudes and angles
                dx = SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac);
                dy = SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac);
                SACS_ALL_TRIAL.visual_amp_x(idx_sac) = dx;
                SACS_ALL_TRIAL.visual_amp_y(idx_sac) = dy;
                SACS_ALL_TRIAL.visual_amp_m(idx_sac) = sqrt(dx^2 + dy^2);
                SACS_ALL_TRIAL.visual_ang(idx_sac)   = atan2d(dy, dx);

                % Compute differences from previous prim saccade
                dx_prim = SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.eye_px_offset(idx_prim);
                dy_prim = SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.eye_py_offset(idx_prim);
                visual_amp_m_prim = sqrt(dx_prim^2 + dy_prim^2);

                % Start and end point differences
                SACS_ALL_TRIAL.diff_start(idx_sac) = ...
                    sqrt((SACS_ALL_TRIAL.eye_px_onset(idx_sac) - SACS_ALL_TRIAL.eye_px_offset(idx_prim))^2 + ...
                         (SACS_ALL_TRIAL.eye_py_onset(idx_sac) - SACS_ALL_TRIAL.eye_py_offset(idx_prim))^2);
                SACS_ALL_TRIAL.diff_finish(idx_sac) = ...
                    sqrt((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac))^2 + ...
                         (SACS_ALL_TRIAL.eye_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac))^2);

                % Angular difference between the saccade vector and expected visual vector
                dot_product = (SACS_ALL_TRIAL.eye_amp_x(idx_sac) * dx_prim) + ...
                              (SACS_ALL_TRIAL.eye_amp_y(idx_sac) * dy_prim);
                angle_diff = abs(acosd(dot_product / ...
                               (SACS_ALL_TRIAL.eye_amp_m(idx_sac) * visual_amp_m_prim)));
                SACS_ALL_TRIAL.diff_ang(idx_sac) = angle_diff;

                % Tagging
                if SACS_ALL_TRIAL.visual_amp_m(idx_sac) > eps && ...
                   SACS_ALL_TRIAL.diff_start(idx_sac) < threshold_pos && ...
                   angle_diff < threshold_ang && ...
                   SACS_ALL_TRIAL.diff_finish(idx_sac) < threshold_pos

                    SACS_ALL_TRIAL.tag(idx_sac) = 4; % 'corr_success'
                    SACS_ALL_TRIAL.flag_last_cue(idx_sac) = flag_last_cue;

                elseif SACS_ALL_TRIAL.diff_start(idx_sac) < threshold_pos
                    SACS_ALL_TRIAL.tag(idx_sac) = 5; % 'corr_fail'
                    SACS_ALL_TRIAL.flag_last_cue(idx_sac) = flag_last_cue;
                end
            end
        end
    end
end
%% BACK_TO_CENTER_SACCADES 6, 7, 8
function SACS_ALL_TRIAL = back_to_center_saccades(TRIAL, SACS_ALL_TRIAL, threshold_pos)

    for counter_cue_pres = 1 : length(TRIAL.time_state_cue_present)
        time_finish_search  = TRIAL.time_state_cue_present(counter_cue_pres);
        if counter_cue_pres==1
            time_start_search = TRIAL.time_start;
        else
            time_start_search = TRIAL.time_state_str_fixation( find(TRIAL.time_state_str_fixation < time_finish_search, 1, 'last') ) - TRIAL.time_pursuit;
        end
        idx_sac = find((SACS_ALL_TRIAL.time_onset > time_start_search) & (SACS_ALL_TRIAL.time_onset < time_finish_search), 1, 'last');
        if ~isempty(idx_sac)
            if (SACS_ALL_TRIAL.tag(idx_sac) >= 1) && (SACS_ALL_TRIAL.tag(idx_sac) <= 3)
                continue;
            end
            if idx_sac ~= 1
                SACS_ALL_TRIAL.time_visual(idx_sac) = SACS_ALL_TRIAL.time_offset(idx_sac-1);
            else
                SACS_ALL_TRIAL.time_visual(idx_sac) = TRIAL.time_state_str_pursuit(idx_sac);
            end
            SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
            SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = SACS_ALL_TRIAL.eye_px_onset(idx_sac);
            SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = SACS_ALL_TRIAL.eye_py_onset(idx_sac);
            SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.start_x;
            SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.start_y;
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

            if (SACS_ALL_TRIAL.diff_finish(idx_sac) < threshold_pos)
                SACS_ALL_TRIAL.tag(idx_sac) = 8; 
            elseif (SACS_ALL_TRIAL.diff_finish(idx_sac) >= threshold_pos)
                SACS_ALL_TRIAL.tag(idx_sac) = 10; 
            end
        end
    end

    idx_back_center = find(SACS_ALL_TRIAL.tag == 8);
    if ~isempty(idx_back_center)
        for counter_back_center = 1 : length(idx_back_center)
            idx_sac = idx_back_center(counter_back_center);
            if (idx_sac == 1) 
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
                SACS_ALL_TRIAL.tag(idx_sac) = 6;
                continue;
            end
            if (SACS_ALL_TRIAL.tag(idx_sac-1) >= 1) && (SACS_ALL_TRIAL.tag(idx_sac-1) <= 5)
                SACS_ALL_TRIAL.time_auditory(idx_sac)      = SACS_ALL_TRIAL.time_offset(idx_sac-1);
                SACS_ALL_TRIAL.time_visual(idx_sac)      = TRIAL.time_state_str_pursuit(find(TRIAL.time_state_str_pursuit < SACS_ALL_TRIAL.time_onset(idx_sac), 1, 'last'));
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
                SACS_ALL_TRIAL.tag(idx_sac) = 7;
                continue;
            end
        end
    end

    idx_prim_corr = find((SACS_ALL_TRIAL.tag >= 1) & (SACS_ALL_TRIAL.tag <= 5));
    if ~isempty(idx_prim_corr)
        for counter_prim_corr = 1 : length(idx_prim_corr)
            idx_sac = idx_prim_corr(counter_prim_corr)+1;
            if idx_sac > length(SACS_ALL_TRIAL.tag)
                break;
            end
            if ~isnan(SACS_ALL_TRIAL.tag(idx_sac))
                continue;
            end
            diff_finish = sqrt( ...
                    ((SACS_ALL_TRIAL.eye_px_offset(idx_sac) - TRIAL.start_x).^2) + ...
                    ((SACS_ALL_TRIAL.eye_py_offset(idx_sac) - TRIAL.start_y).^2) );
            if diff_finish < threshold_pos
                SACS_ALL_TRIAL.time_visual(idx_sac)      = SACS_ALL_TRIAL.time_offset(idx_sac-1);
                SACS_ALL_TRIAL.time_auditory(idx_sac)    = SACS_ALL_TRIAL.time_offset(idx_sac-1);
                if (SACS_ALL_TRIAL.tag(idx_sac-1) >= 1) && (SACS_ALL_TRIAL.tag(idx_sac-1) <= 3)
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
                SACS_ALL_TRIAL.tag(idx_sac) = 7;
            end
        end
    end
end
%% tag_irrelevant_saccades 9, 10
function SACS_ALL_TRIAL = tag_irrelevant_saccades(SACS_ALL_TRIAL, threshold_pos)
    % Compute Euclidean distance between eye and target end points
    diff_finish = sqrt( ...
        ((SACS_ALL_TRIAL.eye_px_offset - SACS_ALL_TRIAL.tgt_px_offset).^2) + ...
        ((SACS_ALL_TRIAL.eye_py_offset - SACS_ALL_TRIAL.tgt_py_offset).^2) );

    % Tag saccades as 'target_irrelev' if they are untagged and close to the target
    idx_target_irrelev = isnan(SACS_ALL_TRIAL.tag) & (diff_finish < threshold_pos);
    SACS_ALL_TRIAL.visual_px_onset( idx_target_irrelev) = SACS_ALL_TRIAL.eye_px_onset(idx_target_irrelev);
    SACS_ALL_TRIAL.visual_py_onset( idx_target_irrelev) = SACS_ALL_TRIAL.eye_py_onset(idx_target_irrelev); 
    SACS_ALL_TRIAL.visual_px_offset(idx_target_irrelev) = SACS_ALL_TRIAL.tgt_px_offset(idx_target_irrelev);
    SACS_ALL_TRIAL.visual_py_offset(idx_target_irrelev) = SACS_ALL_TRIAL.tgt_py_offset(idx_target_irrelev);
    SACS_ALL_TRIAL.diff_finish(     idx_target_irrelev) = diff_finish(idx_target_irrelev);
    SACS_ALL_TRIAL.tag(             idx_target_irrelev) = 9; % target_irrelev

    % Tag remaining untagged saccades as 'other_irrelev'
    idx_other_irrelev = isnan(SACS_ALL_TRIAL.tag);
    SACS_ALL_TRIAL.visual_px_onset( idx_other_irrelev) = SACS_ALL_TRIAL.eye_px_onset( idx_other_irrelev);
    SACS_ALL_TRIAL.visual_py_onset( idx_other_irrelev) = SACS_ALL_TRIAL.eye_py_onset( idx_other_irrelev);
    SACS_ALL_TRIAL.visual_px_offset(idx_other_irrelev) = SACS_ALL_TRIAL.eye_px_offset(idx_other_irrelev);
    SACS_ALL_TRIAL.visual_py_offset(idx_other_irrelev) = SACS_ALL_TRIAL.eye_py_offset(idx_other_irrelev);
    SACS_ALL_TRIAL.diff_finish(     idx_other_irrelev) = 0.0;
    SACS_ALL_TRIAL.tag(             idx_other_irrelev) = 10; % other_irrelev
end
%% Tag 11
function SACS_ALL_TRIAL = tag_primary_no_correction_saccades(SACS_ALL_TRIAL, TRIAL, threshold_pos)
% TAG_PRIMARY_NO_CORRECTION_SACCADES
% Identifies and tags primary saccades that did not lead to a corrective saccade.
% Adds a new saccade entry with tag 11 for these cases.

    num_saccades = length(SACS_ALL_TRIAL.validity); % number of detected saccades
    idx_prim_success = find(SACS_ALL_TRIAL.tag == 1); % indices of 'prim_success'

    if ~isempty(idx_prim_success)
        for counter_prim = 1 : length(idx_prim_success)
            idx_prim = idx_prim_success(counter_prim);
            next_idx_prim = idx_prim + 1;

            diff_finish = sqrt( ...
                ((SACS_ALL_TRIAL.tgt_px_offset(idx_prim) - SACS_ALL_TRIAL.eye_px_offset(idx_prim)).^2) + ...
                ((SACS_ALL_TRIAL.tgt_py_offset(idx_prim) - SACS_ALL_TRIAL.eye_py_offset(idx_prim)).^2) );

            if (((next_idx_prim > num_saccades) || ...
                ((next_idx_prim <= num_saccades) && ~(SACS_ALL_TRIAL.tag(next_idx_prim) == 4))) && ...
                (diff_finish > threshold_pos))

                fields_to_copy = {
                    'validity', 'trial_num', 'count', 'flag_last_cue', ...
                    'start_x', 'start_y', 'cue_x', 'cue_y', 'end_x', 'end_y', ...
                    'time_auditory', 'time_vmax', 'time_visual', 'time_offset', ...
                    'reaction', 'diff_start', 'diff_finish', 'diff_ang', ...
                    'duration', 'eye_vm_max', 'eye_px_onset', 'eye_py_onset', ...
                    'eye_px_offset', 'eye_py_offset', 'eye_amp_x', 'eye_amp_y', ...
                    'eye_amp_m', 'eye_ang', 'tgt_px_onset', 'tgt_px_offset', ...
                    'tgt_py_onset', 'tgt_py_offset', 'tgt_num', 'cue_x_high_rew', ...
                    'cue_y_high_rew', 'cue_x_low_rew', 'cue_y_low_rew', ...
                    'high_rew_tgt_num', 'low_rew_tgt_num', 'task_cond', 'choice', ...
                    'rew_cond', 'jump_cond', 'tgt_cond', 'stim_flag', ...
                    'state_stim_t_start'};

                for i = 1:length(fields_to_copy)
                    fld = fields_to_copy{i};
                    SACS_ALL_TRIAL.(fld)(end+1) = SACS_ALL_TRIAL.(fld)(idx_prim);
                end

                % Special handling for time_onset
                if next_idx_prim > num_saccades
                    SACS_ALL_TRIAL.time_onset(end+1) = TRIAL.time_end;
                else
                    SACS_ALL_TRIAL.time_onset(end+1) = SACS_ALL_TRIAL.time_onset(next_idx_prim);
                end

                % Use offset of prim as visual onset
                SACS_ALL_TRIAL.visual_px_onset(end+1) = SACS_ALL_TRIAL.eye_px_offset(idx_prim); 
                SACS_ALL_TRIAL.visual_py_onset(end+1) = SACS_ALL_TRIAL.eye_py_offset(idx_prim);
                SACS_ALL_TRIAL.visual_px_offset(end+1) = SACS_ALL_TRIAL.tgt_px_offset(idx_prim);
                SACS_ALL_TRIAL.visual_py_offset(end+1) = SACS_ALL_TRIAL.tgt_py_offset(idx_prim);

                % Compute visual amplitude and angle
                if next_idx_prim <= num_saccades
                    amp_x = SACS_ALL_TRIAL.visual_px_offset(next_idx_prim) - SACS_ALL_TRIAL.visual_px_onset(next_idx_prim);
                    amp_y = SACS_ALL_TRIAL.visual_py_offset(next_idx_prim) - SACS_ALL_TRIAL.visual_py_onset(next_idx_prim);
                else
                    amp_x = NaN;
                    amp_y = NaN;
                end
                SACS_ALL_TRIAL.visual_amp_x(end+1) = amp_x;
                SACS_ALL_TRIAL.visual_amp_y(end+1) = amp_y;
                SACS_ALL_TRIAL.visual_amp_m(end+1) = sqrt(amp_x^2 + amp_y^2);
                SACS_ALL_TRIAL.visual_ang(end+1) = atan2d(amp_y, amp_x);

                % Assign tag 11 for primary with no correction
                SACS_ALL_TRIAL.tag(end+1) = 11;

                % Optional backstep experimental fields
                if isfield(TRIAL, 'time_state_dwell_on')
                    SACS_ALL_TRIAL.time_state_dwell_on(end+1) = SACS_ALL_TRIAL.time_state_dwell_on(idx_prim);
                    SACS_ALL_TRIAL.time_dwell(end+1) = SACS_ALL_TRIAL.time_dwell(idx_prim);
                    SACS_ALL_TRIAL.time_state_dwell_off(end+1) = SACS_ALL_TRIAL.time_state_dwell_off(idx_prim);
                end
            end
        end
    end
end
%% fill_nan_visual_and_set_irrelev_reaction
function SACS_ALL_TRIAL = fill_nan_visual_and_set_irrelev_reaction(SACS_ALL_TRIAL, TRIAL)
% Fills missing visual saccade metrics and sets reaction times
% for 'target_irrelev' (tag==9) and 'other_irrelev' (tag==10) saccades.

    % Fill in NaN values for visual saccade amplitudes and angles
    idx_nan_visual_values = isnan(SACS_ALL_TRIAL.visual_amp_x);

    SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values) = ...
        SACS_ALL_TRIAL.visual_px_offset(idx_nan_visual_values) - ...
        SACS_ALL_TRIAL.visual_px_onset(idx_nan_visual_values);

    SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values) = ...
        SACS_ALL_TRIAL.visual_py_offset(idx_nan_visual_values) - ...
        SACS_ALL_TRIAL.visual_py_onset(idx_nan_visual_values);

    SACS_ALL_TRIAL.visual_amp_m(idx_nan_visual_values) = ...
        sqrt(SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values).^2 + ...
             SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values).^2);

    SACS_ALL_TRIAL.visual_ang(idx_nan_visual_values) = ...
        atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values), ...
               SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values));

    SACS_ALL_TRIAL.diff_start(idx_nan_visual_values) = ...
        sqrt((SACS_ALL_TRIAL.eye_px_onset(idx_nan_visual_values) - ...
              SACS_ALL_TRIAL.visual_px_onset(idx_nan_visual_values)).^2 + ...
             (SACS_ALL_TRIAL.eye_py_onset(idx_nan_visual_values) - ...
              SACS_ALL_TRIAL.visual_py_onset(idx_nan_visual_values)).^2);

    SACS_ALL_TRIAL.diff_ang(idx_nan_visual_values) = abs(acosd( ...
        (SACS_ALL_TRIAL.eye_amp_x(idx_nan_visual_values) .* ...
         SACS_ALL_TRIAL.visual_amp_x(idx_nan_visual_values) + ...
         SACS_ALL_TRIAL.eye_amp_y(idx_nan_visual_values) .* ...
         SACS_ALL_TRIAL.visual_amp_y(idx_nan_visual_values)) ./ ...
        (SACS_ALL_TRIAL.eye_amp_m(idx_nan_visual_values) .* ...
         SACS_ALL_TRIAL.visual_amp_m(idx_nan_visual_values)) ));

    % Set reaction time for irrelevant saccades (tag 9 and 10)
    idx_irrelev = find(SACS_ALL_TRIAL.tag == 9 | SACS_ALL_TRIAL.tag == 10);

    if ~isempty(idx_irrelev)
        for counter_target_irrelev = 1:length(idx_irrelev)
            idx_sac = idx_irrelev(counter_target_irrelev);
            if idx_sac == 1
                SACS_ALL_TRIAL.time_visual(idx_sac) = TRIAL.time_state_str_pursuit(idx_sac);
            else
                SACS_ALL_TRIAL.time_visual(idx_sac) = SACS_ALL_TRIAL.time_offset(idx_sac - 1);
            end
            SACS_ALL_TRIAL.reaction(idx_sac) = ...
                (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;
        end
    end
end
%% tag 12
function SACS_ALL_TRIAL = tag_post_correction_saccades(SACS_ALL_TRIAL, threshold_pos)
% TAG_POST_CORRECTION_SACCADES
% This function scans through all saccades tagged as 'irrelevant' (tag 9 or 10)
% and re-tags some of them as 'post-correction return saccades' (tag 12)
% based on spatial and sequential constraints relative to the previous saccade.
%

% OUTPUT:
%   - SACS_ALL_TRIAL: updated with some tag=12 assignments

% Find all saccades currently marked as irrelevant (tags 9 or 10)
idx_irrelev = find((SACS_ALL_TRIAL.tag == 9) | (SACS_ALL_TRIAL.tag == 10));

% Only proceed if such saccades exist
if ~isempty(idx_irrelev)
    for counter_irrelev = 1:length(idx_irrelev)
        idx_sac = idx_irrelev(counter_irrelev);
        prev_idx_sac = idx_sac - 1;  % get the previous saccade index

        if prev_idx_sac >= 1
            % Distance from current saccade's start to previous saccade's end
            start_pos_drift = sqrt( ...
                (SACS_ALL_TRIAL.eye_px_onset(idx_sac) - SACS_ALL_TRIAL.eye_px_offset(prev_idx_sac)).^2 + ...
                (SACS_ALL_TRIAL.eye_py_onset(idx_sac) - SACS_ALL_TRIAL.eye_py_offset(prev_idx_sac)).^2);

            % Distance from current saccade's endpoint to the target location of the previous saccade
            sac_offset_dist_to_prev_sac_tgt = sqrt( ...
                (SACS_ALL_TRIAL.eye_px_offset(idx_sac) - SACS_ALL_TRIAL.tgt_px_offset(prev_idx_sac)).^2 + ...
                (SACS_ALL_TRIAL.eye_py_offset(idx_sac) - SACS_ALL_TRIAL.tgt_py_offset(prev_idx_sac)).^2);

            % Distance between previous saccade's endpoint and its intended target
            diff_finish = sqrt( ...
                (SACS_ALL_TRIAL.eye_px_offset(prev_idx_sac) - SACS_ALL_TRIAL.tgt_px_offset(prev_idx_sac)).^2 + ...
                (SACS_ALL_TRIAL.eye_py_offset(prev_idx_sac) - SACS_ALL_TRIAL.tgt_py_offset(prev_idx_sac)).^2);

            % Apply tagging logic:
            if  (SACS_ALL_TRIAL.tag(prev_idx_sac) == 4) && ... % prev saccade was corrective and successful
                (sac_offset_dist_to_prev_sac_tgt < threshold_pos) && ... % current saccade landed near previous target
                (start_pos_drift < 0.5) && ... % eye hasn't drifted much from end of last saccade
                (diff_finish > 1.0)           % last saccade ended away from its target
                SACS_ALL_TRIAL.tag(idx_sac) = 12;
            end
        end
    end
end
end
%% fix_time_visual_with_photodiode
function sac_data = fix_time_visual_with_photodiode(sac_data, event_data)
% FIX_TIME_VISUAL_WITH_PHOTODIODE Corrects visual event times using photodiode data
%   sac_data = fix_time_visual_with_photodiode(sac_data, event_data) adjusts
%   the time_visual field in sac_data based on photodiode alignment data
%   from event_data for specific tags (1,2,3,6,7)

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
            % 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3 
            % 'back_center_success' tag 6 % 'back_center_prim' tag 7
            
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
                reaction_time = sac_onset_ind - visual_ind_converted; % 'visual_ind_converted' is in BEHAVE time domain
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
%% Tag 13
function SACS_ALL_TRIAL = tag_unreached_corrective_saccades(SACS_ALL_TRIAL, TRIAL)
% TAG_UNREACHED_CORRECTIVE_SACCADES
% Tags saccades that followed a corrective success but failed to reach the target (tag 13).
% 
% Inputs:
%   - SACS_ALL_TRIAL: struct containing all saccade info for a trial
%   - TRIAL: struct with trial metadata (e.g., end time)
%
% Output:
%   - SACS_ALL_TRIAL: updated struct with tag 13 saccades appended

    num_saccades = length(SACS_ALL_TRIAL.validity);
    idx_corr_success = find(SACS_ALL_TRIAL.tag == 4); % Tag 4: corrective success

    if ~isempty(idx_corr_success)
        for counter_corr = 1:length(idx_corr_success)
            idx_corr = idx_corr_success(counter_corr);
            next_idx_corr = idx_corr + 1;

            % Distance between eye offset and target offset (finish error)
            diff_finish = sqrt( ...
                (SACS_ALL_TRIAL.tgt_px_offset(idx_corr) - SACS_ALL_TRIAL.eye_px_offset(idx_corr))^2 + ...
                (SACS_ALL_TRIAL.tgt_py_offset(idx_corr) - SACS_ALL_TRIAL.eye_py_offset(idx_corr))^2 );

            % Condition: (last saccade or next one not a double-corrective) and miss > 1 deg
            if (((next_idx_corr > num_saccades) || ...
                (SACS_ALL_TRIAL.tag(next_idx_corr) ~= 12)) && ...
                (diff_finish > 1.0))

                % Append a new entry with tag 13 and fields copied from idx_corr
                SACS_ALL_TRIAL = append_tag13_entry(SACS_ALL_TRIAL, idx_corr, next_idx_corr, TRIAL);
            end
        end
    end
end
%%
function SACS_ALL_TRIAL = append_tag13_entry(SACS_ALL_TRIAL, idx_corr, next_idx_corr, TRIAL)
% Helper function to append a new entry tagged as 13

    fields_to_copy = {
        'validity','trial_num','count','flag_last_cue','start_x','start_y','cue_x','cue_y', ...
        'end_x','end_y','time_auditory','time_vmax','time_visual','time_offset', ...
        'visual_px_onset','visual_py_onset','visual_px_offset','visual_py_offset','reaction', ...
        'diff_start','diff_finish','diff_ang','duration','eye_vm_max','eye_px_onset','eye_py_onset', ...
        'eye_px_offset','eye_py_offset','eye_amp_x','eye_amp_y','eye_amp_m','eye_ang', ...
        'tgt_px_onset','tgt_px_offset','tgt_py_onset','tgt_py_offset','tgt_num', ...
        'cue_x_high_rew','cue_y_high_rew','cue_x_low_rew','cue_y_low_rew','high_rew_tgt_num', ...
        'low_rew_tgt_num','task_cond','choice','rew_cond','jump_cond','tgt_cond','stim_flag','state_stim_t_start'
    };

    for f = fields_to_copy
        field = f{1};
        SACS_ALL_TRIAL.(field)(end+1) = SACS_ALL_TRIAL.(field)(idx_corr);
    end

    % Time onset depends on whether there is a next saccade
    if next_idx_corr > length(SACS_ALL_TRIAL.time_onset)
        SACS_ALL_TRIAL.time_onset(end+1) = TRIAL.time_end;
    else
        SACS_ALL_TRIAL.time_onset(end+1) = SACS_ALL_TRIAL.time_onset(next_idx_corr);
    end

    % Tag the new saccade as 13
    SACS_ALL_TRIAL.tag(end+1) = 13;

    % Compute visual amplitude and angle based on eye offset → target offset
    SACS_ALL_TRIAL.visual_amp_x(end+1) = SACS_ALL_TRIAL.visual_px_offset(next_idx_corr) - SACS_ALL_TRIAL.visual_px_onset(next_idx_corr);
    SACS_ALL_TRIAL.visual_amp_y(end+1) = SACS_ALL_TRIAL.visual_py_offset(next_idx_corr) - SACS_ALL_TRIAL.visual_py_onset(next_idx_corr);
    SACS_ALL_TRIAL.visual_amp_m(end+1) = hypot(SACS_ALL_TRIAL.visual_amp_x(end), SACS_ALL_TRIAL.visual_amp_y(end));
    SACS_ALL_TRIAL.visual_ang(end+1) = atan2d(SACS_ALL_TRIAL.visual_amp_y(end), SACS_ALL_TRIAL.visual_amp_x(end));

    % Backstep dwell timing, if it exists
    if isfield(SACS_ALL_TRIAL, 'time_state_dwell_on')
        SACS_ALL_TRIAL.time_state_dwell_on(end+1) = SACS_ALL_TRIAL.time_state_dwell_on(idx_corr);
        SACS_ALL_TRIAL.time_dwell(end+1) = SACS_ALL_TRIAL.time_dwell(idx_corr);
        SACS_ALL_TRIAL.time_state_dwell_off(end+1) = SACS_ALL_TRIAL.time_state_dwell_off(idx_corr);
    end
end
%% Tag 6
function SACS_ALL_TRIAL = tag_back_to_center_success(SACS_ALL_TRIAL, TRIAL, idx_sac, threshold_pos)
% TAG_BACK_TO_CENTER_SUCCESS Tags a saccade as 'back to center success' (tag = 6) if it ends near the start position.

    % Compute distance from saccade offset to starting eye position
    diff_finish = sqrt( ...
        (SACS_ALL_TRIAL.eye_px_offset(idx_sac) - TRIAL.start_x)^2 + ...
        (SACS_ALL_TRIAL.eye_py_offset(idx_sac) - TRIAL.start_y)^2 );

    % If end of saccade is near the trial's start location
    if diff_finish < threshold_pos

        % Assign time of visual cue (based on pursuit start)
        SACS_ALL_TRIAL.time_visual(idx_sac)      = TRIAL.time_state_str_pursuit(idx_sac);

        % Set visual movement vector (onset: where eyes were, offset: target start)
        SACS_ALL_TRIAL.visual_px_onset(idx_sac)  = SACS_ALL_TRIAL.eye_px_onset(idx_sac);
        SACS_ALL_TRIAL.visual_py_onset(idx_sac)  = SACS_ALL_TRIAL.eye_py_onset(idx_sac);
        SACS_ALL_TRIAL.visual_px_offset(idx_sac) = TRIAL.start_x;
        SACS_ALL_TRIAL.visual_py_offset(idx_sac) = TRIAL.start_y;

        % Compute reaction time in ms
        SACS_ALL_TRIAL.reaction(idx_sac)         = (SACS_ALL_TRIAL.time_onset(idx_sac) - SACS_ALL_TRIAL.time_visual(idx_sac)) * 1000.0;

        % Visual vector amplitude and angle
        SACS_ALL_TRIAL.visual_amp_x(idx_sac)     = SACS_ALL_TRIAL.visual_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac);
        SACS_ALL_TRIAL.visual_amp_y(idx_sac)     = SACS_ALL_TRIAL.visual_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac);
        SACS_ALL_TRIAL.visual_amp_m(idx_sac)     = hypot(SACS_ALL_TRIAL.visual_amp_x(idx_sac), SACS_ALL_TRIAL.visual_amp_y(idx_sac));
        SACS_ALL_TRIAL.visual_ang(idx_sac)       = atan2d(SACS_ALL_TRIAL.visual_amp_y(idx_sac), SACS_ALL_TRIAL.visual_amp_x(idx_sac));

        % Difference between saccade and visual vectors (start, end, angle)
        SACS_ALL_TRIAL.diff_start(idx_sac) = hypot( ...
            SACS_ALL_TRIAL.eye_px_onset(idx_sac) - SACS_ALL_TRIAL.visual_px_onset(idx_sac), ...
            SACS_ALL_TRIAL.eye_py_onset(idx_sac) - SACS_ALL_TRIAL.visual_py_onset(idx_sac) );

        SACS_ALL_TRIAL.diff_finish(idx_sac) = hypot( ...
            SACS_ALL_TRIAL.eye_px_offset(idx_sac) - SACS_ALL_TRIAL.visual_px_offset(idx_sac), ...
            SACS_ALL_TRIAL.eye_py_offset(idx_sac) - SACS_ALL_TRIAL.visual_py_offset(idx_sac) );

        SACS_ALL_TRIAL.diff_ang(idx_sac) = abs(acosd( ...
            (SACS_ALL_TRIAL.eye_amp_x(idx_sac) * SACS_ALL_TRIAL.visual_amp_x(idx_sac) + ...
             SACS_ALL_TRIAL.eye_amp_y(idx_sac) * SACS_ALL_TRIAL.visual_amp_y(idx_sac)) / ...
            (SACS_ALL_TRIAL.eye_amp_m(idx_sac) * SACS_ALL_TRIAL.visual_amp_m(idx_sac)) ));

        % Assign tag 6: back to center success
        SACS_ALL_TRIAL.tag(idx_sac) = 6;
    end
end
%% tag_saccade_counting
function SACS_ALL_TRIAL = tag_saccade_counting(SACS_ALL_TRIAL)
%TAG_SACCADE_COUNTING Assign count values to each category of tagged saccades
%
% This function sets the 'count' field for each group of tagged saccades
% (e.g., primary, corrective, back-to-center, irrelevant) in the structure SACS_ALL_TRIAL.

    % Tag 1: prim_success, Tag 2: prim_attempt
    prim_tags = (SACS_ALL_TRIAL.tag == 1) | (SACS_ALL_TRIAL.tag == 2);
    if sum(prim_tags) > 0
        SACS_ALL_TRIAL.count(prim_tags) = 1 : sum(prim_tags);
    end

    % Tag 11: prim_no_corr
    prim_no_corr_tags = SACS_ALL_TRIAL.tag == 11;
    if sum(prim_no_corr_tags) > 0
        SACS_ALL_TRIAL.count(prim_no_corr_tags) = 1 : sum(prim_no_corr_tags);
    end

    % Tag 4: corr_success, Tag 5: corr_fail
    corr_tags = (SACS_ALL_TRIAL.tag == 4) | (SACS_ALL_TRIAL.tag == 5);
    if sum(corr_tags) > 0
        SACS_ALL_TRIAL.count(corr_tags) = 1 : sum(corr_tags);
    end

    % Tag 12: db_corr_sac
    db_corr_tags = SACS_ALL_TRIAL.tag == 12;
    if sum(db_corr_tags) > 0
        SACS_ALL_TRIAL.count(db_corr_tags) = 1 : sum(db_corr_tags);
    end

    % Tag 6: back_center_success, Tag 7: back_center_prim
    back_center_tags = (SACS_ALL_TRIAL.tag == 6) | (SACS_ALL_TRIAL.tag == 7);
    if sum(back_center_tags) > 0
        SACS_ALL_TRIAL.count(back_center_tags) = 1 : sum(back_center_tags);
    end

    % Tag 3: prim_fail, Tag 8: back_center_irrelev,
    % Tag 9: target_irrelev, Tag 10: other_irrelev
    irrelev_tags = (SACS_ALL_TRIAL.tag == 3) | ...
                   (SACS_ALL_TRIAL.tag == 8) | ...
                   (SACS_ALL_TRIAL.tag == 9) | ...
                   (SACS_ALL_TRIAL.tag == 10);
    if sum(irrelev_tags) > 0
        SACS_ALL_TRIAL.count(irrelev_tags) = 1 : sum(irrelev_tags);
    end
end
%% Tag 14
function SACS_ALL_TRIAL = tag_virtual_visual_event_after_irrelevant_saccade(SACS_ALL_TRIAL)
% TAG_VIRTUAL_VISUAL_EVENT_AFTER_IRRELEVANT_SACCADE - Creates synthetic visual events (tag 14)
% after 'other irrelevant' saccades (tag 10), using the end of the irrelevant saccade and
% the beginning of the next saccade to define a visual transition.
%
% INPUT:
%   SACS_ALL_TRIAL : struct with all saccade trial fields (e.g., eye data, target info, tags)
%
% OUTPUT:
%   SACS_ALL_TRIAL : same struct with new rows appended (tag 14 events)

    SACS_fields = fieldnames(SACS_ALL_TRIAL);
    idx_other_irrelev = find(SACS_ALL_TRIAL.tag == 10);

    if ~isempty(idx_other_irrelev)
        for counter_other_irrelev = 1 : length(idx_other_irrelev)
            sac_idx = idx_other_irrelev(counter_other_irrelev);
            next_sac_idx = sac_idx + 1;

            % 1. Duplicate all fields of the irrelevant saccade
            for counter_field = 1 : length(SACS_fields)
                field_name = SACS_fields{counter_field};
                SACS_ALL_TRIAL.(field_name)(end+1) = SACS_ALL_TRIAL.(field_name)(sac_idx);
            end

            % 2. Modify the duplicated entry to reflect a synthetic visual event (tag 14)
            new_tag_idx = length(SACS_ALL_TRIAL.tag);
            SACS_ALL_TRIAL.tag(new_tag_idx) = 14;

            % Set visual timing from current offset to next saccade onset
            SACS_ALL_TRIAL.time_visual(new_tag_idx) = SACS_ALL_TRIAL.time_offset(sac_idx);
            SACS_ALL_TRIAL.time_onset(new_tag_idx) = SACS_ALL_TRIAL.time_onset(next_sac_idx);

            % Set synthetic visual target as transition from eye to target
            SACS_ALL_TRIAL.visual_px_offset(new_tag_idx) = SACS_ALL_TRIAL.tgt_px_offset(sac_idx);
            SACS_ALL_TRIAL.visual_py_offset(new_tag_idx) = SACS_ALL_TRIAL.tgt_py_offset(sac_idx);
            SACS_ALL_TRIAL.visual_px_onset(new_tag_idx)  = SACS_ALL_TRIAL.eye_px_offset(sac_idx);
            SACS_ALL_TRIAL.visual_py_onset(new_tag_idx)  = SACS_ALL_TRIAL.eye_py_offset(sac_idx);

            % Compute visual vector amplitude and angle
            dx = SACS_ALL_TRIAL.visual_px_offset(new_tag_idx) - SACS_ALL_TRIAL.visual_px_onset(new_tag_idx);
            dy = SACS_ALL_TRIAL.visual_py_offset(new_tag_idx) - SACS_ALL_TRIAL.visual_py_onset(new_tag_idx);
            SACS_ALL_TRIAL.visual_amp_x(new_tag_idx) = dx;
            SACS_ALL_TRIAL.visual_amp_y(new_tag_idx) = dy;
            SACS_ALL_TRIAL.visual_amp_m(new_tag_idx) = sqrt(dx^2 + dy^2);
            SACS_ALL_TRIAL.visual_ang(new_tag_idx)   = atan2d(dy, dx);

            % Use the eye position after the saccade as the new onset point
            SACS_ALL_TRIAL.eye_px_onset(new_tag_idx) = SACS_ALL_TRIAL.eye_px_offset(sac_idx);
            SACS_ALL_TRIAL.eye_py_onset(new_tag_idx) = SACS_ALL_TRIAL.eye_py_offset(sac_idx);

            % Reaction time: time until next saccade
            SACS_ALL_TRIAL.reaction(new_tag_idx) = ...
                (SACS_ALL_TRIAL.time_onset(next_sac_idx) - SACS_ALL_TRIAL.time_offset(sac_idx)) * 1000;

            % Spatial errors from target at start and end
            SACS_ALL_TRIAL.diff_start(new_tag_idx) = sqrt( ...
                (SACS_ALL_TRIAL.tgt_px_onset(sac_idx) - SACS_ALL_TRIAL.eye_px_onset(sac_idx))^2 + ...
                (SACS_ALL_TRIAL.tgt_py_onset(sac_idx) - SACS_ALL_TRIAL.eye_py_onset(sac_idx))^2 );
            SACS_ALL_TRIAL.diff_finish(new_tag_idx) = sqrt( ...
                (SACS_ALL_TRIAL.tgt_px_offset(sac_idx) - SACS_ALL_TRIAL.eye_px_offset(sac_idx))^2 + ...
                (SACS_ALL_TRIAL.tgt_py_offset(sac_idx) - SACS_ALL_TRIAL.eye_py_offset(sac_idx))^2 );

            % Eye movement fields set to zero (since this is a virtual, not real, saccade)
            SACS_ALL_TRIAL.diff_ang(new_tag_idx) = 0;
            SACS_ALL_TRIAL.eye_amp_x(new_tag_idx) = 0;
            SACS_ALL_TRIAL.eye_amp_y(new_tag_idx) = 0;
            SACS_ALL_TRIAL.eye_amp_m(new_tag_idx) = 0;
            SACS_ALL_TRIAL.eye_ang(new_tag_idx)   = 0;
        end
    end
end
%% Tag 15
function SACS_ALL_TRIAL = tag_back_to_center_virtual_event(SACS_ALL_TRIAL, TRIAL)
%TAG_BACK_TO_CENTER_VIRTUAL_EVENT Adds virtual saccades tagged as 15.
%
% This function finds all real saccades tagged as 'back to center irrelevant' (tag 8),
% and for each one, it creates a new **virtual visual event** (tag 15).
% This virtual saccade represents a **predicted visual consequence** of the return-to-center movement.
%
% INPUTS:
%   - SACS_ALL_TRIAL: struct containing all saccade data
%   - TRIAL: struct with trial metadata, including starting fixation point (start_x, start_y)
%
% OUTPUT:
%   - SACS_ALL_TRIAL: updated with new entries for each tag 15 virtual saccade

    SACS_fields = fieldnames(SACS_ALL_TRIAL);
    idx_back_center_irrelev = find(SACS_ALL_TRIAL.tag == 8); % tag 8: real back-to-center saccade

    if ~isempty(idx_back_center_irrelev)
        for counter = 1 : length(idx_back_center_irrelev)
            sac_idx = idx_back_center_irrelev(counter);
            prev_sac_idx = sac_idx - 1;

            %--- Copy the original saccade's data as a new row ---
            for f = 1 : length(SACS_fields)
                field_name = SACS_fields{f};
                SACS_ALL_TRIAL.(field_name)(end+1) = SACS_ALL_TRIAL.(field_name)(sac_idx);
            end
            new_idx = length(SACS_ALL_TRIAL.tag); % new row index

            %--- Modify the new row to represent the virtual event (tag 15) ---
            SACS_ALL_TRIAL.tag(new_idx) = 15;

            % Virtual visual onset = where the eye landed before the return-to-center
            SACS_ALL_TRIAL.visual_px_onset(new_idx)  = SACS_ALL_TRIAL.eye_px_offset(prev_sac_idx);
            SACS_ALL_TRIAL.visual_py_onset(new_idx)  = SACS_ALL_TRIAL.eye_py_offset(prev_sac_idx);

            % Virtual visual offset = fixation center of the trial
            SACS_ALL_TRIAL.visual_px_offset(new_idx) = TRIAL.start_x;
            SACS_ALL_TRIAL.visual_py_offset(new_idx) = TRIAL.start_y;

            % Compute visual amplitude and angle of the virtual vector
            dx = SACS_ALL_TRIAL.visual_px_offset(new_idx) - SACS_ALL_TRIAL.visual_px_onset(new_idx);
            dy = SACS_ALL_TRIAL.visual_py_offset(new_idx) - SACS_ALL_TRIAL.visual_py_onset(new_idx);
            SACS_ALL_TRIAL.visual_amp_x(new_idx) = dx;
            SACS_ALL_TRIAL.visual_amp_y(new_idx) = dy;
            SACS_ALL_TRIAL.visual_amp_m(new_idx) = sqrt(dx^2 + dy^2);
            SACS_ALL_TRIAL.visual_ang(new_idx)   = atan2d(dy, dx);

            % Eye onset location = same as where eye previously landed
            SACS_ALL_TRIAL.eye_px_onset(new_idx) = SACS_ALL_TRIAL.eye_px_offset(prev_sac_idx);
            SACS_ALL_TRIAL.eye_py_onset(new_idx) = SACS_ALL_TRIAL.eye_py_offset(prev_sac_idx);

            % diff_start: distance between eye and target at start (copied from previous real saccade)
            SACS_ALL_TRIAL.diff_start(new_idx) = sqrt( ...
                (SACS_ALL_TRIAL.tgt_px_onset(sac_idx) - SACS_ALL_TRIAL.eye_px_onset(sac_idx))^2 + ...
                (SACS_ALL_TRIAL.tgt_py_onset(sac_idx) - SACS_ALL_TRIAL.eye_py_onset(sac_idx))^2 );

            % Eye amplitude and angle placeholders (0 because this is a virtual event)
            SACS_ALL_TRIAL.eye_amp_x(new_idx) = 0;
            SACS_ALL_TRIAL.eye_amp_y(new_idx) = 0;
            SACS_ALL_TRIAL.eye_amp_m(new_idx) = 0;
            SACS_ALL_TRIAL.eye_ang(new_idx)   = 0;
            SACS_ALL_TRIAL.diff_ang(new_idx)  = 0;
        end
    end
end
%% Plot debugging
function plot_trial_debug(flag_plot_trial, counter_trial, TRIAL, SACS_ALL_TRIAL, time_device, meta_data)
    if ~flag_plot_trial
        return;
    end

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
    dist_from_start_tgt  = sqrt((TRIAL.start_x - eye_px_filt).^2 + (TRIAL.start_y - eye_py_filt).^2);

    hAx_(1) = subplot(2,1,1); hold on;
    p_eye_px = plot(time_device, eye_px_filt, '-r', 'LineWidth',1);
    p_eye_py = plot(time_device, eye_py_filt, '-b', 'LineWidth',1);
    plot(time_device, tgt_px, '--r');
    plot(time_device, tgt_py, '--b');
    
   if  TRIAL.choice == 1
        p_dist_to_cue = plot(time_device, dist_to_cue_high_rew, 'k', 'LineWidth', 2);
    else
        p_dist_to_cue = plot(time_device, dist_to_cue_low_rew, 'k', 'LineWidth', 2);
    end

    p_dist_to_end = plot(time_device, dist_to_end, 'm', 'LineWidth',2);
    p_dist_from_start_tgt = plot(time_device, dist_from_start_tgt, 'g', 'LineWidth',0.5);
    ylabel('pos (deg)');
    xlabel('time (ms)');

    hAx_(2) = subplot(2,1,2); hold on;
    ylabel('angle diff (deg)');

    tgt_vec = [TRIAL.cue_x - TRIAL.start_x, TRIAL.cue_y - TRIAL.start_y];
    unit_tgt_vec = tgt_vec / norm(tgt_vec);
    
    sac_vec_p = [eye_px_filt - TRIAL.start_x, eye_py_filt - TRIAL.start_y];
    unit_sac_p = sac_vec_p ./ vecnorm(sac_vec_p')';

    sac_vec_v = [eye_vx_filt, eye_vy_filt];
    unit_sac_v = sac_vec_v ./ vecnorm(sac_vec_v')';

    ang_diff_p = acosd(unit_sac_p * unit_tgt_vec');
    ang_diff_v = acosd(unit_sac_v * unit_tgt_vec');

    p_ang_diff_p = plot(time_device, ang_diff_p, '-r', 'LineWidth',0.5);
    p_ang_diff_v = plot(time_device, ang_diff_v, '-b', 'LineWidth',0.25);

    xline(TRIAL.time_state_cue_present,'-', 'cue present');
    xline(TRIAL.time_state_sac_onset,'-', 'sac start');

    for i = 1:length(SACS_ALL_TRIAL.validity)
        tag_ = SACS_ALL_TRIAL.tag(i);
        if tag_ == 14, continue; end

        color_ = get_tag_color(tag_);
        label_ = [meta_data.sac_tag_list{tag_} '_' num2str(SACS_ALL_TRIAL.count(i))];

        if ismember(tag_, [11,13])
            x_1 = xline(hAx_(1), SACS_ALL_TRIAL.time_offset(i), '-', 'LineWidth',2, 'Color',color_);
            x_2 = xline(hAx_(2), SACS_ALL_TRIAL.time_offset(i), '-', 'LineWidth',2, 'Color',color_);
            legend_text = meta_data.sac_tag_list{tag_};
            legend(x_1, legend_text); legend(hAx_(1), 'boxoff');
            legend(x_2, legend_text); legend(hAx_(2), 'boxoff');
        else
            xline(hAx_(1), SACS_ALL_TRIAL.time_onset(i), '-', label_, 'Color', color_, 'Interpreter', 'none', 'FontSize',12);
            xline(hAx_(2), SACS_ALL_TRIAL.time_onset(i), '-', label_, 'Color', color_, 'Interpreter', 'none', 'FontSize',12);
        end
    end

    subplot(2,1,1);
    legend([p_eye_px, p_eye_py, p_dist_to_cue, p_dist_to_end, p_dist_from_start_tgt], ...
        {'x','y','dist to cue','dist to end','dist from start'});
    legend('boxoff');

    subplot(2,1,2);
    legend([p_ang_diff_p, p_ang_diff_v], {'p','v'}); legend('boxoff');

    sgtitle_str = sprintf('%s, %s', ...
        ternary(TRIAL.jump_cond == 1, 'jump', 'no jump'), ...
        ternary(TRIAL.task_cond == 1, 'forced', 'choice'));
    sgtitle(sgtitle_str);

    linkaxes(hAx_, 'x');
    ESN_Beautify_Plot(hFig_, [16,9], 9);
    pause;
end
%% plot colors
function color_ = get_tag_color(tag_)
    switch tag_
        case 1,  color_ = 'r';
        case 2,  color_ = [0.6350, 0.0780, 0.1840];
        case 3,  color_ = [0.4, 0.4, 0.4];
        case 4,  color_ = 'b';
        case 5,  color_ = [0, 0.4470, 0.7410];
        case {6,7}, color_ = [0.4940, 0.1840, 0.5560];
        case {8,9,10}, color_ = [0.6, 0.6, 0.6];
        case 11, color_ = [0.4660, 0.6740, 0.1880];
        case 13, color_ = [0.3010, 0.7450, 0.9330];
        otherwise, color_ = [0.6, 0.6, 0.6];
    end
end
%% ternary
function out = ternary(cond, val_true, val_false)
    if cond
        out = val_true;
    else
        out = val_false;
    end
end
%%