%% SET PARENT PATH
data_path = 'Y:\Ephys'; 
save_path = 'C:\Users\Jafar\Documents\reward';
JSP_global_variables(1)
global event_type_list amp_edges vel_edges length_trace tag_name_list min_rxn_time max_rxn_time ang_values ang_edges inds_span
%% Housekeeping variables
animal_name          = {'132F'};
num_animal           = length(animal_name);
tag_local_list       = [1,10];
tag_list             = {[1,4,6],[10]};
local_amp_edges      = [0.5,10]; % consider only saccades of this amp. range
task_cond_list       = [1,0];    % 1:force, 0:choice
tgt_cond_list        = [1,0];    % 1:high, 0:low
rew_cond_list        = [1,0];    % 1:high, 0:low
jump_cond_list       = [1,0];    % 1:yes, 0:No
tag_1_amp_threshold  = 3.0;        % sacccade apmplitude threshold
tag_4_amp_threshold  = 2.0;
local_RT_edges       = 50:10:600;
tag_fit_idx          = 1; % 1 for tag 1,4,6 and 10 for tag 10
win_size             = 40; % window size for tag 4 saccade
file_path = fullfile(save_path,'population_figs','behavior');
%% PERFORMANCE

% init variable
performance = struct;
sac_num     = struct;

% Loop animals
month_dir = dir2(fullfile(data_path,sprintf('data_%s','132F')));  
month_dir = month_dir([month_dir.isdir]);
counter_day = 0; % consider each day as a session
  % Loop months
 num_month = [16 17 18];
 for idx_m = num_month 
     month_path = fullfile(month_dir(idx_m).folder,month_dir(idx_m).name);
     day_dir = dir2(month_path);
     day_dir = day_dir([day_dir.isdir]);
    % Loop days
    num_day = size(day_dir,1);

    for d_idx = 1:num_day
        counter_day = counter_day + 1;
        rec_path = fullfile(day_dir(d_idx).folder, day_dir(d_idx).name);
        rec_dir = dir2(rec_path);
        rec_dir = rec_dir([rec_dir.isdir]);
        % Loop recs
        num_rec = size(rec_dir,1);
        temp_performance = nan(num_rec,1);
        tem_sac          = nan(num_rec,1);
      
        for r_idx =1:num_rec
            fprintf('loading session: %s.........\n', rec_dir(r_idx).name);
    
            data_dir = dir2(fullfile(rec_dir(r_idx).folder,rec_dir(r_idx).name, 'raw_data','*.mat*'));
            % Load data; if not available, skip the rec
            if isempty(data_dir)
               continue
            end
            loaded_data = load(fullfile(data_dir.folder, data_dir.name));
            data_fieldnames = fieldnames(loaded_data.data);
            num_trials = sum(cell2mat(arrayfun(@(x) strfind(x,'trial_'),data_fieldnames))==1);
            
            tem_sac(r_idx) = num_trials;

            chose_tgt = zeros(num_trials,1);
            task_cond = cell(num_trials,1);

            for jj = 1:num_trials
            %     rew_cond   = strcat(data.(sprintf('trial_%d',jj)).rew_cond);
                choice_cnd = strcat(loaded_data.data.(sprintf('trial_%d',jj)).choice);
                task_cond{jj,1}  = strcat(strcat(loaded_data.data.(sprintf('trial_%d',jj)).task_cond));
            
                if strcmp(choice_cnd,'h')
                    chose_tgt(jj) = 1;
                end
            end
            
            sum_choice_trials = sum(cell2mat(cellfun(@(x) strcmp(x,'choice'),task_cond,'UniformOutput',false)));
            temp_performance(r_idx) = sum(chose_tgt)/sum_choice_trials;
        end
         performance.(sprintf('month_%d',idx_m)).(sprintf('day_%d',d_idx)) = mean(temp_performance,'omitnan');
         sac_num.(sprintf('month_%d',idx_m)).(sprintf('day_%d',d_idx))     = sum(tem_sac,"omitnan");
    end
 end


day_performance =  [cell2mat(struct2cell(performance.month_16));cell2mat(struct2cell(performance.month_17));cell2mat(struct2cell(performance.month_18))];
day_performance(isnan(day_performance)) = [];

num_saccades = [cell2mat(struct2cell(sac_num.month_16));cell2mat(struct2cell(sac_num.month_17));cell2mat(struct2cell(sac_num.month_18))];
num_saccades(num_saccades==0) = [];

fig = figure;
bar(1:length(day_performance),day_performance);
xlabel('Days')
ylabel('Performance of choice trials (%)') 
ESN_Beautify_Plot(fig, [20, 8], 12);

file_path = fullfile(save_path,'population_figs','behavior');
file_name = 'Performance';
saveas(fig,fullfile(file_path,file_name), 'pdf');

% number of saccades
fig = figure;
h = bar(1:length(num_saccades),num_saccades);
h.FaceColor = [0.4940 0.1840 0.5560];
xlabel('Days')
ylabel('Number of saccades') 
ESN_Beautify_Plot(fig, [20, 8], 12);

file_path = fullfile(save_path,'population_figs','behavior');
file_name = 'Saccasde_number';
saveas(fig,fullfile(file_path,file_name), 'pdf');
%% FITING Logarithmic model

% Separate via animal 
% Go thru behavior data at rec-level


% Init. data
amp = struct;
vel = struct;
for counter_animal = 1 : num_animal
    for tag_idx = 1:length(tag_local_list)
        amp.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx)) = [];
        vel.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx)) = [];
    end
end
% Loop animals
for counter_animal = 1 : num_animal
    month_dir = dir2(fullfile(data_path,sprintf('data_%s',animal_name{counter_animal})));
    % Loop months
    num_month = size(month_dir,1);
    for m_idx =  [17 18 19]
        month_path = fullfile(month_dir(m_idx).folder,month_dir(m_idx).name);
        day_dir = dir2(month_path);
       
        % Loop days
        num_day = size(day_dir,1);
        for d_idx = 1:num_day
            rec_path = fullfile(day_dir(d_idx).folder, day_dir(d_idx).name);
            rec_dir = dir2(rec_path);
            rec_dir = rec_dir([rec_dir.isdir],:);
            % Loop recs
            num_rec = size(rec_dir,1)-1;
            for r_idx = 1:num_rec
                
                data_dir = dir2(fullfile(rec_dir(r_idx).folder,rec_dir(r_idx).name, 'analyzed_data\behavior_data\eye','*ANALYZED*'));
                % Load data; if not available, skip the rec
                if isempty(data_dir)
                    continue
                end
                loaded_data = load(fullfile(data_dir.folder, data_dir.name));
                SACS_ALL_DATA = loaded_data.sac_data;
                SACS_amp_bin  = discretize(SACS_ALL_DATA.eye_amp_m,  local_amp_edges);
                % If var. not available, skip the rec
                try
                    idx_validity = SACS_ALL_DATA.validity;
%                     idx_fix      = SACS_ALL_DATA.fix_validity(1,:) & SACS_ALL_DATA.fix_validity(2,:);
                    idx_reaction = (SACS_ALL_DATA.reaction > min_rxn_time) & (SACS_ALL_DATA.reaction < max_rxn_time);
                    idx_amp       = SACS_amp_bin == 1;
                catch
                    continue
                end
                fprintf('%s: %s\n', animal_name{counter_animal}, rec_dir(r_idx).name);
                % Loop over tags
               
                for tag_idx = 1:length(tag_local_list)
                    idx_tag = ismember(SACS_ALL_DATA.tag , tag_list{tag_idx});
                    if ~ismember( tag_list{tag_idx},10)
                        idx_tuned = idx_tag & idx_amp & idx_validity & idx_reaction;
                    else 
                        idx_tuned = idx_tag & idx_amp & idx_validity;
                    end
                    
%                     % Check fixation validity
%                     if ismember(tag_idx,[6,10])
%                         idx_tuned = idx_tuned & idx_fix;
%                     end 
                    temp_amp = SACS_ALL_DATA.eye_amp_m(idx_tuned);
                    temp_vel = SACS_ALL_DATA.eye_vm_max(idx_tuned);
                    % Concatenate data
                    amp.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx)) = [amp.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx)),temp_amp ];
                    vel.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx)) = [vel.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx)),temp_vel ];
                end
            end
        end
    end
end
% Fit eqn.
clearvars fit_obj
for counter_animal = 1 : num_animal
    for tag_idx = 1:length(tag_local_list)
        fit_obj.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx)) = fitlm(log10(amp.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx))),log10(vel.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx))));
    end
end

% Plot
fig = figure;

num_row_fig = 1;
num_col_fig = 3;

x = log10(0.5:0.01:10)';

for counter_animal = 1 : num_animal
    amp.(sprintf('m_%s', animal_name{counter_animal})).tag_tgt = [];
    vel.(sprintf('m_%s', animal_name{counter_animal})).tag_tgt = [];
    counter_tag = 1;
    for tag_idx = 1:length(tag_local_list)
        subplot(num_row_fig, num_col_fig,(counter_animal-1)*num_col_fig+counter_tag);
        scatter(log10(amp.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx))),log10(vel.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx))),...
            5,'k');
        hold on; 
        y = predict(fit_obj.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx)),x);
        if counter_animal == 1
            color_ = 'r';
        else
            color_ = 'b';
        end
        plot(x,y,color_,'LineWidth',1.5);
    
        if counter_animal == 1 && counter_tag == 1
            xlabel('Amp (deg)');
        elseif counter_animal == 1 && counter_tag == 2
            xlabel('Amp (deg)');
        end
        if counter_tag == 1
            ylabel_ = {animal_name{counter_animal},'Velocity (deg/s)'};
            ylabel(ylabel_);
        end
        
        if counter_animal == 1 && tag_idx==1
            title_ = {sprintf('tag %d,%d,%d', 1,4,6)};
        
        elseif counter_animal == 1 && tag_idx==2
             tag_idx_label = 10;
            title_ = {sprintf('tag %d', tag_idx_label)};
        else
            title_ = sprintf('alpha = %.3f, beta = %.3f', fit_params(1),fit_params(2));
        end
        
        title(title_);

        if counter_animal == 2
            subplot(num_row_fig, num_col_fig, counter_tag);
            plot(x,y,color_,'LineWidth',1.5);
        end

        counter_tag = counter_tag + 1;

        if ismember(tag_idx,[1,4,6])
            amp.(sprintf('m_%s', animal_name{counter_animal})).tag_tgt = [amp.(sprintf('m_%s', animal_name{counter_animal})).tag_tgt, amp.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx))];
            vel.(sprintf('m_%s', animal_name{counter_animal})).tag_tgt = [vel.(sprintf('m_%s', animal_name{counter_animal})).tag_tgt, vel.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',tag_idx))];
        end
    end
end

% Targeted tags
for counter_animal = 1 : num_animal
    fit_obj.(sprintf('m_%s', animal_name{counter_animal})).tag_tgt = fitlm(log10(amp.(sprintf('m_%s', animal_name{counter_animal})).tag_1'),log10(vel.(sprintf('m_%s', animal_name{counter_animal})).tag_1'));
    subplot(num_row_fig, num_col_fig,(counter_animal-1)*num_col_fig + 3);
    hold on;
    fit_params = fit_obj.(sprintf('m_%s', animal_name{counter_animal})).tag_1.Coefficients{:,1};
   
    y = predict(fit_obj.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',1)),x);
    if counter_animal == 1
        color_ = 'r';
    else
        color_ = 'b';
    end
    p_1 = plot(x,y,color_,'LineWidth',1.5);
    title_ = {'targeted saccades', sprintf('b = %.3f, a = %.3f', fit_params(1),fit_params(2))};
    
    % Other saccades
%     fit_params = fit_obj.(sprintf('m_%s', animal_name{counter_animal})).tag_2.Coefficients{:,1};
    y = predict(fit_obj.(sprintf('m_%s', animal_name{counter_animal})).(sprintf('tag_%d',2)),x);
    p_2 = plot(x,y,'k','LineWidth',1.5);

    title(title_);

    legend([p_1,p_2],{'tgt saccades','others'});
    

end
sgtitle('fit eqn: y =  10 ^ a*Log(x) + b')
ESN_Beautify_Plot(fig, [20, 8], 12);

file_path = fullfile(save_path,'population_figs','behavior');
file_name = 'saccade_vigor';
if ~exist(file_path, 'dir')
    mkdir(file_path);
end
saveas(fig,fullfile(file_path,file_name), 'pdf');
% Save the models
file_name = 'vigor_model';
% save(fullfile(data_path,'data_population',sprintf('%s.mat',file_name)), 'fit_obj', '-v7.3');
save(fullfile(save_path,'data_population',sprintf('%s.mat',file_name)), 'fit_obj', '-v7.3');
%% VIGOR VS. ACURACY PRIMARY SACCADE
% 1. End point error magnitude
% 2. Generalized variance (determinant of covariance matrix)
tag_idx   = 1; % only applicable for primary saccade
task_id   = 1; % for forced trials
debug_fig = false;
% Load vigor model
loaded_data = load('C:\Users\Jafar\Documents\reward\data_population\vigor_model.mat');
fit_obj     = loaded_data.fit_obj.(sprintf('m_%s', '132F')).(sprintf('tag_%d',tag_idx));

% Init. data
vigor            = struct;
end_error        = struct;
px_offset        = struct;
py_offset        = struct;
reaction_times   = struct;
eye_vm_vel_onset = struct;
eye_m_vel_visual = struct;
eye_m_vmax       = struct;

% Loop animals
month_dir = dir2(fullfile(data_path,sprintf('data_%s','132F')));
    
% Get vigor model per animal
% fit_params = fit_obj.(sprintf('m_%s', '132F')).tag_1.Coefficients{:,1};
% fit_params = coeffvalues(fit_obj.(sprintf('m_%s', '132F')).(sprintf('tag_%d', tag_fit_idx)));
    
counter_day = 0; % consider each day as a session
  % Loop months
 num_month = size(month_dir,1);
 for m_idx =  17 %1:num_month
     month_path = fullfile(month_dir(m_idx).folder,month_dir(m_idx).name);
     day_dir = dir2(month_path);
    % Loop days
     num_day = size(day_dir,1);
    for d_idx = 1:num_day
        counter_day = counter_day + 1;
        rec_path = fullfile(day_dir(d_idx).folder, day_dir(d_idx).name);
        rec_dir = dir2(rec_path);
        rec_dir = rec_dir([rec_dir.isdir]);
        % Loop recs
        num_rec = size(rec_dir,1);
        for r_idx =1:num_rec
            fprintf('loading session: %s.........\n', rec_dir(r_idx).name);
    
            data_dir = dir2(fullfile(rec_dir(r_idx).folder,rec_dir(r_idx).name, 'analyzed_data\behavior_data\eye','*ANALYZED*'));
            % Load data; if not available, skip the rec
            if isempty(data_dir)
               continue
            end
            loaded_data = load(fullfile(data_dir.folder, data_dir.name));
            SACS_ALL_DATA = loaded_data.sac_data;
            SACS_amp_bin  = discretize(SACS_ALL_DATA.eye_amp_m,  local_amp_edges);
            
            for counter_tgt = 1:length(tgt_cond_list)
                tgt_idx = tgt_cond_list(counter_tgt);

             fprintf('Target %d...\n',  tgt_idx);
             % Init. var
             vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             eye_m_vel_visual.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];

             idx_validity  = SACS_ALL_DATA.validity;
             idx_reaction  = (SACS_ALL_DATA.reaction > min_rxn_time) & (SACS_ALL_DATA.reaction < max_rxn_time);
             idx_tag       = SACS_ALL_DATA.tag == tag_idx;
             idx_task_cond = SACS_ALL_DATA.task_cond == task_id;
             idx_tgt       = SACS_ALL_DATA.tgt_cond == tgt_cond_list(counter_tgt);

             if ismember(tag_idx,[1,6])
                idx_amp  = SACS_ALL_DATA.eye_amp_m > tag_1_amp_threshold;
             elseif tag_idx == 4
                 idx_amp = SACS_ALL_DATA.eye_amp_m > tag_4_amp_threshold;
             else
                 idx_amp = SACS_ALL_DATA.eye_amp_m;
             end
            idx_tuned     = idx_tag & idx_validity & idx_reaction & idx_task_cond & idx_amp & idx_tgt;
            % Get relevant saccade variables; easier to debug this way
            eye_amp_m = SACS_ALL_DATA.eye_amp_m(idx_tuned);
            eye_vm_max = SACS_ALL_DATA.eye_vm_max(idx_tuned);
            visual_px_offset = SACS_ALL_DATA.cue_x(idx_tuned) - SACS_ALL_DATA.start_x(idx_tuned); % centering saccade
            visual_py_offset = SACS_ALL_DATA.cue_y(idx_tuned) - SACS_ALL_DATA.start_y(idx_tuned); % centering saccade
            eye_px_onset   = 0; % centering saccade
            eye_py_onset   = 0; % centering saccade
            eye_px_offset = SACS_ALL_DATA.eye_px_offset(idx_tuned) - SACS_ALL_DATA.start_x(idx_tuned); % centering saccade
            eye_py_offset = SACS_ALL_DATA.eye_py_offset(idx_tuned) - SACS_ALL_DATA.start_y(idx_tuned); % centering saccade
            % Compute vigor
            predict_vel   = predict(fit_obj,log10(eye_amp_m'));
            vigor_temp = log10(eye_vm_max')./predict_vel;
            % Compute saccade offset error
            offset_error = sqrt((SACS_ALL_DATA.cue_x(idx_tuned) - SACS_ALL_DATA.eye_px_offset(idx_tuned)).^2 +...
            (SACS_ALL_DATA.cue_y(idx_tuned) - SACS_ALL_DATA.eye_py_offset(idx_tuned)).^2);
        
            % Compute nominal direction of saccades (deg)
            delta_x = visual_px_offset - eye_px_onset;
            delta_y = visual_py_offset - eye_py_onset;
            visual_ang      = wrapTo360(atan2d(delta_y, delta_x));
            visual_ang_list = unique(visual_ang);
            
            % Compute nominal amp. of saccade & scaling factor to
            % normalize amp.
            nominal_amp = sqrt(delta_x.^2 + delta_y.^2);
%                     amp_scaling = norm_amp./nominal_amp;
            amp_scaling = 1;
            % Compute rotational matrix to rotate saccade vector to
            % right
            rot_eye_px_offset = zeros(1, length(visual_ang));
            rot_eye_py_offset = zeros(1, length(visual_ang));

            for visual_ang_temp = visual_ang_list
                idx_ang = visual_ang == visual_ang_temp;

                R = [cosd(-visual_ang_temp), -sind(-visual_ang_temp);
                     sind(-visual_ang_temp), cosd(-visual_ang_temp)];
              rot_eye_p_offset = R*[eye_px_offset; eye_py_offset];
               
              rot_eye_px_offset(idx_ang) = rot_eye_p_offset(1,idx_ang);
              rot_eye_py_offset(idx_ang) = rot_eye_p_offset(2,idx_ang);
            end
              rot_eye_px_offset = rot_eye_px_offset.*amp_scaling;
              rot_eye_py_offset = rot_eye_py_offset.*amp_scaling;
              % get the reaction times
            temp_react = SACS_ALL_DATA.reaction(idx_tuned);
            
            % vmax of the tuned saccade
            eye_m_vmax_temp = SACS_ALL_DATA.eye_vm_max(idx_tuned);

            % compute velocity of the saccades
            vm_data   = loaded_data.trials_data.eye_vm_filt;
            time_data = loaded_data.trials_data.time_1K;

            sac_onset_times    = SACS_ALL_DATA.time_onset(idx_tuned);
            visual_onset_times = SACS_ALL_DATA.time_visual(idx_tuned);  
            trials_nums        = SACS_ALL_DATA.trial_num(idx_tuned);
          
            eye_vm_onset_temp    = nan(length_trace,length(sac_onset_times));
            for ii = 1:length(trials_nums)
                id_trial = trials_nums(ii);
                [~,idx1] = min(abs(time_data{1,id_trial} - sac_onset_times(ii)));
                eye_vm_onset_temp(:,ii) = vm_data{1,id_trial}(idx1-length_trace/2:idx1+length_trace/2-1); 
            end
            
            eye_vm_visual_temp   = nan(length_trace*2,length(visual_onset_times));
            for ii = 1:length(trials_nums)
                id_trial = trials_nums(ii);
                [~,idx2] = min(abs(time_data{1,id_trial} - visual_onset_times(ii)));
                if idx2  <= length_trace || idx2 >= length(vm_data{1,id_trial})-length_trace
                    continue
                else
                      eye_vm_visual_temp(:,ii) = vm_data{1,id_trial}(idx2-length_trace:idx2+length_trace-1);
                end
            end

             %  Validation for rotation (optional)
                if debug_fig 
                    fig = figure;
                    clf;
                    subplot(1,2,1);
                    plot(eye_px_offset, eye_py_offset,'Marker','.','LineStyle','none');
                    xlim([-10,10]);
                    ylim([-10,10]);
                    title('Before');
                    subplot(1,2,2);
                    plot(rot_eye_px_offset, rot_eye_py_offset,'Marker','.','LineStyle','none');
                    xlim([-10,10]);
                    ylim([-10,10]);
                    title('After'); 
                end
                % Concatenate data
                vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), vigor_temp];
                end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), offset_error];
                px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), rot_eye_px_offset];
                py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), rot_eye_py_offset]; 
                reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), temp_react];
                eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), mean(eye_vm_onset_temp,2,'omitnan')]; 
                eye_m_vel_visual.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [eye_m_vel_visual.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), mean(eye_vm_visual_temp,2,'omitnan')];
                eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), mean(eye_m_vmax_temp,'omitnan')];  
                
            end % tgt loop
         end % rec loop
     end % day loop
end % month loop


vigor_edges   = 0.5:0.1:1.5;
num_vigor_bin = length(vigor_edges) - 1;
vigor_x_axis = vigor_edges(1:end-1) + diff(vigor_edges)./2;      
day_fieldnames = fieldnames(vigor);
num_day        = length(day_fieldnames);
% Init .var
num_sac              = struct;
mean_end_error       = struct;
gen_var              = struct;
var_x                = struct;
end_error_y_axis     = struct;
end_error_sem_y_axis = struct;
gen_var_y_axis       = struct;
gen_var_sem_y_axis   = struct;
var_x_y_axis         = struct;
var_x_sem_y_axis     = struct;

for tag_idx = 1
    for counter_tgt = 1:length(tgt_cond_list)
        tgt_idx = tgt_cond_list(counter_tgt);
         for counter_vigor = 1 : num_vigor_bin
             num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))        = zeros(num_day,1);
             mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor)) = nan(num_day,1);
             gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))        = nan(num_day,1);
             var_x.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))          = nan(num_day,1);
         end
    end
end

% Loop thru sessions/days
for counter_day = 1 : num_day
    for tag_idx = 1
      for counter_tgt = 1:length(tgt_cond_list)
          tgt_idx = tgt_cond_list(counter_tgt);
            vigor_bin = discretize(vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), vigor_edges);
    
             for counter_vigor = 1 : num_vigor_bin
                 idx_vigor = vigor_bin == counter_vigor;
                 end_error_temp = end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_vigor);
                 px_offset_temp = px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_vigor);
                 py_offset_temp = py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_vigor);
    
                 num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1)        = length(end_error_temp);
                 mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1) = mean(end_error_temp);
                 gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1)        = det(cov(px_offset_temp, py_offset_temp));
                 var_x.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1)          = var(px_offset_temp);
             end
       end
    end
end

 % Loop thru vigor
for tag_idx = 1
    for counter_tgt = 1:length(tgt_cond_list)
        tgt_idx = tgt_cond_list(counter_tgt);
        for counter_vigor = 1 : num_vigor_bin
            num_sac_temp = num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor));
            num_sac_temp_weight = num_sac_temp ./ sum(num_sac_temp);

            mean_end_error_temp = mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor));
    
            end_error_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1, counter_vigor) = sum(mean_end_error_temp.*num_sac_temp_weight, 'omitnan');
            end_error_sem_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1,counter_vigor) = sqrt(var(mean_end_error_temp,num_sac_temp_weight,'omitnan'))./sqrt(sum(num_sac_temp>0));

            gen_var_temp   = gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor));
            var_x_temp     = var_x.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor));
            gen_var_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1,counter_vigor) = sum(gen_var_temp.*num_sac_temp_weight, 'omitnan');
            gen_var_sem_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1,counter_vigor) = sqrt(var(gen_var_temp,num_sac_temp_weight,'omitnan'))./sqrt(sum(num_sac_temp>0));

            var_x_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1,counter_vigor) = sum(var_x_temp.*num_sac_temp_weight, 'omitnan');
            var_x_sem_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1,counter_vigor) = sqrt(var(var_x_temp,num_sac_temp_weight,'omitnan'))./sqrt(sum(num_sac_temp>0));
        end
   end
end


% vigor and reaction times for High and Low rewards
vigor_prim_sac          = struct;
reaction_prim_sac       = struct;
vmax_prim_sac           = struct;
velocity_onset_prim_sac = struct;
velocity_visual_prim_sac = struct;

for counter_tgt = 1:length(tgt_cond_list)
   tgt_idx = tgt_cond_list(counter_tgt);
    for idx_d = 1:num_day
      vigor_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) =  nan(num_day,1); 
      reaction_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) =  nan(num_day,1);
      vmax_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) =  nan(num_day,1);
      velocity_onset_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) =  nan(length_trace,num_day);
      velocity_visual_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) =  nan(2*length_trace,num_day);
    end
end

for counter_tgt = 1:length(tgt_cond_list)
   tgt_idx = tgt_cond_list(counter_tgt);
    for idx_d = 1:num_day
      vigor_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_d,1)           =  median(vigor.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)),'omitnan'); 
      reaction_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_d,1)        =  median(reaction_times.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)),'omitnan');
      vmax_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_d,1)            =   mean(eye_m_vmax.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)),'omitnan');
      velocity_onset_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(:,idx_d)  =   mean(eye_vm_vel_onset.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)),2,'omitnan');
      velocity_visual_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(:,idx_d) =   mean(eye_m_vel_visual.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)),2,'omitnan');
    end
end


colors_ = [0.6350 0.0780 0.1840;0 0.4470 0.7410];
labels = {'High reward';'Low reward'};
fig = figure;
vigor_prim = [vigor_prim_sac.tag_1.tgt_1 vigor_prim_sac.tag_1.tgt_0 ];
violinplot(vigor_prim, labels, 'ViolinColor',colors_);
box off
ylabel('Saccade vigor')
ylim([0.7 1.2])
title('primary saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_1_vigor_median';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% Distribution of vigor
edge_vigor = 0:0.01:2;
h = histogram(vigor_prim_sac.tag_1.tgt_1,edge_vigor);
h.FaceColor = [0.5020 0.0941 0.1686];
hold on
h = histogram(vigor_prim_sac.tag_1.tgt_0,edge_vigor);
h.FaceColor = [0.2157 0.6039 0.8588];
ylabel('Frequency')
xlabel('Vigor')
box off
xlim([0.7 1.2])
title('primary saccade')
file_name = 'tag_1_vigor_pdf';
saveas(h,fullfile(file_path,file_name));
close gcf
fig = openfig(fullfile(file_path,file_name));
ESN_Beautify_Plot(fig, [8, 8], 11);
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

 % average of reaction times
RT = histogram(reaction_prim_sac.tag_1.tgt_1,local_RT_edges);
RT.FaceColor = [0.5020 0.0941 0.1686];
hold on
RT = histogram(reaction_prim_sac.tag_1.tgt_0,local_RT_edges);
RT.FaceColor = [0.2157 0.6039 0.8588];
ylabel('Frequency')
xlabel('Reaction time (ms)')
title('primary saccade')
xlim([ 50 300])
file_name = 'tag_1_reaction_time_median';
saveas(RT,fullfile(file_path,file_name));
close gcf
fig = openfig(fullfile(file_path,file_name));
ESN_Beautify_Plot(fig, [8, 8], 11);
saveas(fig,fullfile(file_path,file_name),'pdf');
close all


react_prim = [reaction_prim_sac.tag_1.tgt_1 reaction_prim_sac.tag_1.tgt_0 ];
fig = figure;
violinplot(react_prim, labels, 'ViolinColor',colors_);
ylabel('Reaction Times (ms)')
title('primary saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_1_reaction_time_pdf';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% MAX velocity

mn_Vmax  = [vmax_prim_sac.tag_1.tgt_1 vmax_prim_sac.tag_1.tgt_0];
fig = figure;
violinplot(mn_Vmax, labels, 'ViolinColor',colors_);
ylabel('Peak velocity (deg/s)')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_1_peak_velocity';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% average velocity saccade onset
speed_high_onset  = mean(velocity_onset_prim_sac.tag_1.tgt_1,2,'omitnan');
sem_sp_high_onset = std(velocity_onset_prim_sac.tag_1.tgt_1,0,2,'omitnan')/sqrt(size(velocity_onset_prim_sac.tag_1.tgt_1,2));

speed_low_onset  = mean(velocity_onset_prim_sac.tag_1.tgt_0,2,'omitnan');
sem_sp_low_onset = std(velocity_onset_prim_sac.tag_1.tgt_0,0,2,'omitnan')/sqrt(size(velocity_onset_prim_sac.tag_1.tgt_0,2));

fig=figure;
vel_x_axis = -length_trace/2:length_trace/2-1;
boundedline(vel_x_axis, speed_high_onset, sem_sp_high_onset,'color',[0.5020 0.0941 0.1686],'LineWidth',2)
hold on
boundedline(vel_x_axis, speed_low_onset, sem_sp_low_onset,'color',[0.2157 0.6039 0.8588],'LineWidth',2)
xlim([-50 100])
xlabel('Time from saccade onset (ms)')
ylabel('Velocity (deg/sec)')
title('primary saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_1_velocity_onset';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% average of eye velocity visual onset
speed_high_visual  = mean(velocity_visual_prim_sac.tag_1.tgt_1,2,'omitnan');
sem_sp_high_visual = std(velocity_visual_prim_sac.tag_1.tgt_1,0,2,'omitnan')/sqrt(size(velocity_visual_prim_sac.tag_1.tgt_1,2));

speed_low_visual  = mean(velocity_visual_prim_sac.tag_1.tgt_0,2,'omitnan');
sem_sp_low_visual = std(velocity_visual_prim_sac.tag_1.tgt_0,0,2,'omitnan')/sqrt(size(velocity_visual_prim_sac.tag_1.tgt_0,2));

fig=figure;
vel_x_axis = -length_trace:length_trace-1;
boundedline(vel_x_axis, speed_high_visual, sem_sp_high_visual,'color',[0.5020 0.0941 0.1686],'LineWidth',2)
hold on
boundedline(vel_x_axis, speed_low_visual, sem_sp_low_visual,'color',[0.2157 0.6039 0.8588],'LineWidth',2)
xlabel('Time from visual onset (ms)')
ylabel('Velocity (deg/sec)')
xlim([-100 300])
title('primary saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_1_velocity_visual';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% 1- High reward 
fig=figure;
Avg_high = end_error_y_axis.tag_1.tgt_1;
sem_high= end_error_sem_y_axis.tag_1.tgt_1;
boundedline(vigor_x_axis, Avg_high, sem_high,'color',[0.5020 0.0941 0.1686],'LineWidth',2)
hold on

% 2-Low reward
Avg_low= end_error_y_axis.tag_1.tgt_0;
Avg_low(Avg_low == 0) = nan;
sem_low = end_error_sem_y_axis.tag_1.tgt_0;

boundedline(vigor_x_axis, Avg_low, sem_low,'color',[0.2157 0.6039 0.8588],'LineWidth',2)
ylim([0.5 2])
xlabel('Saccade vigor');
ylabel('End point error (deg)');
title('primary saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);

file_name = 'tag_1_endpoint_error';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% Generalized variance
fig=figure;
Avg_high = gen_var_y_axis.tag_1.tgt_1;
sem_high= gen_var_sem_y_axis.tag_1.tgt_1;
boundedline(vigor_x_axis, Avg_high, sem_high,'color',[0.5020 0.0941 0.1686],'LineWidth',2)
hold on


Avg_low= gen_var_y_axis.tag_1.tgt_0;
Avg_low(Avg_low == 0) = nan;
sem_low = gen_var_sem_y_axis.tag_1.tgt_0;

boundedline(vigor_x_axis, Avg_low, sem_low,'color',[0.2157 0.6039 0.8588],'LineWidth',2)

xlabel('Saccade vigor');
ylabel('Generalized variance (deg^4)')
title('primary saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_1_gen_variance';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% Variance along x axis(saccade direction)

fig=figure;
Avg_high = var_x_y_axis.tag_1.tgt_1;
sem_high = var_x_sem_y_axis.tag_1.tgt_1;
boundedline(vigor_x_axis, Avg_high, sem_high,'color',[0.5020 0.0941 0.1686],'LineWidth',2)
hold on

Avg_low= var_x_y_axis.tag_1.tgt_0;
Avg_low(Avg_low == 0) = nan;
sem_low = var_x_sem_y_axis.tag_1.tgt_0;

boundedline(vigor_x_axis, Avg_low, sem_low,'color',[0.2157 0.6039 0.8588],'LineWidth',2)

xlabel('Saccade vigor');
ylabel('variance along saccade direction (deg^2)')
title('primary saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_1_variance_x_axis';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all
%% ACCURACY VS. VIGOR PRIMARY WITH RPE
% 1. End point error magnitude
% 2. Generalized variance (determinant of covariance matrix)
tag_idx   = 1; % only applicable for primary saccade
task_id   = 1; % for forced trials
debug_fig = false;
% Load vigor model
loaded_data = load('C:\Users\Jafar\Documents\reward\data_population\vigor_model.mat');
fit_obj     = loaded_data.fit_obj.(sprintf('m_%s', '132F')).(sprintf('tag_%d', tag_fit_idx));

% Init. data
vigor            = struct;
end_error        = struct;
px_offset        = struct;
py_offset        = struct;
reaction_times   = struct;
eye_vm_vel_onset = struct;
eye_m_vel_visual = struct;
eye_m_vmax       = struct;

% Loop animals
month_dir = dir2(fullfile(data_path,sprintf('data_%s','132F')));
       
counter_day = 0; % consider each day as a session
  % Loop months
 num_month = size(month_dir,1);
 for m_idx =  [17 18 19] %1:num_month
     month_path = fullfile(month_dir(m_idx).folder,month_dir(m_idx).name);
     day_dir = dir2(month_path);
     day_dir = day_dir([day_dir.isdir]);
    % Loop days
     num_day = size(day_dir,1);
    for d_idx = 1:num_day
        counter_day = counter_day + 1;
        rec_path = fullfile(day_dir(d_idx).folder, day_dir(d_idx).name);
        rec_dir = dir2(rec_path);
        rec_dir = rec_dir([rec_dir.isdir]);
        % Loop recs
        num_rec = size(rec_dir,1);
        for r_idx =1:num_rec
            fprintf('loading session: %s.........\n', rec_dir(r_idx).name);
    
            data_dir = dir2(fullfile(rec_dir(r_idx).folder,rec_dir(r_idx).name, 'analyzed_data\behavior_data\eye','*ANALYZED*'));
            % Load data; if not available, skip the rec
            if isempty(data_dir)
               continue
            end
            loaded_data = load(fullfile(data_dir.folder, data_dir.name));
            SACS_ALL_DATA = loaded_data.sac_data;
            SACS_amp_bin  = discretize(SACS_ALL_DATA.eye_amp_m,  local_amp_edges);
            
            for counter_tgt = 1:length(tgt_cond_list)
                tgt_idx = tgt_cond_list(counter_tgt);
                 fprintf('Target %d...\n',  tgt_idx);
                for counter_rew = 1:length(rew_cond_list)
                    rew_idx = rew_cond_list(counter_rew);
                    fprintf('Reward %d...\n',  rew_idx);
            
             % Init. var
             vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
             end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
             px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
             py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
             reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
             eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
             eye_m_vel_visual.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
             eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];

             idx_validity  = SACS_ALL_DATA.validity;
             idx_reaction  = (SACS_ALL_DATA.reaction > min_rxn_time) & (SACS_ALL_DATA.reaction < max_rxn_time);
             idx_tag       = SACS_ALL_DATA.tag == tag_idx;
             idx_task_cond = SACS_ALL_DATA.task_cond == task_id;
             idx_tgt       = SACS_ALL_DATA.tgt_cond == tgt_cond_list(counter_tgt);
             idx_rew       = SACS_ALL_DATA.rew_cond == rew_cond_list(counter_rew);

             if ismember(tag_idx,[1,6])
                idx_amp  = SACS_ALL_DATA.eye_amp_m > tag_1_amp_threshold;
             elseif tag_idx == 4
                 idx_amp = SACS_ALL_DATA.eye_amp_m > tag_4_amp_threshold;
             else
                 idx_amp = SACS_ALL_DATA.eye_amp_m;
             end
            idx_tuned     = idx_tag & idx_validity & idx_reaction & idx_task_cond & idx_amp & idx_tgt & idx_rew;
            % Get relevant saccade variables; easier to debug this way
            eye_amp_m = SACS_ALL_DATA.eye_amp_m(idx_tuned);
            eye_vm_max = SACS_ALL_DATA.eye_vm_max(idx_tuned);
            visual_px_offset = SACS_ALL_DATA.cue_x(idx_tuned) - SACS_ALL_DATA.start_x(idx_tuned); % centering saccade
            visual_py_offset = SACS_ALL_DATA.cue_y(idx_tuned) - SACS_ALL_DATA.start_y(idx_tuned); % centering saccade
            eye_px_onset   = 0; % centering saccade
            eye_py_onset   = 0; % centering saccade
            eye_px_offset = SACS_ALL_DATA.eye_px_offset(idx_tuned) - SACS_ALL_DATA.start_x(idx_tuned); % centering saccade
            eye_py_offset = SACS_ALL_DATA.eye_py_offset(idx_tuned) - SACS_ALL_DATA.start_y(idx_tuned); % centering saccade
            % Compute vigor
            mean_vel   = predict(fit_obj,log10(eye_amp_m'));
            vigor_temp = log10(eye_vm_max')./mean_vel;
            % Compute saccade offset error
            offset_error = sqrt((SACS_ALL_DATA.cue_x(idx_tuned) - SACS_ALL_DATA.eye_px_offset(idx_tuned)).^2 +...
            (SACS_ALL_DATA.cue_y(idx_tuned) - SACS_ALL_DATA.eye_py_offset(idx_tuned)).^2);
        
            % Compute nominal direction of saccades (deg)
            delta_x = visual_px_offset - eye_px_onset;
            delta_y = visual_py_offset - eye_py_onset;
            visual_ang      = wrapTo360(atan2d(delta_y, delta_x));
            visual_ang_list = unique(visual_ang);
            
            % Compute nominal amp. of saccade & scaling factor to
            % normalize amp.
            nominal_amp = sqrt(delta_x.^2 + delta_y.^2);
%                     amp_scaling = norm_amp./nominal_amp;
            amp_scaling = 1;
            % Compute rotational matrix to rotate saccade vector to
            % right
            rot_eye_px_offset = zeros(1, length(visual_ang));
            rot_eye_py_offset = zeros(1, length(visual_ang));

            for visual_ang_temp = visual_ang_list
                idx_ang = visual_ang == visual_ang_temp;

                R = [cosd(-visual_ang_temp), -sind(-visual_ang_temp);
                     sind(-visual_ang_temp), cosd(-visual_ang_temp)];
              rot_eye_p_offset = R*[eye_px_offset; eye_py_offset];
               
              rot_eye_px_offset(idx_ang) = rot_eye_p_offset(1,idx_ang);
              rot_eye_py_offset(idx_ang) = rot_eye_p_offset(2,idx_ang);
            end
              rot_eye_px_offset = rot_eye_px_offset.*amp_scaling;
              rot_eye_py_offset = rot_eye_py_offset.*amp_scaling;
              % get the reaction times
            temp_react = SACS_ALL_DATA.reaction(idx_tuned);
            
            % vmax of the tuned saccade
            eye_m_vmax_temp = SACS_ALL_DATA.eye_vm_max(idx_tuned);

            % compute velocity of the saccades
            vm_data   = loaded_data.trials_data.eye_vm_filt;
            time_data = loaded_data.trials_data.time_1K;

            sac_onset_times    = SACS_ALL_DATA.time_onset(idx_tuned);
            visual_onset_times = SACS_ALL_DATA.time_visual(idx_tuned);  
            trials_nums        = SACS_ALL_DATA.trial_num(idx_tuned);
          
            eye_vm_onset_temp    = nan(length_trace,length(sac_onset_times));
            for ii = 1:length(trials_nums)
                id_trial = trials_nums(ii);
                [~,idx1] = min(abs(time_data{1,id_trial} - sac_onset_times(ii)));
                eye_vm_onset_temp(:,ii) = vm_data{1,id_trial}(idx1-length_trace/2:idx1+length_trace/2-1); 
            end
            
            eye_vm_visual_temp   = nan(length_trace*2,length(visual_onset_times));
            for ii = 1:length(trials_nums)
                id_trial = trials_nums(ii);
                [~,idx2] = min(abs(time_data{1,id_trial} - visual_onset_times(ii)));
                if idx2  <= length_trace || idx2 >= length(vm_data{1,id_trial})-length_trace
                    continue
                else
                      eye_vm_visual_temp(:,ii) = vm_data{1,id_trial}(idx2-length_trace:idx2+length_trace-1);
                end
            end

             %  Validation for rotation (optional)
                if debug_fig 
                    fig = figure;
                    clf;
                    subplot(1,2,1);
                    plot(eye_px_offset, eye_py_offset,'Marker','.','LineStyle','none');
                    xlim([-10,10]);
                    ylim([-10,10]);
                    title('Before');
                    subplot(1,2,2);
                    plot(rot_eye_px_offset, rot_eye_py_offset,'Marker','.','LineStyle','none');
                    xlim([-10,10]);
                    ylim([-10,10]);
                    title('After'); 
                end
                % Concatenate data
                vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), vigor_temp];
                end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), offset_error];
                px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), rot_eye_px_offset];
                py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), rot_eye_py_offset]; 
                reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), temp_react];
                eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), mean(eye_vm_onset_temp,2,'omitnan')]; 
                eye_m_vel_visual.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [eye_m_vel_visual.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), mean(eye_vm_visual_temp,2,'omitnan')];
                eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), mean(eye_m_vmax_temp,'omitnan')];  
                end %  rew loop  
            end % tgt loop
         end % rec loop
     end % day loops
end % month loop



% vigor_edges   = 0.5:0.2:3.0;
vigor_edges   = 0.9:0.05:1.1;
num_vigor_bin = length(vigor_edges) - 1;
vigor_x_axis = vigor_edges(1:end-1) + diff(vigor_edges)./2;      
day_fieldnames = fieldnames(vigor);
num_day        = length(day_fieldnames);
% Init .var
end_error_y_axis     = struct;
end_error_sem_y_axis = struct;
gen_var_y_axis       = struct;
gen_var_sem_y_axis   = struct;



for counter_tgt = 1:length(tgt_cond_list)
    tgt_idx = tgt_cond_list(counter_tgt);
    for counter_rew = 1:length(rew_cond_list)
        rew_idx = rew_cond_list(counter_rew);
        for counter_vigor = 1 : num_vigor_bin
            num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor))        = zeros(num_day,1);
            mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor)) = nan(num_day,1);
            gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor))        = nan(num_day,1);
        end
    end
end

% Loop thru sessions/days
for counter_day = 1 : num_day
   for counter_tgt = 1:length(tgt_cond_list)
     tgt_idx = tgt_cond_list(counter_tgt);
      for counter_rew = 1:length(rew_cond_list)
          rew_idx = rew_cond_list(counter_rew);
            vigor_bin = discretize(vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), vigor_edges);

         for counter_vigor = 1 : num_vigor_bin
             idx_vigor = vigor_bin == counter_vigor;
             end_error_temp = end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_vigor);
             px_offset_temp = px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_vigor);
             py_offset_temp = py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_vigor);

             num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1)        = length(end_error_temp);
             mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1) = mean(end_error_temp);
             gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1)        = det(cov(px_offset_temp, py_offset_temp));
         end
      end
   end
end

 % Loop thru vigor
for counter_tgt = 1:length(tgt_cond_list)
    tgt_idx = tgt_cond_list(counter_tgt);
    for counter_rew = 1:length(rew_cond_list)
        rew_idx = rew_cond_list(counter_rew);
        for counter_vigor = 1 : num_vigor_bin
            num_sac_temp = num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor));
            num_sac_temp_weight = num_sac_temp ./ sum(num_sac_temp);

            mean_end_error_temp = mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor));
    
            end_error_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(1, counter_vigor) = sum(mean_end_error_temp.*num_sac_temp_weight, 'omitnan');
            end_error_sem_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(1,counter_vigor) = sqrt(var(mean_end_error_temp,num_sac_temp_weight,'omitnan'))./sqrt(sum(num_sac_temp>0));

            gen_var_temp = gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor));
            gen_var_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(1,counter_vigor) = sum(gen_var_temp.*num_sac_temp_weight, 'omitnan');
            gen_var_sem_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(1,counter_vigor) = sqrt(var(gen_var_temp,num_sac_temp_weight,'omitnan'))./sqrt(sum(num_sac_temp>0));
        end  
   end
end

% vigor and reaction times for High and Low rewards
vigor_prim_sac          = struct;
reaction_prim_sac       = struct;
vmax_prim_sac           = struct;
velocity_onset_prim_sac = struct;

for counter_tgt = 1:length(tgt_cond_list)
    tgt_idx = tgt_cond_list(counter_tgt);
    for counter_rew = 1:length(rew_cond_list)
        rew_idx = rew_cond_list(counter_rew);
        for idx_d = 1:num_day
              vigor_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =  nan(num_day,1); 
              reaction_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =  nan(num_day,1);
              vmax_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =  nan(num_day,1);
              velocity_onset_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =  nan(length_trace,num_day);
        end
    end
end

for counter_tgt = 1:length(tgt_cond_list)
    tgt_idx = tgt_cond_list(counter_tgt);
    for counter_rew = 1:length(rew_cond_list)
        rew_idx = rew_cond_list(counter_rew);
        for idx_d = 1:num_day
              vigor_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_d,1) =...
                  mean(vigor.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)),'omitnan');
              reaction_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_d,1) =...
                  mean(reaction_times.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)),'omitnan');
              vmax_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_d,1) =...
                  mean(eye_m_vmax.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)),'omitnan');
              velocity_onset_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(:,idx_d) =...
                  eye_vm_vel_onset.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx));
        end
    end
end


% plot
labels = {'CTRL (HH)';'-RPE';'+RPE';'CTRL (LL)'};
vig_prim = [vigor_prim_sac.tag_1.tgt_1.rew_1 vigor_prim_sac.tag_1.tgt_1.rew_0 vigor_prim_sac.tag_1.tgt_0.rew_1 vigor_prim_sac.tag_1.tgt_0.rew_0];
colors_ = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.9290 0.6940 0.1250;0.4660 0.6740 0.1880];
% vigor
fig = figure;
violinplot(vig_prim, labels, 'ViolinColor',colors_);
ylabel('Vigor')
xtickangle(45)
box off
title('Primary saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_1_vigor';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all


% Reaction times
reaction_corr = [reaction_prim_sac.tag_1.tgt_1.rew_1 reaction_prim_sac.tag_1.tgt_1.rew_0 reaction_prim_sac.tag_1.tgt_0.rew_1 reaction_prim_sac.tag_1.tgt_0.rew_0];
fig = figure;
violinplot(reaction_corr, labels, 'ViolinColor',colors_);
ylabel('Reaction Times (ms)')
xtickangle(45)
title('Primary saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_1_reaction_time';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

vmax_corr = [vmax_prim_sac.tag_1.tgt_1.rew_1 vmax_prim_sac.tag_1.tgt_1.rew_0  vmax_prim_sac.tag_1.tgt_0.rew_1  vmax_prim_sac.tag_1.tgt_0.rew_0];
fig = figure;
violinplot(vmax_corr, labels, 'ViolinColor',colors_);
xtickangle(45)
ylabel('Peak Velocity (deg/sec)')
title('Primary saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_1_peak_velocity';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% velocity
vel_onset_HH    = mean(velocity_onset_prim_sac.tag_1.tgt_1.rew_1,2,'omitnan');
sem_onset_HH    = std(velocity_onset_prim_sac.tag_1.tgt_1.rew_1,0,2,'omitnan')/sqrt(size(velocity_onset_prim_sac.tag_1.tgt_1.rew_1,2));
vel_onset_HL    = mean(velocity_onset_prim_sac.tag_1.tgt_1.rew_0,2,'omitnan');
sem_onset_HL    = std(velocity_onset_prim_sac.tag_1.tgt_1.rew_0,0,2,'omitnan')/sqrt(size(velocity_onset_prim_sac.tag_1.tgt_1.rew_0,2));
vel_onset_LH    = mean(velocity_onset_prim_sac.tag_1.tgt_0.rew_1,2,'omitnan');
sem_onset_LH    = std(velocity_onset_prim_sac.tag_1.tgt_0.rew_1,0,2,'omitnan')/sqrt(size(velocity_onset_prim_sac.tag_1.tgt_0.rew_1,2));
vel_onset_LL    = mean(velocity_onset_prim_sac.tag_1.tgt_0.rew_0,2,'omitnan');
sem_onset_LL    = std(velocity_onset_prim_sac.tag_1.tgt_0.rew_0,0,2,'omitnan')/sqrt(size(velocity_onset_prim_sac.tag_1.tgt_0.rew_0,2));

 
fig = figure;
hold on
boundedline(inds_span, vel_onset_HH, sem_onset_HH,'color',colors_(1,:),'LineWidth',2)
boundedline(inds_span, vel_onset_HL, sem_onset_HL,'color',colors_(2,:),'LineWidth',2)
boundedline(inds_span, vel_onset_LH, sem_onset_LH,'color',colors_(3,:),'LineWidth',2)
boundedline(inds_span, vel_onset_LL, sem_onset_LL,'color',colors_(4,:),'LineWidth',2)
xlabel('Time from saccade onset (ms)')
ylabel('Velocity (deg/sec)')
title('Primary saccade')
xlim([-50 100])
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_1_velocity';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all 

 % end point error
end_error_prim = [end_error_y_axis.tag_1.tgt_1.rew_1; end_error_y_axis.tag_1.tgt_1.rew_0;  end_error_y_axis.tag_1.tgt_0.rew_1;  end_error_y_axis.tag_1.tgt_0.rew_0]; 
end_error_prim(end_error_prim == 0)=nan;
fig = figure;
hold on
boundedline(vigor_x_axis, end_error_prim(1,:),end_error_sem_y_axis.tag_1.tgt_1.rew_1 ,'color',colors_(1,:),'LineWidth',2)
boundedline(vigor_x_axis, end_error_prim(2,:),end_error_sem_y_axis.tag_1.tgt_1.rew_0 ,'color',colors_(2,:),'LineWidth',2)
boundedline(vigor_x_axis, end_error_prim(3,:),end_error_sem_y_axis.tag_1.tgt_0.rew_1 ,'color',colors_(3,:),'LineWidth',2)
boundedline(vigor_x_axis, end_error_prim(4,:),end_error_sem_y_axis.tag_1.tgt_0.rew_0 ,'color',colors_(4,:),'LineWidth',2)
xlabel('Saccade vigor');
ylabel('End point error (deg)');
ESN_Beautify_Plot(fig, [8, 8], 11);
title('Primary saccade')
file_name = 'tag_1_RPE_endpoint_error';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all   

% plot generalized variance  
fig = figure;
hold on
boundedline(vigor_x_axis, gen_var_y_axis.tag_1.tgt_1.rew_1,gen_var_sem_y_axis.tag_1.tgt_1.rew_1 ,'color',colors_(1,:),'LineWidth',2)
boundedline(vigor_x_axis, gen_var_y_axis.tag_1.tgt_1.rew_0,gen_var_sem_y_axis.tag_1.tgt_1.rew_0 ,'color',colors_(2,:),'LineWidth',2)
boundedline(vigor_x_axis, gen_var_y_axis.tag_1.tgt_0.rew_1,gen_var_sem_y_axis.tag_1.tgt_0.rew_1 ,'color',colors_(3,:),'LineWidth',2)
boundedline(vigor_x_axis, gen_var_y_axis.tag_1.tgt_0.rew_0,gen_var_sem_y_axis.tag_1.tgt_0.rew_0 ,'color',colors_(4,:),'LineWidth',2)
xlabel('Saccade vigor');
ylabel('Generalized variance (deg^4)')
ESN_Beautify_Plot(fig, [8, 8], 11);
title('primary saccade')
file_name = 'tag_1_RPE_gen_variance';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all  


%% ACCURACY VS. VIGOR CORRECTIVE SACCADE
% 1. End point error magnitude
% 2. Generalized variance (determinant of covariance matrix)
fit_idx  = 1; % 1 for tag 1,4,6
tag_idx  = 4;
task_id  = 1; % for forced trials

debug_fig = false; % whether to plot to make sure rotation was done correctly
% Load vigor model
loaded_data = load('C:\Users\Jafar\Documents\reward\data_population\vigor_model.mat');
fit_obj     = loaded_data.fit_obj.(sprintf('m_%s', '132F')).(sprintf('tag_%d', tag_fit_idx));

% Init. data
vigor            = struct;
end_error        = struct;
px_offset        = struct;
py_offset        = struct;
reaction_times   = struct;
eye_vm_vel_onset = struct;
eye_m_vmax       = struct;


% Loop animals
month_dir = dir2(fullfile(data_path,sprintf('data_%s','132F')));
    
 % Get vigor model per animal
% fit_params = coeffvalues(fit_obj.(sprintf('m_%s', '132F')).(sprintf('tag_%d', fit_idx)));
    
counter_day = 0; % consider each day as a session
% Loop months
 num_month = size(month_dir,1);
 for m_idx =  [17 18 19] %1:num_month
     month_path = fullfile(month_dir(m_idx).folder,month_dir(m_idx).name);
     day_dir = dir2(month_path);
     day_dir = day_dir([day_dir.isdir]);
    % Loop days
     num_day = size(day_dir,1);
    for d_idx = 1:num_day
        counter_day = counter_day + 1;
        rec_path = fullfile(day_dir(d_idx).folder, day_dir(d_idx).name);
        rec_dir = dir2(rec_path);
        rec_dir = rec_dir([rec_dir.isdir]);
        % Loop recs
        num_rec = size(rec_dir,1);
        for r_idx =1:num_rec
            fprintf('loading session: %s.........\n', rec_dir(r_idx).name);
    
            data_dir = dir2(fullfile(rec_dir(r_idx).folder,rec_dir(r_idx).name, 'analyzed_data\behavior_data\eye','*ANALYZED*'));
            % Load data; if not available, skip the rec
            if isempty(data_dir)
               continue
            end
            loaded_data = load(fullfile(data_dir.folder, data_dir.name));
            SACS_ALL_DATA = loaded_data.sac_data;
            SACS_amp_bin  = discretize(SACS_ALL_DATA.eye_amp_m,  local_amp_edges);

            for counter_tgt = 1:length(tgt_cond_list)
                tgt_idx = tgt_cond_list(counter_tgt);
                fprintf('Target %d...\n',  tgt_idx);
                for counter_rew = 1:length(rew_cond_list)
                    rew_idx = rew_cond_list(counter_rew);
                  fprintf('Reward condition %d...\n',  rew_idx)
                     % Init. var
                     vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
                     end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
                     px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
                     py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
                     reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
                     eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
                     eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = [];
    
                     idx_validity  = SACS_ALL_DATA.validity;
                     idx_reaction  = (SACS_ALL_DATA.reaction > min_rxn_time) & (SACS_ALL_DATA.reaction < max_rxn_time);
                     idx_tag       = SACS_ALL_DATA.tag == tag_idx;
                     idx_task_cond = SACS_ALL_DATA.task_cond == task_id;
                     idx_tgt       = SACS_ALL_DATA.tgt_cond == tgt_cond_list(counter_tgt);
                     idx_rew       = SACS_ALL_DATA.rew_cond == rew_cond_list(counter_rew);
              

                    if ismember(tag_idx,[1,6])
                        idx_amp  = SACS_ALL_DATA.eye_amp_m > tag_1_amp_threshold;
                    elseif tag_idx == 4
                        idx_amp = SACS_ALL_DATA.eye_amp_m > tag_4_amp_threshold;
                    else
                      idx_amp = SACS_ALL_DATA.eye_amp_m;
                    end

                    idx_tuned = idx_tag & idx_validity & idx_reaction & idx_task_cond & idx_tgt & idx_rew;
       
                    % Get relevant saccade variables; easier to debug this way
                    eye_amp_m = SACS_ALL_DATA.eye_amp_m(idx_tuned);
                    eye_vm_max = SACS_ALL_DATA.eye_vm_max(idx_tuned);
                    visual_px_offset = SACS_ALL_DATA.end_x(idx_tuned) - SACS_ALL_DATA.cue_x(idx_tuned); % centering saccade
                    visual_py_offset = SACS_ALL_DATA.end_y(idx_tuned) - SACS_ALL_DATA.cue_y(idx_tuned); % centering saccade
                    eye_px_onset   = 0; % centering saccade
                    eye_py_onset   = 0; % centering saccade
                    eye_px_offset = SACS_ALL_DATA.eye_px_offset(idx_tuned) - SACS_ALL_DATA.cue_x(idx_tuned); % centering saccade
                    eye_py_offset = SACS_ALL_DATA.eye_py_offset(idx_tuned) - SACS_ALL_DATA.cue_y(idx_tuned); % centering saccade
                    % Compute vigor
                   mean_vel   = predict(fit_obj,log10(eye_amp_m'));
                   vigor_temp = log10(eye_vm_max')./mean_vel;
                    % Compute saccade offset error
                    offset_error = sqrt((SACS_ALL_DATA.end_x(idx_tuned) - SACS_ALL_DATA.eye_px_offset(idx_tuned)).^2 +...
                        (SACS_ALL_DATA.end_y(idx_tuned) - SACS_ALL_DATA.eye_py_offset(idx_tuned)).^2);
                
%                   Compute nominal direction of saccades (deg)
                    delta_x = visual_px_offset - eye_px_onset;
                    delta_y = visual_py_offset - eye_py_onset;
                    visual_ang      = wrapTo360(atan2d(delta_y, delta_x));
                    visual_ang_list = unique(visual_ang);
                    % Compute nominal amp. of saccade & scaling factor to
                    % normalize amp.
                    nominal_amp = sqrt(delta_x.^2 + delta_y.^2);
%                   amp_scaling = norm_amp./nominal_amp;
                    amp_scaling = 1;
                    % Compute rotational matrix to rotate saccade vector to
                    % right
                    rot_eye_px_offset = zeros(1, length(visual_ang));
                    rot_eye_py_offset = zeros(1, length(visual_ang));

                    for visual_ang_temp = visual_ang_list
                        idx_ang = visual_ang == visual_ang_temp;
    
                        R = [cosd(-visual_ang_temp), -sind(-visual_ang_temp);
                             sind(-visual_ang_temp), cosd(-visual_ang_temp)];
                      rot_eye_p_offset = R*[eye_px_offset; eye_py_offset];
                       
                      rot_eye_px_offset(idx_ang) = rot_eye_p_offset(1,idx_ang);
                      rot_eye_py_offset(idx_ang) = rot_eye_p_offset(2,idx_ang);
                    end
                      rot_eye_px_offset = rot_eye_px_offset.*amp_scaling;
                      rot_eye_py_offset = rot_eye_py_offset.*amp_scaling;
                    
                    temp_react = SACS_ALL_DATA.reaction(idx_tuned);

%                     Validation for rotation (optional)
                        if debug_fig 
                            fig = figure;
                            clf;
                            subplot(1,2,1);
                            plot(eye_px_offset, eye_py_offset,'Marker','.','LineStyle','none');
                            xlim([-10,10]);
                            ylim([-10,10]);
                            title('Before');
                            subplot(1,2,2);
                            plot(rot_eye_px_offset, rot_eye_py_offset,'Marker','.','LineStyle','none');
                            xlim([-10,10]);
                            ylim([-10,10]);
                            title('After');
                            
                        end
            
                       % compute velocity of the saccades
                        vm_data   = loaded_data.trials_data.eye_vm_filt;
                        time_data = loaded_data.trials_data.time_1K;
            
                        sac_onset_times    = SACS_ALL_DATA.time_onset(idx_tuned);
                        trials_nums        = SACS_ALL_DATA.trial_num(idx_tuned);
                      
                        eye_vm_onset_temp    = nan(2*win_size,length(sac_onset_times));
                        for ii = 1:length(trials_nums)
                            id_trial = trials_nums(ii);
                            [~,idx1] = min(abs(time_data{1,id_trial} - sac_onset_times(ii)));
                            if idx1 + win_size > length(time_data{1,id_trial})
                                continue
                            else
                                 eye_vm_onset_temp(:,ii) = vm_data{1,id_trial}(idx1-win_size:idx1+win_size-1); 
                            end
                        end

                        % Concatenate data
                    vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =...
                        [vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), vigor_temp];

                    end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =...
                        [end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), offset_error];

                    px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =...
                        [px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), rot_eye_px_offset];

                    py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =...
                        [py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), rot_eye_py_offset];
                    
                    reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) = ...
                        [reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), temp_react];
                     
                    eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =...
                        [eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), mean(eye_vm_onset_temp,2,'omitnan')]; 

                   eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =...
                       [eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), mean(eye_vm_max,'omitnan')];  
               end % rew_amt loop      
            end % tgt loop
         end % rec loop
     end % day loop
end % month loop


vigor_edges   = 0.9:0.05:1.1;
num_vigor_bin = length(vigor_edges) - 1;
vigor_x_axis = vigor_edges(1:end-1) + diff(vigor_edges)./2;      
day_fieldnames = fieldnames(vigor);
num_day        = length(day_fieldnames);
% Init .var
end_error_y_axis     = struct;
end_error_sem_y_axis = struct;
gen_var_y_axis       = struct;
gen_var_sem_y_axis   = struct;


for counter_tgt = 1:length(tgt_cond_list)
    tgt_idx = tgt_cond_list(counter_tgt);
    for counter_rew = 1:length(rew_cond_list)
        rew_idx = rew_cond_list(counter_rew);
        for counter_vigor = 1 : num_vigor_bin
            num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor))        = zeros(num_day,1);
            mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor)) = nan(num_day,1);
            gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor))        = nan(num_day,1);
        end
    end
end

% Loop thru sessions/days
for counter_day = 1 : num_day
   for counter_tgt = 1:length(tgt_cond_list)
     tgt_idx = tgt_cond_list(counter_tgt);
      for counter_rew = 1:length(rew_cond_list)
          rew_idx = rew_cond_list(counter_rew);
            vigor_bin = discretize(vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)), vigor_edges);

         for counter_vigor = 1 : num_vigor_bin
             idx_vigor = vigor_bin == counter_vigor;
             end_error_temp = end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_vigor);
             px_offset_temp = px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_vigor);
             py_offset_temp = py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_vigor);

             num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1)        = length(end_error_temp);
             mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1) = mean(end_error_temp);
             gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1)        = det(cov(px_offset_temp, py_offset_temp));
         end
      end
   end
end

 % Loop thru vigor
for counter_tgt = 1:length(tgt_cond_list)
    tgt_idx = tgt_cond_list(counter_tgt);
    for counter_rew = 1:length(rew_cond_list)
        rew_idx = rew_cond_list(counter_rew);
        for counter_vigor = 1 : num_vigor_bin
            num_sac_temp = num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor));
            num_sac_temp_weight = num_sac_temp ./ sum(num_sac_temp);

            mean_end_error_temp = mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor));
    
            end_error_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(1, counter_vigor) = sum(mean_end_error_temp.*num_sac_temp_weight, 'omitnan');
            end_error_sem_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(1,counter_vigor) = sqrt(var(mean_end_error_temp,num_sac_temp_weight,'omitnan'))./sqrt(sum(num_sac_temp>0));

            gen_var_temp = gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)).(sprintf('bin_%d',counter_vigor));
            gen_var_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(1,counter_vigor) = sum(gen_var_temp.*num_sac_temp_weight, 'omitnan');
            gen_var_sem_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(1,counter_vigor) = sqrt(var(gen_var_temp,num_sac_temp_weight,'omitnan'))./sqrt(sum(num_sac_temp>0));
        end  
   end
end

% vigor and reaction times for High and Low rewards
vigor_corr_sac          = struct;
reaction_corr_sac       = struct;
vmax_corr_sac           = struct;
velocity_onset_corr_sac = struct;

for counter_tgt = 1:length(tgt_cond_list)
    tgt_idx = tgt_cond_list(counter_tgt);
    for counter_rew = 1:length(rew_cond_list)
        rew_idx = rew_cond_list(counter_rew);
        for idx_d = 1:num_day
              vigor_corr_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =  nan(num_day,1); 
              reaction_corr_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =  nan(num_day,1);
              vmax_corr_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =  nan(num_day,1);
              velocity_onset_corr_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)) =  nan(2*win_size,num_day);
        end
    end
end

for counter_tgt = 1:length(tgt_cond_list)
    tgt_idx = tgt_cond_list(counter_tgt);
    for counter_rew = 1:length(rew_cond_list)
        rew_idx = rew_cond_list(counter_rew);
        for idx_d = 1:num_day
              vigor_corr_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_d,1) =...
                  mean(vigor.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)),'omitnan');
              reaction_corr_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_d,1) =...
                  mean(reaction_times.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)),'omitnan');
              vmax_corr_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(idx_d,1) =...
                  mean(eye_m_vmax.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx)),'omitnan');
              velocity_onset_corr_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx))(:,idx_d) =...
                  eye_vm_vel_onset.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('rew_%d',rew_idx));
        end
    end
end


labels = {'CTRL (HH)';'-RPE';'+RPE';'CTRL (LL)'};
vig_corr = [vigor_corr_sac.tag_4.tgt_1.rew_1 vigor_corr_sac.tag_4.tgt_1.rew_0 vigor_corr_sac.tag_4.tgt_0.rew_1 vigor_corr_sac.tag_4.tgt_0.rew_0];
colors_ = [0.6350 0.0780 0.1840;0 0.4470 0.7410;0.9290 0.6940 0.1250;0.4660 0.6740 0.1880];

% vigor
fig = figure;
violinplot(vig_corr, labels, 'ViolinColor',colors_);
ylabel('Vigor')
xtickangle(45)
box off
title('corrective saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_4_vigor';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% Reaction times
reaction_corr = [reaction_corr_sac.tag_4.tgt_1.rew_1 reaction_corr_sac.tag_4.tgt_1.rew_0 reaction_corr_sac.tag_4.tgt_0.rew_1 reaction_corr_sac.tag_4.tgt_0.rew_0];
fig = figure;
violinplot(reaction_corr, labels, 'ViolinColor',colors_);
xtickangle(45)
ylabel('Reaction Times (ms)')
title('corrective saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_4_reaction_time';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

vmax_corr = [vmax_corr_sac.tag_4.tgt_1.rew_1 vmax_corr_sac.tag_4.tgt_1.rew_0  vmax_corr_sac.tag_4.tgt_0.rew_1  vmax_corr_sac.tag_4.tgt_0.rew_0];
fig = figure;
violinplot(vmax_corr, labels, 'ViolinColor',colors_);
xtickangle(45)
ylabel('Peak Velocity (deg/sec)')
title('corrective saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_4_peak_velocity';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% velocity
vel_onset_HH    = mean(velocity_onset_corr_sac.tag_4.tgt_1.rew_1,2,'omitnan');
sem_onset_HH    = std(velocity_onset_corr_sac.tag_4.tgt_1.rew_1,0,2,'omitnan')/sqrt(size(velocity_onset_corr_sac.tag_4.tgt_1.rew_1,2));
vel_onset_HL    = mean(velocity_onset_corr_sac.tag_4.tgt_1.rew_0,2,'omitnan');
sem_onset_HL    = std(velocity_onset_corr_sac.tag_4.tgt_1.rew_0,0,2,'omitnan')/sqrt(size(velocity_onset_corr_sac.tag_4.tgt_1.rew_0,2));
vel_onset_LH    = mean(velocity_onset_corr_sac.tag_4.tgt_0.rew_1,2,'omitnan');
sem_onset_LH    = std(velocity_onset_corr_sac.tag_4.tgt_0.rew_1,0,2,'omitnan')/sqrt(size(velocity_onset_corr_sac.tag_4.tgt_0.rew_1,2));
vel_onset_LL    = mean(velocity_onset_corr_sac.tag_4.tgt_0.rew_0,2,'omitnan');
sem_onset_LL    = std(velocity_onset_corr_sac.tag_4.tgt_0.rew_0,0,2,'omitnan')/sqrt(size(velocity_onset_corr_sac.tag_4.tgt_0.rew_0,2));

time_axis = -win_size: 1:win_size-1;   
fig = figure;
hold on
[h1,h2] = boundedline(time_axis, vel_onset_HH, sem_onset_HH,'color',colors_(1,:),'LineWidth',2);
[p1,p2] = boundedline(time_axis, vel_onset_HL, sem_onset_HL,'color',colors_(2,:),'LineWidth',2);
[g1,g2] = boundedline(time_axis, vel_onset_LH, sem_onset_LH,'color',colors_(3,:),'LineWidth',2);
[f1,f2] = boundedline(time_axis, vel_onset_LL, sem_onset_LL,'color',colors_(4,:),'LineWidth',2);
xlabel('Time from corrective saccade (ms)')
ylabel('Velocity (deg/sec)')
title('corrective saccade')
legend([h1, p1, g1, f1], {'HH','HL','LH','LL'});
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'tag_4_velocity';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all   

% endpoint error
fig = figure;
hold on
boundedline(vigor_x_axis, end_error_y_axis.tag_4.tgt_1.rew_1,end_error_sem_y_axis.tag_4.tgt_1.rew_1 ,'color',colors_(1,:),'LineWidth',2)
boundedline(vigor_x_axis, end_error_y_axis.tag_4.tgt_1.rew_0,end_error_sem_y_axis.tag_4.tgt_1.rew_0 ,'color',colors_(2,:),'LineWidth',2)
boundedline(vigor_x_axis, end_error_y_axis.tag_4.tgt_0.rew_1,end_error_sem_y_axis.tag_4.tgt_0.rew_1 ,'color',colors_(3,:),'LineWidth',2)
boundedline(vigor_x_axis, end_error_y_axis.tag_4.tgt_0.rew_0,end_error_sem_y_axis.tag_4.tgt_0.rew_0 ,'color',colors_(4,:),'LineWidth',2)
xlabel('Saccade vigor');
ylabel('End point error (deg)');
ESN_Beautify_Plot(fig, [8, 8], 11);
title('corrective saccade')
file_name = 'tag_4_endpoint_error';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all   


% plot generalized variance  
fig = figure;
hold on
boundedline(vigor_x_axis, gen_var_y_axis.tag_4.tgt_1.rew_1,gen_var_sem_y_axis.tag_4.tgt_1.rew_1 ,'color',colors_(1,:),'LineWidth',2)
boundedline(vigor_x_axis, gen_var_y_axis.tag_4.tgt_1.rew_0,gen_var_sem_y_axis.tag_4.tgt_1.rew_0 ,'color',colors_(2,:),'LineWidth',2)
boundedline(vigor_x_axis, gen_var_y_axis.tag_4.tgt_0.rew_1,gen_var_sem_y_axis.tag_4.tgt_0.rew_1 ,'color',colors_(3,:),'LineWidth',2)
boundedline(vigor_x_axis, gen_var_y_axis.tag_4.tgt_0.rew_0,gen_var_sem_y_axis.tag_4.tgt_0.rew_0 ,'color',colors_(4,:),'LineWidth',2)
xlabel('Saccade vigor');
ylabel('Generalized variance (deg^4)')
ESN_Beautify_Plot(fig, [8, 8], 11);
title('corrective saccade')
file_name = 'tag_4_gen_variance';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all  
%% VIGOR CHOICE TRIALS
tag_idx   = 1; % only applicable for primary saccade
task_id   = 0; % for choice trials
debug_fig = false;
tag_fit_idx = 1;
% Load vigor model
loaded_data = load('C:\Users\Jafar\Documents\reward\data_population\vigor_model.mat');
fit_obj     = loaded_data.fit_obj.(sprintf('m_%s', '132F')).(sprintf('tag_%d', tag_fit_idx));
% Init. data
vigor            = struct;
end_error        = struct;
px_offset        = struct;
py_offset        = struct;
reaction_times   = struct;
eye_vm_vel_onset = struct;
eye_m_vel_visual = struct;
eye_m_vmax       = struct;

% Loop animals
month_dir = dir2(fullfile(data_path,sprintf('data_%s','132F')));
    
% Get vigor model per animal
% fit_params = coeffvalues(fit_obj.(sprintf('m_%s', '132F')).(sprintf('tag_%d', tag_fit_idx)));
    
counter_day = 0; % consider each day as a session
  % Loop months
 num_month = size(month_dir,1);
 for m_idx =  [17 18 19] %1:num_month
     month_path = fullfile(month_dir(m_idx).folder,month_dir(m_idx).name);
     day_dir = dir2(month_path);
     day_dir = day_dir([day_dir.isdir]);
    % Loop days
     num_day = size(day_dir,1);
    
    for d_idx = 1:num_day
        counter_day = counter_day + 1;
        rec_path = fullfile(day_dir(d_idx).folder, day_dir(d_idx).name);
        rec_dir = dir2(rec_path);
        rec_dir = rec_dir([rec_dir.isdir]);
        % Loop recs
        num_rec = size(rec_dir,1);
        for r_idx =1:num_rec
            fprintf('loading session: %s.........\n', rec_dir(r_idx).name);
    
            data_dir = dir2(fullfile(rec_dir(r_idx).folder,rec_dir(r_idx).name, 'analyzed_data\behavior_data\eye','*ANALYZED*'));
            % Load data; if not available, skip the rec
            if isempty(data_dir)
               continue
            end
            loaded_data = load(fullfile(data_dir.folder, data_dir.name));
            SACS_ALL_DATA = loaded_data.sac_data;
            SACS_amp_bin  = discretize(SACS_ALL_DATA.eye_amp_m,  local_amp_edges);
            
            for counter_tgt = 1:length(tgt_cond_list)
                tgt_idx = tgt_cond_list(counter_tgt);

             fprintf('choice %d...\n',  tgt_idx);
             % Init. var
             vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             eye_m_vel_visual.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];
             eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [];

             idx_validity  = SACS_ALL_DATA.validity;
             idx_reaction  = (SACS_ALL_DATA.reaction > min_rxn_time) & (SACS_ALL_DATA.reaction < max_rxn_time);
             idx_tag       = SACS_ALL_DATA.tag == tag_idx;
             idx_task_cond = SACS_ALL_DATA.task_cond == task_id;
             idx_tgt       = SACS_ALL_DATA.choice == tgt_cond_list(counter_tgt);

             if ismember(tag_idx,[1,6])
                idx_amp  = SACS_ALL_DATA.eye_amp_m > tag_1_amp_threshold;
             elseif tag_idx == 4
                 idx_amp = SACS_ALL_DATA.eye_amp_m > tag_4_amp_threshold;
             else
                 idx_amp = SACS_ALL_DATA.eye_amp_m;
             end

           % get the direcion of the eye
           visual_px_offset_low_rew = SACS_ALL_DATA.cue_x_low_rew - SACS_ALL_DATA.start_x; % centering saccade
           visual_py_offset_low_rew = SACS_ALL_DATA.cue_y_low_rew - SACS_ALL_DATA.start_y; % centering saccade
           
           visual_px_offset_high_rew = SACS_ALL_DATA.cue_x_high_rew - SACS_ALL_DATA.start_x; % centering saccade
           visual_py_offset_high_rew = SACS_ALL_DATA.cue_y_high_rew - SACS_ALL_DATA.start_y; % centering saccade
           
           eye_px_offset_t = SACS_ALL_DATA.eye_px_offset - SACS_ALL_DATA.start_x; % centering saccade
           eye_py_offset_t = SACS_ALL_DATA.eye_py_offset - SACS_ALL_DATA.start_y; % centering saccade

           % Distance from targets
           dist_to_tgt_low_rew = sqrt((visual_px_offset_low_rew - eye_px_offset_t).^2 + (visual_py_offset_low_rew - eye_py_offset_t).^2); 
           dist_to_tgt_high_rew = sqrt((visual_px_offset_high_rew - eye_px_offset_t).^2 + (visual_py_offset_high_rew - eye_py_offset_t).^2);
            
           idx_dist_high = dist_to_tgt_high_rew <= 1.5;
           idx_dist_low = dist_to_tgt_low_rew <= 1.5;
            if tgt_idx == 0
                idx_tuned = idx_tag & idx_validity & idx_reaction & idx_task_cond & idx_amp & idx_tgt & idx_dist_low;
            elseif tgt_idx == 1
                idx_tuned = idx_tag & idx_validity & idx_reaction & idx_task_cond & idx_amp & idx_tgt & idx_dist_high;
            end

          
            if tgt_idx == 0
                visual_px_offset = SACS_ALL_DATA.cue_x_low_rew(idx_tuned) - SACS_ALL_DATA.start_x(idx_tuned); % centering saccade
                visual_py_offset = SACS_ALL_DATA.cue_y_low_rew(idx_tuned) - SACS_ALL_DATA.start_y(idx_tuned); % centering saccade

                  % Compute saccade offset error
                offset_error = sqrt((SACS_ALL_DATA.cue_x_low_rew(idx_tuned) - SACS_ALL_DATA.eye_px_offset(idx_tuned)).^2 +...
                (SACS_ALL_DATA.cue_y_low_rew(idx_tuned) - SACS_ALL_DATA.eye_py_offset(idx_tuned)).^2);

            elseif tgt_idx == 1
                visual_px_offset = SACS_ALL_DATA.cue_x_high_rew(idx_tuned) - SACS_ALL_DATA.start_x(idx_tuned); % centering saccade
                visual_py_offset = SACS_ALL_DATA.cue_y_high_rew(idx_tuned) - SACS_ALL_DATA.start_y(idx_tuned); % centering saccade
                  % Compute saccade offset error
                offset_error = sqrt((SACS_ALL_DATA.cue_x_high_rew(idx_tuned) - SACS_ALL_DATA.eye_px_offset(idx_tuned)).^2 +...
                (SACS_ALL_DATA.cue_y_high_rew(idx_tuned) - SACS_ALL_DATA.eye_py_offset(idx_tuned)).^2);
            end
           
            eye_px_onset   = 0; % centering saccade
            eye_py_onset   = 0; % centering saccade
            eye_px_offset = SACS_ALL_DATA.eye_px_offset(idx_tuned) - SACS_ALL_DATA.start_x(idx_tuned); % centering saccade
            eye_py_offset = SACS_ALL_DATA.eye_py_offset(idx_tuned) - SACS_ALL_DATA.start_y(idx_tuned); % centering saccade


            % Get relevant saccade variables; easier to debug this way
            eye_amp_m = SACS_ALL_DATA.eye_amp_m(idx_tuned);
            eye_vm_max = SACS_ALL_DATA.eye_vm_max(idx_tuned);
            % Compute vigor
            mean_vel   = predict(fit_obj,log10(eye_amp_m'));
            vigor_temp = log10(eye_vm_max')./mean_vel;
          
        
            % Compute nominal direction of saccades (deg)
            delta_x = visual_px_offset - eye_px_onset;
            delta_y = visual_py_offset - eye_py_onset;
            visual_ang      = wrapTo360(atan2d(delta_y, delta_x));
            visual_ang_list = unique(visual_ang);
            
            % Compute nominal amp. of saccade & scaling factor to
            % normalize amp.
            nominal_amp = sqrt(delta_x.^2 + delta_y.^2);
%                     amp_scaling = norm_amp./nominal_amp;
            amp_scaling = 1;
            % Compute rotational matrix to rotate saccade vector to
            % right
            rot_eye_px_offset = zeros(1, length(visual_ang));
            rot_eye_py_offset = zeros(1, length(visual_ang));

            for visual_ang_temp = visual_ang_list
                idx_ang = visual_ang == visual_ang_temp;

                R = [cosd(-visual_ang_temp), -sind(-visual_ang_temp);
                     sind(-visual_ang_temp), cosd(-visual_ang_temp)];
              rot_eye_p_offset = R*[eye_px_offset; eye_py_offset];
               
              rot_eye_px_offset(idx_ang) = rot_eye_p_offset(1,idx_ang);
              rot_eye_py_offset(idx_ang) = rot_eye_p_offset(2,idx_ang);
            end
              rot_eye_px_offset = rot_eye_px_offset.*amp_scaling;
              rot_eye_py_offset = rot_eye_py_offset.*amp_scaling;
                %  Validation for rotation (optional)
                if debug_fig 
                    fig = figure;
                    clf;
                    subplot(1,2,1);
                    plot(eye_px_offset, eye_py_offset,'Marker','.','LineStyle','none');
                    xlim([-10,10]);
                    ylim([-10,10]);
                    title('Before');
                    subplot(1,2,2);
                    plot(rot_eye_px_offset, rot_eye_py_offset,'Marker','.','LineStyle','none');
                    xlim([-10,10]);
                    ylim([-10,10]);
                    title('After'); 
                end
              % get the reaction times
            temp_react = SACS_ALL_DATA.reaction(idx_tuned);
            
            % vmax of the tuned saccade
            eye_m_vmax_temp = SACS_ALL_DATA.eye_vm_max(idx_tuned);

            % compute velocity of the saccades
            vm_data   = loaded_data.trials_data.eye_vm_filt;
            time_data = loaded_data.trials_data.time_1K;

            sac_onset_times    = SACS_ALL_DATA.time_onset(idx_tuned);
            visual_onset_times = SACS_ALL_DATA.time_visual(idx_tuned);  
            trials_nums        = SACS_ALL_DATA.trial_num(idx_tuned);
          
            eye_vm_onset_temp    = nan(length_trace,length(sac_onset_times));
            for ii = 1:length(trials_nums)
                id_trial = trials_nums(ii);
                [~,idx1] = min(abs(time_data{1,id_trial} - sac_onset_times(ii)));
                eye_vm_onset_temp(:,ii) = vm_data{1,id_trial}(idx1-length_trace/2:idx1+length_trace/2-1); 
            end
            
            eye_vm_visual_temp   = nan(length_trace*2,length(visual_onset_times));
            for ii = 1:length(trials_nums)
                id_trial = trials_nums(ii);
                [~,idx2] = min(abs(time_data{1,id_trial} - visual_onset_times(ii)));
                if idx2  <= length_trace || idx2 >= length(vm_data{1,id_trial})-length_trace
                    continue
                else
                      eye_vm_visual_temp(:,ii) = vm_data{1,id_trial}(idx2-length_trace:idx2+length_trace-1);
                end
            end

          
                % Concatenate data
                vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), vigor_temp];
                end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), offset_error];
                px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), rot_eye_px_offset];
                py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), rot_eye_py_offset]; 
                reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [reaction_times.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), temp_react];
                eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [eye_vm_vel_onset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), mean(eye_vm_onset_temp,2,'omitnan')]; 
                eye_m_vel_visual.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [eye_m_vel_visual.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), mean(eye_vm_visual_temp,2,'omitnan')];
                eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) = [eye_m_vmax.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), mean(eye_m_vmax_temp,'omitnan')];  
                
            end % tgt loop
         end % rec loop
     end % day loop
end % month loop


vigor_edges   = 0.9:0.05:1.1;
num_vigor_bin = length(vigor_edges) - 1;
vigor_x_axis = vigor_edges(1:end-1) + diff(vigor_edges)./2;      
day_fieldnames = fieldnames(vigor);
num_day        = length(day_fieldnames);
% Init .var
num_sac              = struct;
mean_end_error       = struct;
gen_var              = struct;
var_x                = struct;
end_error_y_axis     = struct;
end_error_sem_y_axis = struct;
gen_var_y_axis       = struct;
gen_var_sem_y_axis   = struct;
var_x_y_axis         = struct;
var_x_sem_y_axis     = struct;

for tag_idx = 1
    for counter_tgt = 1:length(tgt_cond_list)
        tgt_idx = tgt_cond_list(counter_tgt);
         for counter_vigor = 1 : num_vigor_bin
             num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))        = zeros(num_day,1);
             mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor)) = nan(num_day,1);
             gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))        = nan(num_day,1);
             var_x.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))          = nan(num_day,1);
         end
    end
end

% Loop thru sessions/days
for counter_day = 1 : num_day
    for tag_idx = 1
      for counter_tgt = 1:length(tgt_cond_list)
          tgt_idx = tgt_cond_list(counter_tgt);
            vigor_bin = discretize(vigor.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)), vigor_edges);
    
             for counter_vigor = 1 : num_vigor_bin
                 idx_vigor = vigor_bin == counter_vigor;
                 end_error_temp = end_error.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_vigor);
                 px_offset_temp = px_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_vigor);
                 py_offset_temp = py_offset.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_vigor);
    
                 num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1)        = length(end_error_temp);
                 mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1) = mean(end_error_temp);
                 gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1)        = det(cov(px_offset_temp, py_offset_temp));
                 var_x.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor))(counter_day,1)          = var(px_offset_temp);
             end
       end
    end
end

 % Loop thru vigor
for tag_idx = 1
    for counter_tgt = 1:length(tgt_cond_list)
        tgt_idx = tgt_cond_list(counter_tgt);
        for counter_vigor = 1 : num_vigor_bin
            num_sac_temp = num_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor));  
            num_sac_temp_weight = num_sac_temp ./ sum(num_sac_temp);

            mean_end_error_temp = mean_end_error.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor));
    
            end_error_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1, counter_vigor) = mean(mean_end_error_temp.*num_sac_temp_weight, 'omitnan');
            end_error_sem_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1,counter_vigor) = sqrt(var(mean_end_error_temp,num_sac_temp_weight,'omitnan'))./sqrt(sum(num_sac_temp>0));

            gen_var_temp   = gen_var.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor));
            var_x_temp     = var_x.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)).(sprintf('bin_%d',counter_vigor));
            gen_var_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1,counter_vigor) = sum(gen_var_temp.*num_sac_temp_weight, 'omitnan');
            gen_var_sem_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1,counter_vigor) = sqrt(var(gen_var_temp,num_sac_temp_weight,'omitnan'))./sqrt(sum(num_sac_temp>0));

            var_x_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1,counter_vigor) = sum(var_x_temp.*num_sac_temp_weight, 'omitnan');
            var_x_sem_y_axis.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(1,counter_vigor) = sqrt(var(var_x_temp,num_sac_temp_weight,'omitnan'))./sqrt(sum(num_sac_temp>0));
        end
   end
end


% vigor and reaction times for High and Low rewards
vigor_prim_sac          = struct;
reaction_prim_sac       = struct;
vmax_prim_sac           = struct;
velocity_onset_prim_sac = struct;
velocity_visual_prim_sac = struct;

for counter_tgt = 1:length(tgt_cond_list)
   tgt_idx = tgt_cond_list(counter_tgt);
    for idx_d = 1:num_day
      vigor_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) =  nan(num_day,1); 
      reaction_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) =  nan(num_day,1);
      vmax_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) =  nan(num_day,1);
      velocity_onset_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) =  nan(length_trace,num_day);
      velocity_visual_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)) =  nan(2*length_trace,num_day);
    end
end

for counter_tgt = 1:length(tgt_cond_list)
   tgt_idx = tgt_cond_list(counter_tgt);
    for idx_d = 1:num_day
      vigor_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_d,1)           =  median(vigor.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)),'omitnan'); 
      reaction_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_d,1)        =  median(reaction_times.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)),'omitnan');
      vmax_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(idx_d,1)            =   mean(eye_m_vmax.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)),'omitnan');
      velocity_onset_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(:,idx_d)  =   mean(eye_vm_vel_onset.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)),2,'omitnan');
      velocity_visual_prim_sac.(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx))(:,idx_d) =   mean(eye_m_vel_visual.(sprintf('day_%d',idx_d)).(sprintf('tag_%d',tag_idx)).(sprintf('tgt_%d',tgt_idx)),2,'omitnan');
    end
end


colors_ = [0.6350 0.0780 0.1840;0 0.4470 0.7410];
labels = {'High reward';'Low reward'};
fig = figure;
vigor_prim = [vigor_prim_sac.tag_1.tgt_1 vigor_prim_sac.tag_1.tgt_0 ];
violinplot(vigor_prim, labels, 'ViolinColor',colors_);
box off
ylabel('Saccade vigor')
title('Choice saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'choice_vigor_median';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all



 % average of reaction times
RT = histogram(reaction_prim_sac.tag_1.tgt_1,local_RT_edges);
RT.FaceColor = [0.5020 0.0941 0.1686];
hold on
RT = histogram(reaction_prim_sac.tag_1.tgt_0,local_RT_edges);
RT.FaceColor = [0.2157 0.6039 0.8588];
ylabel('Frequency')
xlabel('Reaction time (ms)')
title('Choice saccade')
xlim([ 100 400])
file_name = 'choice_reaction_time_median';
saveas(RT,fullfile(file_path,file_name));
close gcf
fig = openfig(fullfile(file_path,file_name));
ESN_Beautify_Plot(fig, [8, 8], 11);
saveas(fig,fullfile(file_path,file_name),'pdf');
close all


react_prim = [reaction_prim_sac.tag_1.tgt_1 reaction_prim_sac.tag_1.tgt_0 ];
fig = figure;
violinplot(react_prim, labels, 'ViolinColor',colors_);
ylabel('Reaction Times (ms)')
title('Choice saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'choice_reaction_time_pdf';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% MAX velocity

mn_Vmax  = [vmax_prim_sac.tag_1.tgt_1 vmax_prim_sac.tag_1.tgt_0];
fig = figure;
violinplot(mn_Vmax, labels, 'ViolinColor',colors_);
ylabel('Peak velocity (deg/s)')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'choice_peak_velocity';
title('Choice saccade')
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% average velocity saccade onset
speed_high_onset  = mean(velocity_onset_prim_sac.tag_1.tgt_1,2,'omitnan');
sem_sp_high_onset = std(velocity_onset_prim_sac.tag_1.tgt_1,0,2,'omitnan')/sqrt(size(velocity_onset_prim_sac.tag_1.tgt_1,2));

speed_low_onset  = mean(velocity_onset_prim_sac.tag_1.tgt_0,2,'omitnan');
sem_sp_low_onset = std(velocity_onset_prim_sac.tag_1.tgt_0,0,2,'omitnan')/sqrt(size(velocity_onset_prim_sac.tag_1.tgt_0,2));

fig=figure;
vel_x_axis = -length_trace/2:length_trace/2-1;
boundedline(vel_x_axis, speed_high_onset, sem_sp_high_onset,'color',[0.5020 0.0941 0.1686],'LineWidth',2)
hold on
boundedline(vel_x_axis, speed_low_onset, sem_sp_low_onset,'color',[0.2157 0.6039 0.8588],'LineWidth',2)
xlim([-50 100])
xlabel('Time from saccade onset (ms)')
ylabel('Velocity (deg/sec)')
title('Choice saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'choice_velocity_onset';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% average of eye velocity visual onset
speed_high_visual  = mean(velocity_visual_prim_sac.tag_1.tgt_1,2,'omitnan');
sem_sp_high_visual = std(velocity_visual_prim_sac.tag_1.tgt_1,0,2,'omitnan')/sqrt(size(velocity_visual_prim_sac.tag_1.tgt_1,2));

speed_low_visual  = mean(velocity_visual_prim_sac.tag_1.tgt_0,2,'omitnan');
sem_sp_low_visual = std(velocity_visual_prim_sac.tag_1.tgt_0,0,2,'omitnan')/sqrt(size(velocity_visual_prim_sac.tag_1.tgt_0,2));

fig=figure;
vel_x_axis = -length_trace:length_trace-1;
boundedline(vel_x_axis, speed_high_visual, sem_sp_high_visual,'color',[0.5020 0.0941 0.1686],'LineWidth',2)
hold on
boundedline(vel_x_axis, speed_low_visual, sem_sp_low_visual,'color',[0.2157 0.6039 0.8588],'LineWidth',2)
xlabel('Time from visual onset (ms)')
ylabel('Velocity (deg/sec)')
xlim([-100 300])
title('Choice saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'choice_velocity_visual';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% 1- High reward 
fig=figure;
Avg_high = end_error_y_axis.tag_1.tgt_1;
Avg_high(Avg_high==0)=nan;
sem_high= end_error_sem_y_axis.tag_1.tgt_1;
boundedline(vigor_x_axis, Avg_high, sem_high,'color',[0.5020 0.0941 0.1686],'LineWidth',2)
hold on

% 2-Low reward
Avg_low= end_error_y_axis.tag_1.tgt_0;
Avg_low(Avg_low == 0) = nan;
sem_low = end_error_sem_y_axis.tag_1.tgt_0;

boundedline(vigor_x_axis, Avg_low, sem_low,'color',[0.2157 0.6039 0.8588],'LineWidth',2)
xlabel('Saccade vigor');
ylabel('End point error (deg)');
title('choice saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);

file_name = 'choice_endpoint_error';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% Generalized variance
fig=figure;
Avg_high = gen_var_y_axis.tag_1.tgt_1;
Avg_high(Avg_high==0)=nan;
sem_high= gen_var_sem_y_axis.tag_1.tgt_1;
boundedline(vigor_x_axis, Avg_high, sem_high,'color',[0.5020 0.0941 0.1686],'LineWidth',2)
hold on


Avg_low= gen_var_y_axis.tag_1.tgt_0;
Avg_low(Avg_low == 0) = nan;
sem_low = gen_var_sem_y_axis.tag_1.tgt_0;

boundedline(vigor_x_axis, Avg_low, sem_low,'color',[0.2157 0.6039 0.8588],'LineWidth',2)

xlabel('Saccade vigor');
ylabel('Generalized variance (deg^4)')
title('choice saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'choice_gen_variance';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all

% Variance along x axis(saccade direction)

fig=figure;
Avg_high = var_x_y_axis.tag_1.tgt_1;
sem_high = var_x_sem_y_axis.tag_1.tgt_1;
boundedline(vigor_x_axis, Avg_high, sem_high,'color',[0.5020 0.0941 0.1686],'LineWidth',2)
hold on

Avg_low= var_x_y_axis.tag_1.tgt_0;
Avg_low(Avg_low == 0) = nan;
sem_low = var_x_sem_y_axis.tag_1.tgt_0;

boundedline(vigor_x_axis, Avg_low, sem_low,'color',[0.2157 0.6039 0.8588],'LineWidth',2)

xlabel('Saccade vigor');
ylabel('variance along saccade direction (deg^2)')
title('choice saccade')
ESN_Beautify_Plot(fig, [8, 8], 11);
file_name = 'choice_variance_x_axis';
saveas(fig,fullfile(file_path,file_name),'pdf');
close all
%% CHOICE PERFORMANCE
tag_idx   = 1; % only applicable for primary saccade
task_id   = 0; % for choice trials
debug_fig = false;

% Init. data
performace  = struct;

% Loop animals
month_dir = dir2(fullfile(data_path,sprintf('data_%s','132F')));
        
counter_day = 0; % consider each day as a session
  % Loop months
 num_month = size(month_dir,1);
 for m_idx =  17 %1:num_month
     month_path = fullfile(month_dir(m_idx).folder,month_dir(m_idx).name);
     day_dir = dir2(month_path);
    % Loop days
     num_day = size(day_dir,1);
    for d_idx = 1:num_day
        counter_day = counter_day + 1;
        rec_path = fullfile(day_dir(d_idx).folder, day_dir(d_idx).name);
        rec_dir  = dir2(rec_path);
        rec_dir = rec_dir([rec_dir.isdir]);
        % Loop recs
        num_rec = size(rec_dir,1);
        for r_idx =1:num_rec
            fprintf('loading session: %s.........\n', rec_dir(r_idx).name);
    
            data_dir = dir2(fullfile(rec_dir(r_idx).folder,rec_dir(r_idx).name, 'analyzed_data\behavior_data\eye','*ANALYZED*'));
            % Load data; if not available, skip the rec
            if isempty(data_dir)
               continue
            end
            loaded_data   = load(fullfile(data_dir.folder, data_dir.name));
            SACS_ALL_DATA = loaded_data.sac_data;
            SACS_amp_bin  = discretize(SACS_ALL_DATA.eye_amp_m,  local_amp_edges);

             % Init. var
             performace.(sprintf('day_%d',counter_day)).(sprintf('tag_%d',tag_idx)) = [];


             idx_validity  = SACS_ALL_DATA.validity;
             idx_reaction  = (SACS_ALL_DATA.reaction > min_rxn_time) & (SACS_ALL_DATA.reaction < max_rxn_time);
             idx_tag       = SACS_ALL_DATA.tag == tag_idx;
             idx_task_cond = SACS_ALL_DATA.task_cond == task_id;
            
             if ismember(tag_idx,[1,6])
                idx_amp  = SACS_ALL_DATA.eye_amp_m > tag_1_amp_threshold;
             elseif tag_idx == 4
                 idx_amp = SACS_ALL_DATA.eye_amp_m > tag_4_amp_threshold;
             else
                 idx_amp = SACS_ALL_DATA.eye_amp_m;
             end

           % get the direcion of the eye
           visual_px_offset_low_rew = SACS_ALL_DATA.cue_x_low_rew - SACS_ALL_DATA.start_x; % centering saccade
           visual_py_offset_low_rew = SACS_ALL_DATA.cue_y_low_rew - SACS_ALL_DATA.start_y; % centering saccade
           
           visual_px_offset_high_rew = SACS_ALL_DATA.cue_x_high_rew - SACS_ALL_DATA.start_x; % centering saccade
           visual_py_offset_high_rew = SACS_ALL_DATA.cue_y_high_rew - SACS_ALL_DATA.start_y; % centering saccade
           
           eye_px_offset_t = SACS_ALL_DATA.eye_px_offset - SACS_ALL_DATA.start_x; % centering saccade
           eye_py_offset_t = SACS_ALL_DATA.eye_py_offset - SACS_ALL_DATA.start_y; % centering saccade

           % Distance from targets
           dist_to_tgt_low_rew = sqrt((visual_px_offset_low_rew - eye_px_offset_t).^2 + (visual_py_offset_low_rew - eye_py_offset_t).^2); 
           dist_to_tgt_high_rew = sqrt((visual_px_offset_high_rew - eye_px_offset_t).^2 + (visual_py_offset_high_rew - eye_py_offset_t).^2);
            
           idx_dist_high = dist_to_tgt_high_rew <= 1.5;
           idx_dist_low = dist_to_tgt_low_rew <= 1.5;
            if tgt_idx == 0
                idx_tuned = idx_tag & idx_validity & idx_reaction & idx_task_cond & idx_amp & idx_tgt & idx_dist_low;
            elseif tgt_idx == 1
                idx_tuned = idx_tag & idx_validity & idx_reaction & idx_task_cond & idx_amp & idx_tgt & idx_dist_high;
            end

           
            eye_px_onset   = 0; % centering saccade
            eye_py_onset   = 0; % centering saccade
            eye_px_offset = SACS_ALL_DATA.eye_px_offset(idx_tuned) - SACS_ALL_DATA.start_x(idx_tuned); % centering saccade
            eye_py_offset = SACS_ALL_DATA.eye_py_offset(idx_tuned) - SACS_ALL_DATA.start_y(idx_tuned); % centering saccade


            % Get relevant saccade variables; easier to debug this way
            eye_amp_m  = SACS_ALL_DATA.eye_amp_m(idx_tuned);
            eye_vm_max = SACS_ALL_DATA.eye_vm_max(idx_tuned);
          
        
            % Compute nominal direction of saccades (deg)
            delta_x = visual_px_offset - eye_px_onset;
            delta_y = visual_py_offset - eye_py_onset;
            visual_ang      = wrapTo360(atan2d(delta_y, delta_x));
            visual_ang_list = unique(visual_ang);
            
            % Compute nominal amp. of saccade & scaling factor to
            % normalize amp.
            nominal_amp = sqrt(delta_x.^2 + delta_y.^2);
%                     amp_scaling = norm_amp./nominal_amp;
            amp_scaling = 1;
            % Compute rotational matrix to rotate saccade vector to
            % right
            rot_eye_px_offset = zeros(1, length(visual_ang));
            rot_eye_py_offset = zeros(1, length(visual_ang));

            for visual_ang_temp = visual_ang_list
                idx_ang = visual_ang == visual_ang_temp;

                R = [cosd(-visual_ang_temp), -sind(-visual_ang_temp);
                     sind(-visual_ang_temp), cosd(-visual_ang_temp)];
              rot_eye_p_offset = R*[eye_px_offset; eye_py_offset];
               
              rot_eye_px_offset(idx_ang) = rot_eye_p_offset(1,idx_ang);
              rot_eye_py_offset(idx_ang) = rot_eye_p_offset(2,idx_ang);
            end
              rot_eye_px_offset = rot_eye_px_offset.*amp_scaling;
              rot_eye_py_offset = rot_eye_py_offset.*amp_scaling;
                %  Validation for rotation (optional)
                if debug_fig 
                    fig = figure;
                    clf;
                    subplot(1,2,1);
                    plot(eye_px_offset, eye_py_offset,'Marker','.','LineStyle','none');
                    xlim([-10,10]);
                    ylim([-10,10]);
                    title('Before');
                    subplot(1,2,2);
                    plot(rot_eye_px_offset, rot_eye_py_offset,'Marker','.','LineStyle','none');
                    xlim([-10,10]);
                    ylim([-10,10]);
                    title('After'); 
                end
               
         end % rec loop
     end % day loop
end % month loop




















