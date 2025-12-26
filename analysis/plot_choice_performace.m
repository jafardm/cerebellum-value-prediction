function plot_choice_performace(data_path)
    JDM_params_funcs

    user = 'JDM';
    loc  = 'ctx';
    tags = 1:10;

    path_data_monkey_sorted = params.path_data_monkey_sorted;
    animal_list = params.animal_list;
    session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);
    
    num_animal   = numel(session_list);
    data_choice  = cell(num_animal,1);   % session-level performance per animal
    trial_counts = cell(num_animal,1);   % session-level #trials per animal
    chosen_data  = cell(num_animal,1);   % [n_high  n_low] per session

    for counter_animal = 1:num_animal
        disp(animal_list{counter_animal});
        current_path      = path_data_monkey_sorted{counter_animal};
        current_sess_list = session_list{counter_animal};
        num_sess          = numel(current_sess_list);
    
        sess_performance  = nan(num_sess,1);    % choice performance per session
        sess_ntrials      = nan(num_sess,1);    % #choice trials per session
        chosen_trials     = nan(num_sess,2);    % [n_high n_low] per session
      
        parfor counter_sess = 1:num_sess
            current_sess = current_sess_list{counter_sess};
            disp(['     ' current_sess]);
    
            data_eye     = MAF_load_eye_traces(current_path,current_sess);
            units_info   = MAF_extract_cell_metadata(current_path,current_sess);
            cell_num_rec = sum(cell2mat(cellfun(@isempty, ...
                                units_info.cell_ids,'UniformOutput',false)),2);
    
            [~,selected_cell_id] = min(cell_num_rec);
            data = MAF_load_cell(current_path,current_sess,...
                                 units_info.cell_list{selected_cell_id});
    
            rec_flag   = data.rec_info.rec_flag;
            eye_traces = data_eye.eye_traces(rec_flag);
            sacs_tr    = [eye_traces.sac];
            tags_tr    = [sacs_tr.tag]';
            data_rec   = data.data_recordings;   
            eye        = [data_rec.eye];
            sac        = [eye.sac];
            sac_tag_   = [sac.tag]';
            sac_amp_   = [sac.eye_amp]';
    
            assert(all(tags_tr == sac_tag_));
            sac_choice = [sac.choice]';
            task_cond_ = [sac.task_cond]'; 
            sac_rt_    = (([sac.time_onset] - [sac.time_visual])*1e3)'; % ms
       
            idx_rxt      = sac_rt_ > 0 & sac_rt_ < 600;
            tag_1        = sac_tag_ == 1;
            sac_amp_idx  = sac_amp_ > 3.5;
            assert(length(sac_rt_) == length(sac_tag_));

            res          = MAF_find_sac_outliers(sac);
            ind_selected = (ismember(sac_tag_,tags) & ~res.rm_ind' & idx_rxt); % not outliers + tag 1–10
            idx_tuned    = ind_selected & tag_1 & sac_amp_idx;
        
            sess_task_cond = task_cond_(idx_tuned);
            sess_tag       = tag_1(idx_tuned);
            sess_cho       = sac_choice(idx_tuned);
           
            conditions = [sess_task_cond, sess_tag, sess_cho];

            trial_pecr = sum(sess_task_cond == 0) / ...
                (sum(sess_task_cond == 0) +  sum(sess_task_cond == 1));
            fprintf('\n---> Choice percentage %.2f. \n',trial_pecr)

            idx_high = ismember(conditions, [0,1,1], 'rows'); % choice H
            idx_low  = ismember(conditions, [0,1,0], 'rows'); % choice L

            % Count number of high/low choices
            n_high = sum(idx_high);
            n_low  = sum(idx_low);
            n_tot  = n_high + n_low;

            if n_tot > 0
                sess_performance(counter_sess) = n_high / n_tot;
                chosen_trials(counter_sess,:)  = [n_high  n_low];
            else
                sess_performance(counter_sess) = NaN;
                chosen_trials(counter_sess,:)  = [NaN NaN];
            end
            sess_ntrials(counter_sess) = n_tot;

            fprintf('---> Choice performance %.2f%% (N=%d trials).\n', ...
                    100 * sess_performance(counter_sess), n_tot);
        end  

        data_choice{counter_animal}  = sess_performance;
        trial_counts{counter_animal} = sess_ntrials;
        chosen_data{counter_animal}  = chosen_trials;
    end

    %% Plotting
    fig = figure;

    colors = [0.00 0.45 0.74; 
              0.85 0.33 0.10]; 

    % -------- Top row: performance across sessions (spans both columns) ---
    subplot(2,2,[1 2])
    hold on;

    for subj_idx = 1:2
        subj_perf   = data_choice{subj_idx};
        n_sessions  = numel(subj_perf);
        plot(1:n_sessions, subj_perf * 100, 'o-', ...
            'DisplayName', sprintf('%s (n=%d sessions)', ...
                                   animal_list{subj_idx}, n_sessions), ...
            'Color', colors(subj_idx,:), ...
            'LineWidth', 1.5);
    end

    xlabel('Session');
    ylabel('Choice Performance (%)');
    title('Choice Performance Across Sessions');
    legend('Location', 'best');
    ylim([60 100]);
    box off

        % -------- Bottom left: 132F – P(High/Low per session) ----------
    if num_animal >= 1
        subplot(2,2,3); cla; hold on;
        counts = chosen_data{1};    % [n_high n_low] per session

        % Per-session total trials
        n_tot_sess = nansum(counts, 2);

        % Compute probabilities per session
        p_mat = nan(size(counts));  % [P_high P_low] per session
        valid = n_tot_sess > 0;
        p_mat(valid,1) = counts(valid,1) ./ n_tot_sess(valid); % P(High)
        p_mat(valid,2) = counts(valid,2) ./ n_tot_sess(valid); % P(Low)

        % Convert to percent for plotting
        p_mat = p_mat * 100;

        % Violin plot: each column is a distribution across sessions
        violinplot(p_mat, {'High','Low'});
        ylabel('Choice probability per session (%)');
        title(sprintf('%s: P(High / Low | session)', animal_list{1}));
        ylim([0 100]);
        box off

        % ---- Add total N labels (number of trials) ----
        total_high = nansum(counts(:,1));
        total_low  = nansum(counts(:,2));

        y_limits = ylim;
        y_text   = y_limits(2) - 0.05 * range(y_limits);

        text(1, y_text, sprintf('N=%d', total_high), ...
             'HorizontalAlignment','center', 'FontSize', 9);
        text(2, y_text, sprintf('N=%d', total_low), ...
             'HorizontalAlignment','center', 'FontSize', 9);
    end

    % -------- Bottom right: 65F – P(High/Low per session) ----------
    if num_animal >= 2
        subplot(2,2,4); cla; hold on;
        counts = chosen_data{2};    % [n_high n_low] per session

        n_tot_sess = nansum(counts, 2);

        p_mat = nan(size(counts));
        valid = n_tot_sess > 0;
        p_mat(valid,1) = counts(valid,1) ./ n_tot_sess(valid); % P(High)
        p_mat(valid,2) = counts(valid,2) ./ n_tot_sess(valid); % P(Low)

        p_mat = p_mat * 100;

        violinplot(p_mat, {'High','Low'});
        ylabel('Choice probability per session (%)');
        title(sprintf('%s: P(High / Low | session)', animal_list{2}));
        ylim([0 100]);
        box off

        total_high = nansum(counts(:,1));
        total_low  = nansum(counts(:,2));

        y_limits = ylim;
        y_text   = y_limits(2) - 0.05 * range(y_limits);

        text(1, y_text, sprintf('N=%d', total_high), ...
             'HorizontalAlignment','center', 'FontSize', 9);
        text(2, y_text, sprintf('N=%d', total_low), ...
             'HorizontalAlignment','center', 'FontSize', 9);
    end
   
    
    ESN_Beautify_Plot(fig, [12, 8]);
    fname = 'choice performace.pdf';
    img_save_path = fullfile(data_path, 'population_figs', 'behavior');
     print(fig, fullfile(img_save_path, fname), '-dpdf', '-bestfit');
end
