function plot_behave_by_trilas(data_path)

    % Load parameters
    JDM_params_funcs;
    subj_names = params.animal_list;
    num_subjects = numel(subj_names);
    
    % Paths and data
    img_save_path = fullfile(data_path, 'population_figs', 'behavior');
    load(fullfile(data_path, 'population_data', 'Behave_trial_by_trial'));

    % Define distinct colors
    Line_Color = [
        0.00, 0.45, 0.74;  % blue
        0.85, 0.33, 0.10;  % orange
        0.47, 0.67, 0.19;  % green
        0.49, 0.18, 0.56   % purple
    ];
    
    labels = {'HH→HX', 'HL→HX', 'LL→LX', 'LH→LX'};

    fig = figure;
    for iA = 1:num_subjects
        subplot(1, 2, iA)
        hold on
        all_data = [];
        group_labels = [];

        rt_data = data_behav.sac_rt{iA};  % Each is a 1x4 cell of vectors
        
        for i = 1:4
            this_rt = rt_data{i};
            all_data = [all_data; this_rt];
            group_labels = [group_labels; repmat(labels(i), numel(this_rt), 1)];
        end
        
        group_labels = categorical(group_labels, labels, 'Ordinal', true);

        % Violin plot
        v = violinplot(all_data, group_labels);

        for i = 1:length(v)
            v(i).ViolinColor = {Line_Color(i,:)};
        end
        if iA == 1
            ylabel('RT difference (ms)');
        end
        title(sprintf('Subject %s', subj_names{iA}));
        xtickangle(45);
        box off
        hold off
    end
    
    sgtitle('Primary saccade reaction Trial changes');
    ESN_Beautify_Plot(fig, [8, 6]);
%     print(fig, fullfile(img_save_path, 'RT_diff_by_trial.pdf'), '-dpdf', '-bestfit');
    %% ------- Vmax ------

    fig = figure;
    for iA = 1:num_subjects
        subplot(1, 2, iA)
        hold on
        all_data = [];
        group_labels = [];

        data_ = data_behav.vmax{iA};  % Each is a 1x4 cell of vectors
        
        for i = 1:4
            this_rt = data_{i};
            all_data = [all_data; this_rt];
            group_labels = [group_labels; repmat(labels(i), numel(this_rt), 1)];
        end
        
        group_labels = categorical(group_labels, labels, 'Ordinal', true);

        % Violin plot
        v = violinplot(all_data, group_labels);

        for i = 1:length(v)
            v(i).ViolinColor = {Line_Color(i,:)};
        end
        if iA==1
            ylabel('Max velocity (deg/sec)');
        end
        title(sprintf('Subject %s', subj_names{iA}));
        xtickangle(45);
        box off
        hold off
    end
    
    sgtitle('Primary saccade Vmax Trial changes');
    ESN_Beautify_Plot(fig, [8, 6]);
%     print(fig, fullfile(img_save_path, 'vmax_diff_by_trial.pdf'), '-dpdf', '-bestfit');
 %% ------- end-erorr ------
    fig = figure;
    for iA = 1:num_subjects
        subplot(1, 2, iA)
        hold on
        all_data = [];
        group_labels = [];

        data_ = data_behav.sac_error{iA};  % Each is a 1x4 cell of vectors
        
        for i = 1:4
            this_rt = data_{i};
            all_data = [all_data; this_rt];
            group_labels = [group_labels; repmat(labels(i), numel(this_rt), 1)];
        end
        
        group_labels = categorical(group_labels, labels, 'Ordinal', true);

        % Violin plot
        v = violinplot(all_data, group_labels);

        for i = 1:length(v)
            v(i).ViolinColor = {Line_Color(i,:)};
        end
        if iA ==1
           ylabel('End point error(deg)');
        end
        title(sprintf('Subject %s', subj_names{iA}));
        xtickangle(45);
        box off
        hold off
    end
    
    sgtitle('Primary saccade End-error changes');
    ESN_Beautify_Plot(fig, [8, 6]);
%     print(fig, fullfile(img_save_path, 'Enderror_diff_by_trial.pdf'), '-dpdf', '-bestfit');
end
