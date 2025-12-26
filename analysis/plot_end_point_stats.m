function plot_end_point_stats(data_path)

    img_save_path = fullfile(data_path, 'population_figs', 'behavior');
   
    JDM_params_funcs;
    load(fullfile(data_path, 'population_data', 'stat_data_enderror_prim_sac.mat'));

    subj_names   = params.animal_list;
    num_subjects = numel(subj_names);

    %% Collect variance differences (High − Low)
    all_diff_x = {};
    all_diff_y = {};

    for jj = 1:num_subjects
        var_x_high = data_behav.var_x_high{jj, 1};
        var_y_high = data_behav.var_y_high{jj, 1};

        var_x_low  = data_behav.var_x_low{jj, 1};
        var_y_low  = data_behav.var_y_low{jj, 1};

        diff_x_axis = var_x_high - var_x_low;
        diff_y_axis = var_y_high - var_y_low;

        all_diff_x{jj} = diff_x_axis;
        all_diff_y{jj} = diff_y_axis;
    end

    %% Create figure with two subplots
    fig = figure;
    set(gcf, 'Position', [100 100 1000 400]);

    %% ---- X-axis plot ----
    subplot(1,2,1); hold on;
    all_vals_x = [];
    all_labels_x = [];
    x_pvals = zeros(num_subjects,1);
    x_ns = zeros(num_subjects,1);

    for jj = 1:num_subjects
        diff_vals = all_diff_x{jj};
        all_vals_x   = [all_vals_x; diff_vals(:)];
        all_labels_x = [all_labels_x; repmat(subj_names(jj), numel(diff_vals), 1)];
        [p, ~] = signrank(diff_vals);
        x_pvals(jj) = p;
        x_ns(jj)    = numel(diff_vals);
    end

    violinplot(all_vals_x, all_labels_x);
    ylabel('Endpoint variance (High − Low), saccade axis');
    set(gca, 'FontSize', 12);
    yline(0, '--k');

    % Add statistics to title
    stats_txt = cell(1, num_subjects);
    for jj = 1:num_subjects
        if x_pvals(jj) < 0.0001
            p_text = 'p<0.0001';
        else
            p_text = sprintf('p=%.4f', x_pvals(jj));
        end
        stats_txt{jj} = sprintf('%s: %s (n=%d)', subj_names{jj}, p_text, x_ns(jj));
    end
    title({'Endpoint Variance Difference: X-axis'; strjoin(stats_txt, ' | ')});

    %% ---- Y-axis plot ----
    subplot(1,2,2); hold on;
    all_vals_y = [];
    all_labels_y = [];
    y_pvals = zeros(num_subjects,1);
    y_ns = zeros(num_subjects,1);

    for jj = 1:num_subjects
        diff_vals = all_diff_y{jj};
        all_vals_y   = [all_vals_y; diff_vals(:)];
        all_labels_y = [all_labels_y; repmat(subj_names(jj), numel(diff_vals), 1)];
        [p, ~] = signrank(diff_vals);
        y_pvals(jj) = p;
        y_ns(jj)    = numel(diff_vals);
    end

    violinplot(all_vals_y, all_labels_y);
    ylabel('Endpoint variance (High − Low), perpendicular axis');
    set(gca, 'FontSize', 12);
    yline(0, '--k');

    % Add statistics to title
    stats_txt = cell(1, num_subjects);
    for jj = 1:num_subjects
        if y_pvals(jj) < 0.0001
            p_text = 'p<0.0001';
        else
            p_text = sprintf('p=%.4f', y_pvals(jj));
        end
        stats_txt{jj} = sprintf('%s: %s (n=%d)', subj_names{jj}, p_text, y_ns(jj));
    end
    title({'Endpoint Variance Difference: Y-axis'; strjoin(stats_txt, ' | ')});

    %% Save figure
    fname = 'high-vs-low-endpoint-var-stats-prim-sac.pdf';
    ESN_Beautify_Plot(fig, [12, 5]);
    print(fig, fullfile(img_save_path, fname), '-dpdf', '-bestfit');

end
