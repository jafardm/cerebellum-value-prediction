function plot_end_point_error_contour_choice(data_path)

img_save_path = fullfile(data_path, 'population_figs', 'behavior');

JDM_params_funcs;
load(fullfile(data_path, 'population_data', 'sac_kinematic_prim_sac.mat'));

subj_names  = params.animal_list;
num_subjects = numel(subj_names);

group_labels = {'H-chosen', 'L-chosen'};
Line_Color   = [190, 37, 43; 0, 10, 255] / 255;

%% Loop through subjects
for jj = 1:num_subjects
    data_error = tag1_end_point_error{jj,1}(5:6);
    % High = 5, Low = 6
 
    x_lim = [-2.5 1.1];
    y_lim = [-1.5 1.5];

    % KDE grid
    num_p = 200;
    [xi, yi] = meshgrid(linspace(x_lim(1), x_lim(2), num_p), ...
                        linspace(y_lim(1), y_lim(2), num_p));
    dx = xi(1,2) - xi(1,1);
    dy = yi(2,1) - yi(1,1);

    % -------- Shared bandwidth (High & Low pooled) --------
    pooled = [];
    for i = 1:2
        xy = data_error{i};
        if ~isempty(xy)
            pooled = [pooled; xy]; %#ok<AGROW>
        end
    end
    use_bw = false;
    if ~isempty(pooled)
        [~,~,bw] = ksdensity(pooled);   % 1x2 bandwidth for 2-D
        use_bw = true;
    end
    % ------------------------------------------------------

    % Compute KDE -> probability per grid cell
    f_cell_all = cell(1,2);
    max_cell = 0;
    for i = 1:2
        xy = data_error{i};
        if isempty(xy)
            f_cell_all{i} = [];
            continue;
        end

        if use_bw
            f = ksdensity(xy, [xi(:), yi(:)], 'Bandwidth', bw); % density: 1/deg^2
        else
            f = ksdensity(xy, [xi(:), yi(:)]);                  % density: 1/deg^2
        end
        f = reshape(f, size(xi));
        f_cell = f * dx * dy;                                   % probability per grid cell (unitless)
        f_cell_all{i} = f_cell;
        max_cell = max(max_cell, max(f_cell(:)));

        % Sanity check (should be ~1 if the grid covers the mass)
        approx_mass = sum(f_cell(:));
        fprintf('Subject %s - %s mass ~ %.3f\n', subj_names{jj}, group_labels{i}, approx_mass);
    end

    % Figure with shared colorbar
    fig = figure('Position', [100, 100, 800, 400]);
    colormap(parula);

    ax_cb = axes('Position', [0.91 0.22 0.03 0.56]); 
    caxis([0 max_cell]); 
    colorbar(ax_cb);
    ylabel(ax_cb, 'Probability per grid cell');
    axis(ax_cb, 'off');

    for i = 1:2
        xy = data_error{i};
        f_cell = f_cell_all{i};
        if isempty(xy) || isempty(f_cell), continue; end

        subplot(1,2,i);
        set(gca, 'Color', 'w'); hold on; axis equal;
        xlim(x_lim); ylim(y_lim);
        xlabel('X Error (deg)'); ylabel('Y Error (deg)');
        title(sprintf('%s - %s', subj_names{jj}, group_labels{i}), 'FontSize', 12);
        % Red crosshair at (0,0)
        plot([0 0], y_lim, 'r-', 'LineWidth', 1.5);   % vertical red line
        plot(x_lim, [0 0], 'r-', 'LineWidth', 1.5);   % horizontal red line


        % Heatmap: probability per grid cell
        surf(xi, yi, ones(size(f_cell)) * 1e-6, f_cell, ...
             'EdgeColor', 'none', 'FaceAlpha', 0.6);
        caxis([0 max_cell]); view(2);

        % Custom grid overlay
        xt = get(gca,'XTick'); yt = get(gca,'YTick');
        for k = 1:numel(xt)
            plot([xt(k) xt(k)], y_lim, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        end
        for k = 1:numel(yt)
            plot(x_lim, [yt(k) yt(k)], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        end

        % Contours at mass levels (10%, 25%, 50%, 75%, 90%, 95%)
        levels = [0.10 0.25 0.50 0.75 0.90 0.95];
        flat   = sort(f_cell(:), 'descend');
        cum    = cumsum(flat) / sum(flat);
        thresholds = arrayfun(@(q) flat(find(cum >= q, 1, 'first')), levels);
        contour(xi, yi, f_cell, thresholds, 'LineColor', 'k', 'LineWidth', 1.5);

        % Median marker
        med_x = median(xy(:,1), 'omitnan');
        med_y = median(xy(:,2), 'omitnan');
        scatter(med_x, med_y, 8, 'o', 'MarkerFaceColor', 'w', ...
                'MarkerEdgeColor', Line_Color(i,:), 'LineWidth', 2);
    end

    % Save figure
    fname = sprintf('%s-choice-endpoint-density.pdf', subj_names{jj});
    ESN_Beautify_Plot(fig, [10, 4]);
    print(fig, fullfile(img_save_path, fname), '-dpdf', '-bestfit');
end
end
