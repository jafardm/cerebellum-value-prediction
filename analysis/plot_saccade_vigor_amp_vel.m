function plot_saccade_vigor_amp_vel(data_path)
    JDM_params_funcs;
    S = load(fullfile(data_path, 'population_data', 'behave_data.mat'));
    if isfield(S,'data_behav'), data_behav = S.data_behav; else, error('behave_data missing'); end

    alpha = .9;
    sz    = .2;

    amp_edges = params.pop.sac.amp_edges;
    vel_edges = params.pop.sac.vel_edges;
    subject_names = params.animal_list;

    % Two tag groups
    tag_cond = {10, [1,4,6,7]};
    cond_labels = {'Tag 10','Tag 1,4,6,7'}; 
    clr = {'b','r'};  

    fig = figure('Name','Saccade Amplitude–Velocity (log–log)');
    ncol = 2; 
    nrow = ceil(numel(subject_names) / ncol);
    tl = tiledlayout(fig, nrow, ncol, 'TileSpacing','compact', 'Padding','compact');

    xL = [-1, 1.5];              
    yTicks = 1:3;                
    yL = [yTicks(1), yTicks(end)];

    for counter_animal = 1:numel(subject_names)
        nexttile; hold on; box on;

        current_sac_amp = data_behav.sac_amp{counter_animal};
        current_sac_vel = data_behav.sac_vel{counter_animal};
        current_sac_tag = data_behav.sac_tag{counter_animal};

        valid_all = isfinite(current_sac_amp) & isfinite(current_sac_vel) & ...
                    current_sac_amp > 0 & current_sac_vel > 0;

        b_ms = gmregress(log10(current_sac_amp(valid_all)), log10(current_sac_vel(valid_all)));

        % --- Plot tag groups
        for k = 1:numel(tag_cond)
            ind_tag = ismember(current_sac_tag, tag_cond{k}) & valid_all;
            if any(ind_tag)
                scatter(log10(current_sac_amp(ind_tag)), ...
                        log10(current_sac_vel(ind_tag)), ...
                        sz, 'filled', clr{k}, 'MarkerFaceAlpha', alpha, ...
                        'DisplayName', cond_labels{k});
            end
        end

        % --- Regression line (excluded from legend)
        x_ = linspace(xL(1), xL(2), 100);
        plot(x_, x_*b_ms(2) + b_ms(1), 'k', 'LineWidth', 1, 'HandleVisibility','off');

        % --- Bin-edge guides
        if numel(amp_edges) > 2
            xline(log10(amp_edges(2:end-1)), '--k', 'HandleVisibility','off');
        end
        if numel(vel_edges) > 2
            yline(log10(vel_edges(2:end-1)), '--k', 'HandleVisibility','off');
        end

        % --- Axes formatting
        xlim(xL); ylim(yL);
        xticks(xL(1):xL(2));
        yticks(yTicks);
        xticklabels(arrayfun(@(v) sprintf('%g', 10.^v), xticks, 'UniformOutput', false));
        yticklabels(arrayfun(@(v) sprintf('%g', 10.^v), yticks, 'UniformOutput', false));

        xlabel('amp (deg)');
        ylabel('vel (deg/s)');
        title(subject_names{counter_animal});

        % Legend only once per row (or only for first panel)
        if counter_animal == 1
            legend('Location','best'); legend boxoff;
        end
    end

    title(tl, 'Saccade Amplitude–Velocity (log–log)');
    out_dir = fullfile(data_path, 'population_figs', 'behavior');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    ESN_Beautify_Plot(fig, [10, 4]);

    export_fig(fig, fullfile(out_dir,'vigor_amp_vel'), '-pdf', '-painters');
end
