function plot_cs_rate_irrelevant_sac_vel_bins(data_path)
% CS rates for task-irrelevant saccades, aligned to CS-on direction.
% One subplot per velocity bin:
%   Left y-axis: mean ± SEM CS rate for CS-on-framed direction groups
%   Right y-axis: mean eye velocity (area)
%
% Requires: population_data/cs_on_rate_irrelevant_onset_sac.mat
% Fields: cs_on_data, cs_rate_by_vel, cs_vm_by_vel, n_trials_by_vel, vel_edges_out, params

    JDM_params_funcs;
    S = load(fullfile(data_path,'population_data','cs_on_rate_irrelevant_onset_sac'));

    % --- pull from saved file ---
    cs_on_data       = S.cs_on_data;
    cs_rate_by_vel   = S.cs_rate_by_vel;    % [num_cs x 500 x 8 x nVel]
    cs_vm_by_vel     = S.cs_vm_by_vel;      % [num_cs x 500 x 8 x nVel]
    n_trials_by_vel  = S.n_trials_by_vel;   % [num_cs x 8 x nVel]
    if isfield(S,'vel_edges_out')
        vel_edges = S.vel_edges_out;
    else
        vel_edges = params.pop.sac.vel_edges;
    end

    % --- compute CS-on direction for each cell (stronger of vis vs sac) ---
    cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_data, 'UniformOutput', false));
    cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_data, 'UniformOutput', false));
    cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, cs_on_data, 'UniformOutput', false));
    cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, cs_on_data, 'UniformOutput', false));

    num_cs  = numel(cs_on_data);
    use_sac = cs_on_rho_sac > cs_on_rho_vis;
    cs_on_ang_deg         = cs_on_ang_vis;
    cs_on_ang_deg(use_sac)= cs_on_ang_sac(use_sac);
    cs_on_ang = cs_on_ang_deg/180*pi;   % radians, [num_cs x 1]

    % --- build per-cell rotation order so dir #1 = CS-on ---
    ang_bins = params.sac.ang_values(:)'/180*pi; % 1x8 radians
    order8   = nan(num_cs, 8);
    for c = 1:num_cs
        % robust circular difference without angdiff()
        d = abs(atan2(sin(ang_bins - cs_on_ang(c)), cos(ang_bins - cs_on_ang(c)))); % 1x8
        [~, cs_on_bin] = min(d);
        order8(c,:) = circshift(1:8, 1 - cs_on_bin);
    end

    % --- rotate traces to CS-on frame ---
    [~, T, ~, nVel] = size(cs_rate_by_vel);
    cs_rate_by_vel_rot = nan(num_cs, T, 8, nVel);
    cs_vm_by_vel_rot   = nan(num_cs, T, 8, nVel);
    for c = 1:num_cs
        for v = 1:nVel
            cs_rate_by_vel_rot(c,:,:,v) = cs_rate_by_vel(c,:,order8(c,:),v);
            cs_vm_by_vel_rot(c,:,:,v)   = cs_vm_by_vel(c,:,  order8(c,:),v);
        end
    end

    % --- plotting params ---
    dir_groups = { 1, [2 8], [3 7], [4 6], 5 }; % CS-on, ±45, ±90, ±135, 180
    leg_dirs   = {'cs-on','cs±45','cs±90','cs±135','cs±180'};
    colors_    = flipud(cool(5));
    vm_color   = [204, 174, 98]/256;

    vm_lim     = [0 400];
    y_cs_lim   = [0.6 2];
    t_ind_sac  = (-100:150) + 250; % time indices (t=0 at 250)
    t          = t_ind_sac - 250;

    % total trial counts per velocity bin (sum over cells & dirs)
    Nbin_total = squeeze(nansum(nansum(n_trials_by_vel,2),1));  % [nVel x 1]

    % --- figure with one tile per velocity bin ---
    fig1 = figure;
    tl = tiledlayout(1, nVel, 'TileSpacing','compact','Padding','compact');

    % collect legend line handles (only from first panel)
    hLeg = gobjects(numel(dir_groups),1);

    for v = 1:nVel
        ax = nexttile(v); hold(ax,'on'); box(ax,'on');

        % bin label
        if v < nVel
            bin_str = sprintf('%d–%d deg/s', vel_edges(v), vel_edges(v+1));
        else
            bin_str = sprintf('%d–≤%d deg/s', vel_edges(v), vel_edges(v+1));
        end

        % --- plot CS rate groups ---
        for g = numel(dir_groups):-1:1
            idx_dirs = dir_groups{g};
            cur_rate = squeeze(mean(cs_rate_by_vel_rot(:, t_ind_sac, idx_dirs, v), 3, 'omitnan')); % [num_cs x time]
            mu  = mean(cur_rate, 1, 'omitnan');
            N   = sum(~isnan(cur_rate), 1);
            se  = std(cur_rate, 0, 1, 'omitnan') ./ max(1, sqrt(N));

            yyaxis(ax,'left')
            [hl1,hl2] = boundedline(t, mu, se, 'alpha', 'cmap', colors_(g,:));
            if v == 1
                hl1.DisplayName = leg_dirs{g};
                hLeg(g) = hl1; % save handle for legend
            else
                hl1.HandleVisibility = 'off';
            end
            hl2.HandleVisibility = 'off';
        end

        % --- velocity overlay for this bin ---
        cur_vm = squeeze(mean(cs_vm_by_vel_rot(:, t_ind_sac, :, v), 3, 'omitnan')); % [num_cs x time]
        mu_vm  = mean(cur_vm, 1, 'omitnan');

        xline(ax,0,'--','HandleVisibility','off');

        yyaxis(ax,'right')
        area(ax, t, mu_vm, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
            'EdgeColor','none','HandleVisibility','off');
        ylim(ax, vm_lim);
        ylabel(ax,'Velocity (deg/s)');
        ax.YColor = 'k'; % right axis color

        yyaxis(ax,'left')
        ylim(ax, y_cs_lim);
        if v == 1, ylabel(ax,'CS rate (spk/s)'); end
        xlabel(ax,'time from saccade onset (ms)');

        title(sprintf('%s  |  N=%d', bin_str, Nbin_total(v)));
    end

    % shared title + legend
    title(tl, sprintf('Irrelevant saccades — CS cells (n = %d)', num_cs));
    if all(isgraphics(hLeg))
        legend(hLeg, leg_dirs, 'Location','northeast'); legend boxoff;
    end

    ESN_Beautify_Plot(fig1, [max(6, 4.5*nVel), 4]);

    out_dir = fullfile(data_path,'population_figs','ephys');
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    print(fig1, fullfile(out_dir,'cs_rate_irrelevant_sac_byVelocity.pdf'), '-dpdf','-bestfit');
end
