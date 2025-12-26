function plot_cs_rate_irrelevant_sac(data_path)
% CS rates for task-irrelevant saccades, aligned to CS-on direction.
% Left y-axis: mean ± SEM CS rate; Right y-axis: mean eye velocity (area).
%
% Expects file: population_data/cs_on_rate_irrelevant_onset_sac.mat
% providing: cs_on_data, cs_rate [num_cs x T x 8], cs_vm [num_cs x T x 8], params

    JDM_params_funcs;
    S = load(fullfile(data_path,'population_data','cs_on_rate_irrelevant_onset_sac'));
    % pull into workspace
    cs_on_data = S.cs_on_data;
    cs_rate    = S.cs_rate;         % [num_cs x 500 x 8] (onset-aligned)
    cs_vm      = S.cs_vm;           % [num_cs x 500 x 8]
 
    % ----- compute CS-on angle (pick stronger of vis vs sac) -----
    cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_data, 'UniformOutput', false));
    cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_data, 'UniformOutput', false));

    cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, cs_on_data, 'UniformOutput', false));
    cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, cs_on_data, 'UniformOutput', false));

    num_cs = numel(cs_on_data);
    use_sac = cs_on_rho_sac > cs_on_rho_vis;

    cs_on_ang_deg = cs_on_ang_vis;
    cs_on_ang_deg(use_sac) = cs_on_ang_sac(use_sac);    % degrees (0..360)
    cs_on_ang = cs_on_ang_deg/180*pi;                   % radians

    % ----- build per-cell rotation order so dir #1 = CS-on -----
    ang_bins = params.sac.ang_values/180*pi;            % 8 nominal directions (radians)
    [~, cs_on_bin] = min(abs(angdiff(repmat(ang_bins, num_cs, 1), cs_on_ang*ones(size(ang_bins)))), [], 2);

    order8 = nan(num_cs, 8);
    for c = 1:num_cs
        % circular shift so index 1 is the CS-on bin for that cell
        order8(c,:) = circshift(1:8, 1 - cs_on_bin(c));
    end

    % ----- rotate time series to CS-on frame -----
    [~, T, ~] = size(cs_rate);
    cs_rate_rot = nan(num_cs, T, 8);
    cs_vm_rot   = nan(num_cs, T, 8);
    for c = 1:num_cs
        cs_rate_rot(c,:,:) = cs_rate(c,:,order8(c,:));
        cs_vm_rot(c,:,:)   = cs_vm(c,:,  order8(c,:));
    end

    % ----- plotting params -----
    leg_dirs = {'cs-on','cs±45','cs±90','cs±135','cs±180'};  % groups
    colors_  = flipud(cool(5));
    vm_color = [204, 174, 98]/256;
    vm_lim   = [0 400];
    t_ind_sac = (-249:250) + 250;    % indices; t=0 at 250
    t = t_ind_sac - 250;

    % mapping for symmetric direction groups (in CS-on frame)
    dir_groups = {
        1,          ... % cs-on
        [2 8],      ... % ±45
        [3 7],      ... % ±90
        [4 6],      ... % ±135
        5           ... % 180
    };

    fig1 = figure;
    hold on; box on;

    % draw from far angles first so cs-on overlays on top
    for g = numel(dir_groups):-1:1
        idx_dirs = dir_groups{g};

        % average across the symmetric directions first -> [num_cs x T]
        cur_rate = squeeze(mean(cs_rate_rot(:, t_ind_sac, idx_dirs), 3, 'omitnan'));
        % population mean/SEM across cells
        mu  = mean(cur_rate, 1, 'omitnan');
        N   = sum(~isnan(cur_rate), 1);
        se  = std(cur_rate, 0, 1, 'omitnan') ./ max(1, sqrt(N));

        % CS rate (left y-axis)
        yyaxis left
        [hl1,hl2] = boundedline(t, mu, se, 'alpha', 'cmap', colors_(g,:));
        hl1.DisplayName = leg_dirs{g};
        hl2.HandleVisibility = 'off';
        hold on
    end

    % velocity overlay once (mean across all 8 dirs, then across cells)
    cur_vm = squeeze(mean(cs_vm_rot(:, t_ind_sac, :), 3, 'omitnan'));   % [num_cs x time]
    mu_vm  = mean(cur_vm, 1, 'omitnan');

    xline(0,'--','HandleVisibility','off');
    xlabel('time from saccade onset (ms)')

    yyaxis right
    area(t, mu_vm, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
        'EdgeColor', 'none', 'HandleVisibility','off');
    ylim(vm_lim);
    ylabel('Velocity (deg/s)');
    ax = gca; 
    ax.YColor = 'k';

    yyaxis left
    ylim([0.6 2])
    ylabel('CS rate (spk/s)');
    title([ 'Irrelevant saccades- ' sprintf('CS Cells (n = %d)', num_cs) ])
    legend('Location','northeast'); legend boxoff;

    ESN_Beautify_Plot(fig1,[10,4])

    out_dir = fullfile(data_path,'population_figs','ephys');
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    print(fig1, fullfile(out_dir,'cs_rate_irrelevant_sac.pdf'), '-dpdf', '-bestfit');
end
