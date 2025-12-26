function plot_figure_1(data_path)
    % --------------------------
    % 1) COMPUTE EVERYTHING ONCE
    % --------------------------
    CS_on_data = compute_CS_ON(data_path);
  
    make_figC1_all_dirs(CS_on_data, data_path);
    make_figC2_high_low_by_dir(CS_on_data, data_path);
    make_figC3_high_low_heatmap_target(CS_on_data, data_path);
    make_figC4_value_CS_on_only(CS_on_data, data_path);     
    make_figC5_scatter_unity(CS_on_data, data_path);
    make_figC6_pos_neg_by_dir(CS_on_data, data_path);
    make_figC7_pos_group_target(CS_on_data, data_path);
    make_figC8_pos_group_saccade(CS_on_data, data_path);
    make_figC9_sorted_heatmap_target(CS_on_data, data_path);
    make_figC10_sorted_heatmap_saccade(CS_on_data, data_path);
end

% ===================  COMPUTE PHASE  =====================

function C = compute_CS_ON(data_path)
    JDM_params_funcs
    S = load(fullfile(data_path,'population_data','cs_on_rate_spe_rpe_prim_sac'));
    % Expect: S.cs_on_data, S.cs_rate, S.cs_vm, params
    cs_on_data = S.cs_on_data;
    cs_rate    = S.cs_rate;
    cs_vm      = S.cs_vm;

    % ---- constants / config in one place
    C.cfg.leg            = {'cs-on','cs±45','cs±90','cs±135','cs±180'};
    C.cfg.colors_h       = cool(5);
    tmp                  = winter(5);
    C.cfg.colors_l       = tmp(end:-1:1,:);
    C.cfg.vm_color       = [204,174,98]/256;
    C.cfg.vm_lim         = [0 600];
    C.cfg.dirs           = 1:5;         % cs-on and symmetric pairs
    C.cfg.dir_on         = 5;           % by construction after rotation
    C.cfg.t              = -249:250;    % 500 samples
    C.cfg.t_ind_vis      = (-50:200) + 250;
    C.cfg.t_ind_sac      = (-100:100) + 250;
    C.cfg.win_vis_ms     = [50 150];
    C.cfg.win_sac_ms     = [-50 50];
    C.cfg.high_conds     = [1 2];   % forced high, etc.
    C.cfg.low_conds      = [3 4 ];
    C.cfg.rwb_map        = params.rwb_map;
    C.cfg.print_size     = [10 4];

    % ---- basic counts
    C.num_cells = numel(cs_on_data);

    % ---- pick vis vs sac for CS-on direction & build rotation order
    cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_data,'UniformOutput',false));
    cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_data,'UniformOutput',false));
    cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, cs_on_data,'UniformOutput',false));
    cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, cs_on_data,'UniformOutput',false));

    ind_m          = cs_on_rho_sac > cs_on_rho_vis;      % choose stronger vector
    cs_on_ang_rad  = cs_on_ang_vis/180*pi;               % default vis
    cs_on_ang_rad(ind_m) = cs_on_ang_sac(ind_m)/180*pi;

    ang_bins = params.sac.ang_values/180*pi;             % 8 bins
    [~, cs_on_bin] = min(abs(angdiff(repmat(ang_bins, C.num_cells,1), cs_on_ang_rad*ones(size(ang_bins)))),[],2);

    order_ang = nan(C.num_cells,8);
    for i = 1:C.num_cells
        order_ang(i,:) = circshift(1:8, 5 - cs_on_bin(i)); % rotate so cs-on becomes index 5
    end
    C.order_ang = [order_ang, order_ang(:,1)]; % 9th column repeats for wrap

    % ---- rotate rate & velocity to CS-on frame
    C.cs_rate_rot = nan(C.num_cells, 500, 8, 6, 2);
    C.cs_vm_rot   = nan(C.num_cells, 500, 8, 6, 2);
    for i = 1:C.num_cells
        C.cs_rate_rot(i,:,:,:,:) = cs_rate(i,:,C.order_ang(i,1:8),:,1:2); % pick visual & sac in last dim as-is
        C.cs_vm_rot(i,:,:,:,:)   = cs_vm(i,:,  C.order_ang(i,1:8),:,1:2);
    end

    % ---- average high/low across conditions
    C.cs_rate_h = squeeze(mean(C.cs_rate_rot(:,:,:, C.cfg.high_conds, :), 4, 'omitnan')); % [cells x time x dir x epoch]
    C.cs_rate_l = squeeze(mean(C.cs_rate_rot(:,:,:, C.cfg.low_conds,  :), 4, 'omitnan'));
    C.cs_vm_h   = squeeze(mean(C.cs_vm_rot(:,:,:,   C.cfg.high_conds, :), 4, 'omitnan'));
    C.cs_vm_l   = squeeze(mean(C.cs_vm_rot(:,:,:,   C.cfg.low_conds,  :), 4, 'omitnan'));

    % ---- reward classification (using both epochs windows)
    t = C.cfg.t;
    t_vis_idx = find(t >= C.cfg.win_vis_ms(1) & t <= C.cfg.win_vis_ms(2));
    t_sac_idx = find(t >= C.cfg.win_sac_ms(1) & t <= C.cfg.win_sac_ms(2));

    cs_h_vis = squeeze(C.cs_rate_h(:, t_vis_idx, C.cfg.dir_on, 1));
    cs_l_vis = squeeze(C.cs_rate_l(:, t_vis_idx, C.cfg.dir_on, 1));
    cs_h_sac = squeeze(C.cs_rate_h(:, t_sac_idx, C.cfg.dir_on, 2));
    cs_l_sac = squeeze(C.cs_rate_l(:, t_sac_idx, C.cfg.dir_on, 2));

    h_vis_mean = mean(cs_h_vis, 2, 'omitnan');
    l_vis_mean = mean(cs_l_vis, 2, 'omitnan');
    h_sac_mean = mean(cs_h_sac, 2, 'omitnan');
    l_sac_mean = mean(cs_l_sac, 2, 'omitnan');

    C.cs_high_tot = mean([h_vis_mean, h_sac_mean], 2, 'omitnan');
    C.cs_low_tot  = mean([l_vis_mean, l_sac_mean], 2, 'omitnan');

    C.reward_diff = C.cs_high_tot - C.cs_low_tot;
    C.idx_pos     = C.reward_diff > 0;
    C.idx_neg     = C.reward_diff < 0;

    % rho summaries
    C.cs_on_rho_sac_mean_pos = mean(cs_on_rho_sac(C.idx_pos));
    C.cs_on_rho_sac_mean_neg = mean(cs_on_rho_sac(C.idx_neg));

    % ---- time axes used in plots
    C.time_vis = C.cfg.t_ind_vis - 250;
    C.time_sac = C.cfg.t_ind_sac - 250;

   % ---- VALUE (H−L) in CS-on, both epochs (safe onsets)
    val_vis = squeeze( ...
        C.cs_rate_h(:, C.cfg.t_ind_vis, C.cfg.dir_on, 1) - ...
        C.cs_rate_l(:, C.cfg.t_ind_vis, C.cfg.dir_on, 1)); % [cells x time]
    C.value_vis_mean = mean(val_vis, 1, 'omitnan');
    C.value_vis_sem  = std(val_vis, [], 1, 'omitnan') / sqrt(C.num_cells);
    iv = find(C.value_vis_mean > 3.3 * C.value_vis_sem, 1, 'first');
    if isempty(iv)
        C.value_onset_vis_ms = NaN;
    else
        C.value_onset_vis_ms = C.time_vis(iv);
    end
    
    val_sac = squeeze( ...
        C.cs_rate_h(:, C.cfg.t_ind_sac, C.cfg.dir_on, 2) - ...
        C.cs_rate_l(:, C.cfg.t_ind_sac, C.cfg.dir_on, 2));
    C.value_sac_mean = mean(val_sac, 1, 'omitnan');
    C.value_sac_sem  = std(val_sac, [], 1, 'omitnan') / sqrt(C.num_cells);
    is_ = find(C.value_sac_mean > 3.3 * C.value_sac_sem, 1, 'first');
    if isempty(is_)
        C.value_onset_sac_ms = NaN;
    else
        C.value_onset_sac_ms = C.time_sac(is_);
    end

       % ---- VALUE (H−L) in CS-on+180, both epochs
    val_vis_180 = squeeze( ...
        C.cs_rate_h(:, C.cfg.t_ind_vis, 1, 1) - ...   % dir=1 is CS-on+180
        C.cs_rate_l(:, C.cfg.t_ind_vis, 1, 1));
    C.value_vis180_mean = mean(val_vis_180, 1, 'omitnan');
    C.value_vis180_sem  = std(val_vis_180, [], 1, 'omitnan') / sqrt(C.num_cells);
    iv180 = find(C.value_vis180_mean > 3.3 * C.value_vis180_sem, 1, 'first');
    if isempty(iv180)
        C.value_onset_vis180_ms = NaN;
    else
        C.value_onset_vis180_ms = C.time_vis(iv180);
    end

    val_sac_180 = squeeze( ...
        C.cs_rate_h(:, C.cfg.t_ind_sac, 1, 2) - ...
        C.cs_rate_l(:, C.cfg.t_ind_sac, 1, 2));
    C.value_sac180_mean = mean(val_sac_180, 1, 'omitnan');
    C.value_sac180_sem  = std(val_sac_180, [], 1, 'omitnan') / sqrt(C.num_cells);
    is180 = find(C.value_sac180_mean > 3.3 * C.value_sac180_sem, 1, 'first');
    if isempty(is180)
        C.value_onset_sac180_ms = NaN;
    else
        C.value_onset_sac180_ms = C.time_sac(is180);
    end

    % ---- DIRECTIONALITY (CS-on − CS+180), averaged across reward (safe onsets)
    cs_comb = (C.cs_rate_h + C.cs_rate_l) / 2;  % avg across reward
    
    dir_vis = squeeze( ...
        cs_comb(:, C.cfg.t_ind_vis, C.cfg.dir_on, 1) - ...
        cs_comb(:, C.cfg.t_ind_vis, 1,            1));
    C.dir_vis_mean = mean(dir_vis, 1, 'omitnan');
    C.dir_vis_sem  = std(dir_vis, [], 1, 'omitnan') / sqrt(C.num_cells);
    iv2 = find(C.dir_vis_mean > 3.3 * C.dir_vis_sem, 1, 'first');
    if isempty(iv2)
        C.dir_onset_vis_ms = NaN;
    else
        C.dir_onset_vis_ms = C.time_vis(iv2);
    end
    
    dir_sac = squeeze( ...
        cs_comb(:, C.cfg.t_ind_sac, C.cfg.dir_on, 2) - ...
        cs_comb(:, C.cfg.t_ind_sac, 1,            2));
    C.dir_sac_mean = mean(dir_sac, 1, 'omitnan');
    C.dir_sac_sem  = std(dir_sac, [], 1, 'omitnan') / sqrt(C.num_cells);
    is2 = find(C.dir_sac_mean > 3.3 * C.dir_sac_sem, 1, 'first');
    if isempty(is2)
        C.dir_onset_sac_ms = NaN;
    else
        C.dir_onset_sac_ms = C.time_sac(is2);
    end

    % ---- heatmap sorting indices (target & saccade)
    % target sorting by (low - high) in 50-150ms
    tv = find(C.time_vis >= C.cfg.win_vis_ms(1) & C.time_vis <= C.cfg.win_vis_ms(2));
    cs_high_vis = squeeze(C.cs_rate_h(:, C.cfg.t_ind_vis, C.cfg.dir_on, 1));
    cs_low_vis  = squeeze(C.cs_rate_l(:, C.cfg.t_ind_vis, C.cfg.dir_on, 1));
    diff_sort   = mean(cs_low_vis(:, tv), 2, 'omitnan') - mean(cs_high_vis(:, tv), 2, 'omitnan');
    [~, C.sort_idx_target] = sort(diff_sort, 'descend');

    % saccade sorting by (low - high) in -50..50ms
    ts = find(C.time_sac >= C.cfg.win_sac_ms(1) & C.time_sac <= C.cfg.win_sac_ms(2));
    cs_high_sac = squeeze(C.cs_rate_h(:, C.cfg.t_ind_sac, C.cfg.dir_on, 2));
    cs_low_sac  = squeeze(C.cs_rate_l(:, C.cfg.t_ind_sac, C.cfg.dir_on, 2));
    diff_sort_s = mean(cs_low_sac(:, ts), 2, 'omitnan') - mean(cs_high_sac(:, ts), 2, 'omitnan');
    [~, C.sort_idx_saccade] = sort(diff_sort_s, 'descend');
end

% ===================  PLOT PHASE  ========================

function make_figC1_all_dirs(C, data_path)
    num_cs = C.num_cells; dirs = C.cfg.dirs;
    fig = figure;
    colors_t = cool(5);

    % pre-avg across conds already in C.cs_rate_h/l; here do across dirs per loop
    vm_vis_all = squeeze(mean(C.cs_vm_h(:, C.cfg.t_ind_vis, :, 1), 3, 'omitnan'));
    vm_sac_all = squeeze(mean(C.cs_vm_h(:, C.cfg.t_ind_sac, :, 2), 3, 'omitnan'));

    for k = numel(dirs):-1:1
        d = dirs(k);
        if d==1 || d==5
            dset = d;
        else
            dset = [d, 6-d];
        end
        % VIS
        cs_vis = squeeze(mean(C.cs_rate_h(:, C.cfg.t_ind_vis, dset, 1), 3, 'omitnan'));
        subplot(1,2,1); title(sprintf('CS Cells (n = %d)', num_cs))
        colororder({'k','k'}); yyaxis left
        [hl,hp] = boundedline(C.time_vis, mean(cs_vis,1,'omitnan'), std(cs_vis,[],1,'omitnan')/sqrt(num_cs), ...
                              'Color', colors_t(k,:), 'alpha');
        hl.DisplayName = C.cfg.leg{6-k}; hp.HandleVisibility = 'off';
        hold on; xlabel('time from target onset (ms)'); xline(0,'--','HandleVisibility','off'); legend('Box','off')
        if k == 1
            yyaxis right; area(C.time_vis, mean(vm_vis_all,1,'omitnan'), ...
                'FaceColor', C.cfg.vm_color, 'FaceAlpha', .3, 'EdgeColor','none','HandleVisibility','off');
            ylim(C.cfg.vm_lim);
        end
        % SAC
        cs_sac = squeeze(mean(C.cs_rate_h(:, C.cfg.t_ind_sac, dset, 2), 3, 'omitnan'));
        subplot(1,2,2); title(sprintf('CS Cells (n = %d)', num_cs))
        colororder({'k','k'}); yyaxis left
        [hl,hp] = boundedline(C.time_sac, mean(cs_sac,1,'omitnan'), std(cs_sac,[],1,'omitnan')/sqrt(num_cs), ...
                              'Color', colors_t(k,:), 'alpha');
        hl.DisplayName = C.cfg.leg{6-k}; hp.HandleVisibility = 'off';
        hold on; xlabel('time from saccade onset (ms)'); xline(0,'--','HandleVisibility','off');
        if k == 1
            yyaxis right; area(C.time_sac, mean(vm_sac_all,1,'omitnan'), ...
                'FaceColor', C.cfg.vm_color, 'FaceAlpha', .3, 'EdgeColor','none','HandleVisibility','off');
            ylim(C.cfg.vm_lim); ylabel('Velocity (deg/s)');
        end
    end
    ESN_Beautify_Plot(fig, C.cfg.print_size);
    print(fig, fullfile(data_path,'population_figs','ephys','Fig1-C1.pdf'), '-dpdf', '-bestfit');
end

function make_figC2_high_low_by_dir(C, data_path)
    num_cs = C.num_cells; dirs = C.cfg.dirs;
    fig = figure;
    for k = numel(dirs):-1:1
        d = dirs(k);
        if d==1 || d==5
            dset = d;
        else
            dset = [d, 6-d];
        end

        % VIS
        cs_h_vis = squeeze(mean(C.cs_rate_h(:, C.cfg.t_ind_vis, dset, 1), 3));
        cs_l_vis = squeeze(mean(C.cs_rate_l(:, C.cfg.t_ind_vis, dset, 1), 3));
        if k==1
            vm_vis_h = squeeze(mean(C.cs_vm_h(:, C.cfg.t_ind_vis, :, 1), 3));
        end
        subplot(1,2,1); title(sprintf('CS Cells (n = %d)', num_cs))
        colororder({'k','k'}); yyaxis left
        [h1,p1] = boundedline(C.time_vis, mean(cs_h_vis,1,'omitnan'), std(cs_h_vis,[],1,'omitnan')/sqrt(num_cs), ...
                              'Color', C.cfg.colors_h(k,:), 'alpha');
        [h2,p2] = boundedline(C.time_vis, mean(cs_l_vis,1,'omitnan'), std(cs_l_vis,[],1,'omitnan')/sqrt(num_cs), ...
                              'Color', C.cfg.colors_h(k,:), 'alpha');
        set(h2,'LineStyle','--');
        h1.DisplayName = ['High - ' C.cfg.leg{6-k}];
        h2.DisplayName = ['Low - '  C.cfg.leg{6-k}];
        p1.HandleVisibility='off'; p2.HandleVisibility='off';
        hold on; xlabel('time from stimulus onset (ms)'); xline(0,'--','HandleVisibility','off');
        if k==1
            yyaxis right; area(C.time_vis, mean(vm_vis_h,1,'omitnan'), 'FaceColor',C.cfg.vm_color,'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off'); ylim(C.cfg.vm_lim);
        end

        % SAC
        cs_h_sac = squeeze(mean(C.cs_rate_h(:, C.cfg.t_ind_sac, dset, 2), 3));
        cs_l_sac = squeeze(mean(C.cs_rate_l(:, C.cfg.t_ind_sac, dset, 2), 3));
        if k==1
            vm_sac_h = squeeze(mean(C.cs_vm_h(:, C.cfg.t_ind_sac, :, 2), 3));
        end
        subplot(1,2,2); title(sprintf('CS Cells (n = %d)', num_cs))
        colororder({'k','k'}); yyaxis left
        [h3,p3] = boundedline(C.time_sac, mean(cs_h_sac,1,'omitnan'), std(cs_h_sac,[],1,'omitnan')/sqrt(num_cs), ...
                              'Color', C.cfg.colors_h(k,:), 'alpha');
        [h4,p4] = boundedline(C.time_sac, mean(cs_l_sac,1,'omitnan'), std(cs_l_sac,[],1,'omitnan')/sqrt(num_cs), ...
                              'Color', C.cfg.colors_h(k,:), 'alpha');
        set(h4,'LineStyle','--');
        h3.DisplayName = ['High - ' C.cfg.leg{6-k}];
        h4.DisplayName = ['Low - '  C.cfg.leg{6-k}];
        p3.HandleVisibility='off'; p4.HandleVisibility='off';
        hold on; xlabel('time from saccade onset (ms)'); xline(0,'--','HandleVisibility','off');
        if k==1
            yyaxis right; area(C.time_sac, mean(vm_sac_h,1,'omitnan'), 'FaceColor',C.cfg.vm_color,'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off'); ylim(C.cfg.vm_lim); ylabel('Velocity (deg/s)');
        end
        legend('Box','off','Location','northeast');
    end
    ESN_Beautify_Plot(fig, C.cfg.print_size);
    print(fig, fullfile(data_path,'population_figs','ephys','Fig1-C2.pdf'), '-dpdf', '-bestfit');
end

function make_figC3_high_low_heatmap_target(C, data_path)
    num_cs = C.num_cells;
    % Build z-scored matrices (target-aligned, CS-on)
    cs_h = squeeze(C.cs_rate_h(:, C.cfg.t_ind_vis, C.cfg.dir_on, 1));
    cs_l = squeeze(C.cs_rate_l(:, C.cfg.t_ind_vis, C.cfg.dir_on, 1));
    cs_h_z = (cs_h - mean(cs_h,2,'omitnan')) ./ std(cs_h,0,2,'omitnan');
    cs_l_z = (cs_l - mean(cs_l,2,'omitnan')) ./ std(cs_l,0,2,'omitnan');

    % Random ordering (kept from your original)
    rng(42); idx = randperm(num_cs);
    H = cs_h_z(idx,:); L = cs_l_z(idx,:);

    fig = figure;
    subplot(1,2,1); imagesc(C.time_vis, 1:num_cs, H);
    xlabel('Time from target onset (ms)'); ylabel('Neurons (sorted by peak time)'); title('High Value — CS-on'); colormap(C.cfg.rwb_map); caxis([0 4]); colorbar; xline(0,'--w','LineWidth',1.2);
    subplot(1,2,2); imagesc(C.time_vis, 1:num_cs, L);
    xlabel('Time from target onset (ms)'); ylabel('Neurons (same order)'); title('Low Value — CS-on');  colormap(C.cfg.rwb_map); caxis([0 4]); colorbar; xline(0,'--w','LineWidth',1.2);
    ESN_Beautify_Plot(fig, C.cfg.print_size);
    print(fig, fullfile(data_path,'population_figs','ephys','Fig1-C3.pdf'), '-dpdf', '-bestfit');
end

function make_figC4_value_CS_on_only(C, data_path)
    fig = figure; 
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

    % ===== Left: target-aligned =====
    nexttile; hold on
    % Value (H-L, CS-on)
    [hv,pv] = boundedline(C.time_vis, C.value_vis_mean, C.value_vis_sem, 'r', 'alpha');
    pv.HandleVisibility = 'off'; hv.DisplayName = 'Value (H−L, \theta)';

    % Value (H-L, CS-on+180)
    [hv180,pv180] = boundedline(C.time_vis, C.value_vis180_mean, C.value_vis180_sem, 'm', 'alpha');
    pv180.HandleVisibility = 'off'; hv180.DisplayName = 'Value (H−L, \theta+180)';

    % Direction (CS-on − CS+180, avg reward)
    [hd,pd] = boundedline(C.time_vis, C.dir_vis_mean, C.dir_vis_sem, 'b', 'alpha');
    pd.HandleVisibility = 'off'; hd.DisplayName = 'Direction (\theta−\theta+180)';

    xline(0,'--k','HandleVisibility','off')
    if ~isnan(C.value_onset_vis_ms),    xline(C.value_onset_vis_ms,'--r','HandleVisibility','off'); end
    if ~isnan(C.dir_onset_vis_ms),      xline(C.dir_onset_vis_ms,'--b','HandleVisibility','off');   end

    xlabel('Time from target onset (ms)');
    ylabel('CS rate difference (Hz)');
    title('Value & Direction — target-aligned');
    legend([hv hv180 hd], ...
        {'Value (H−L, \theta)','Value (H−L, \theta+180)','Direction (\theta−\theta+180)'}, ...
        'Box','off','Location','northeast')

    % ===== Right: saccade-aligned =====
    nexttile; hold on
    % Value (H-L, CS-on)
    [hv2,pv2] = boundedline(C.time_sac, C.value_sac_mean, C.value_sac_sem, 'r', 'alpha');
    pv2.HandleVisibility = 'off'; hv2.DisplayName = 'Value (H−L, \theta)';

    % Value (H-L, CS-on+180)
    [hv2180,pv2180] = boundedline(C.time_sac, C.value_sac180_mean, C.value_sac180_sem, 'm', 'alpha');
    pv2180.HandleVisibility = 'off'; hv2180.DisplayName = 'Value (H−L, \theta+180)';

    % Direction (CS-on − CS+180, avg reward)
    [hd2,pd2] = boundedline(C.time_sac, C.dir_sac_mean, C.dir_sac_sem, 'b', 'alpha');
    pd2.HandleVisibility = 'off'; hd2.DisplayName = 'Direction (\theta−\theta+180)';

    xline(0,'--k','HandleVisibility','off')
    %if ~isnan(C.value_onset_sac_ms),    xline(C.value_onset_sac_ms,'--r','HandleVisibility','off'); end
    if ~isnan(C.dir_onset_sac_ms),      xline(C.dir_onset_sac_ms,'--b','HandleVisibility','off');   end

    xlabel('Time from saccade onset (ms)');
    ylabel('CS rate difference (Hz)');
    title('Value & Direction — saccade-aligned');
    legend([hv2 hv2180 hd2], ...
        {'Value (H−L, \theta)','Value (H−L, \theta+180)','Direction (\theta−\theta+180)'}, ...
        'Box','off','Location','northeast')

    ESN_Beautify_Plot(fig, C.cfg.print_size);
    print(fig, fullfile(data_path,'population_figs','ephys','Fig1-E.pdf'), '-dpdf', '-bestfit');

    % console reporting
    fprintf('Value onset (target, \theta):      %s ms\n', num2str(C.value_onset_vis_ms));
    fprintf('Direction onset (target):         %s ms\n', num2str(C.dir_onset_vis_ms));
    fprintf('Value onset (saccade, \theta):     %s ms\n', num2str(C.value_onset_sac_ms));
    fprintf('Direction onset (saccade):        %s ms\n', num2str(C.dir_onset_sac_ms));
end



function make_figC5_scatter_unity(C, data_path)
    num_cs = C.num_cells;
    above_unity = sum(C.cs_high_tot > C.cs_low_tot);
    max_val = max([C.cs_low_tot; C.cs_high_tot],[],'omitnan');

    fig = figure;
    subplot(1,2,1)
    scatter(C.cs_low_tot, C.cs_high_tot, 25, 'filled', 'MarkerFaceAlpha', .6, 'MarkerEdgeAlpha', .2); hold on
    plot([0 max_val],[0 max_val],'k--','LineWidth',1.5)
    xlabel('Low reward CS rate'); ylabel('High reward CS rate');
    title(sprintf('Reward Modulation of CS (N = %d, above unity = %d)', num_cs, above_unity));
    axis equal; xlim([0 max_val]); ylim([0 max_val]);

    subplot(1,2,2)
    scatter(C.cs_low_tot(C.idx_pos), C.cs_high_tot(C.idx_pos), 30, [0.2 0.6 1], 'filled', 'MarkerFaceAlpha', .6); hold on
    scatter(C.cs_low_tot(C.idx_neg), C.cs_high_tot(C.idx_neg), 30, [1 0.3 0.3], 'filled', 'MarkerFaceAlpha', .6);
    plot([0 max_val],[0 max_val],'k--','LineWidth',1.5)
    xlabel('Low reward CS rate'); ylabel('High reward CS rate');
    legend({'Reward-Positive','Reward-Negative'}, 'Location','southeast','Box','off');
    title(sprintf('Reward Modulation (Pos = %d, Neg = %d)', sum(C.idx_pos), sum(C.idx_neg)));
    axis equal; xlim([0 max_val]); ylim([0 max_val]);

    ESN_Beautify_Plot(fig, C.cfg.print_size);
    print(fig, fullfile(data_path,'population_figs','ephys','Fig1-C5.pdf'), '-dpdf', '-bestfit');
end

function make_figC6_pos_neg_by_dir(C, data_path)
    % Compares reward-positive vs reward-negative across directions (both epochs)
    dirs = C.cfg.dirs;
    n_pos = sum(C.idx_pos); n_neg = sum(C.idx_neg);
    fig = figure;

    for k = numel(dirs):-1:1
        if d==1 || d==5
            dset = d;
        else
            dset = [d, 6-d];
        end

        % POS group
        pos_vis = squeeze(mean(C.cs_rate_rot(C.idx_pos, C.cfg.t_ind_vis, dset, 1, 1), 3));
        pos_sac = squeeze(mean(C.cs_rate_rot(C.idx_pos, C.cfg.t_ind_sac, dset, 2, 1), 3));
        if k==1
            pos_vm_vis = squeeze(mean(C.cs_vm_rot(C.idx_pos, C.cfg.t_ind_vis, :, 1, 1), 3));
            pos_vm_sac = squeeze(mean(C.cs_vm_rot(C.idx_pos, C.cfg.t_ind_sac, :, 2, 1), 3));
        end
        % NEG group
        neg_vis = squeeze(mean(C.cs_rate_rot(C.idx_neg, C.cfg.t_ind_vis, dset, 1, 1), 3));
        neg_sac = squeeze(mean(C.cs_rate_rot(C.idx_neg, C.cfg.t_ind_sac, dset, 2, 1), 3));

        % VIS panel
        subplot(1,2,1); title(sprintf('rho pos=%.2f, neg=%.2f', C.cs_on_rho_sac_mean_pos, C.cs_on_rho_sac_mean_neg),'FontSize',8);
        colororder({'k','k'}); yyaxis left
        [h1,p1] = boundedline(C.time_vis, mean(pos_vis,1,'omitnan'), std(pos_vis,[],1,'omitnan')/sqrt(n_pos), 'Color', C.cfg.colors_h(k,:), 'alpha');
        [h2,p2] = boundedline(C.time_vis, mean(neg_vis,1,'omitnan'), std(neg_vis,[],1,'omitnan')/sqrt(n_neg), 'alpha');
        set(h2,'LineStyle','--','Color', C.cfg.colors_l(k,:));
        h1.DisplayName = ['pos - ' C.cfg.leg{6-k}]; h2.DisplayName = ['neg - ' C.cfg.leg{6-k}];
        p1.HandleVisibility='off'; p2.HandleVisibility='off';
        hold on; xlabel('time from stimulus onset (ms)'); xline(0,'--','HandleVisibility','off');
        if k==1
            yyaxis right; area(C.time_vis, mean(pos_vm_vis,1,'omitnan'), 'FaceColor', C.cfg.vm_color, 'FaceAlpha', .3, 'EdgeColor','none','HandleVisibility','off'); ylim(C.cfg.vm_lim);
        end

        % SAC panel
        subplot(1,2,2);
        colororder({'k','k'}); yyaxis left
        [h3,p3] = boundedline(C.time_sac, mean(pos_sac,1,'omitnan'), std(pos_sac,[],1,'omitnan')/sqrt(n_pos), 'Color', C.cfg.colors_h(k,:), 'alpha');
        [h4,p4] = boundedline(C.time_sac, mean(neg_sac,1,'omitnan'), std(neg_sac,[],1,'omitnan')/sqrt(n_neg), 'alpha');
        set(h4,'LineStyle','--','Color', C.cfg.colors_l(k,:));
        h3.DisplayName = ['pos - ' C.cfg.leg{6-k}]; h4.DisplayName = ['neg - ' C.cfg.leg{6-k}];
        p3.HandleVisibility='off'; p4.HandleVisibility='off';
        hold on; xlabel('time from saccade onset (ms)'); xline(0,'--','HandleVisibility','off');
        if k==1
            yyaxis right; area(C.time_sac, mean(pos_vm_sac,1,'omitnan'), 'FaceColor', C.cfg.vm_color, 'FaceAlpha', .3, 'EdgeColor','none','HandleVisibility','off'); ylim(C.cfg.vm_lim); ylabel('Velocity (deg/s)');
        end
        legend('Box','off','Location','northeast');
    end
    ESN_Beautify_Plot(fig, C.cfg.print_size);
    print(fig, fullfile(data_path,'population_figs','ephys','Fig1-C6.pdf'), '-dpdf', '-bestfit');
end

function make_figC7_pos_group_target(C, data_path)
    % Reward-positive vs negative (target-aligned, CS-on)
    pos_h = squeeze(C.cs_rate_h(C.idx_pos, C.cfg.t_ind_vis, C.cfg.dir_on, 1));
    pos_l = squeeze(C.cs_rate_l(C.idx_pos, C.cfg.t_ind_vis, C.cfg.dir_on, 1));
    neg_h = squeeze(C.cs_rate_h(C.idx_neg, C.cfg.t_ind_vis, C.cfg.dir_on, 1));
    neg_l = squeeze(C.cs_rate_l(C.idx_neg, C.cfg.t_ind_vis, C.cfg.dir_on, 1));

    pos_vm = squeeze(C.cs_vm_h(C.idx_pos, C.cfg.t_ind_vis, C.cfg.dir_on, 1));
    neg_vm = squeeze(C.cs_vm_h(C.idx_neg, C.cfg.t_ind_vis, C.cfg.dir_on, 1));

    fig = figure;
    subplot(1,2,1); hold on; colororder({'k','k'}); yyaxis left
    [h1,p1] = boundedline(C.time_vis, mean(pos_h,1,'omitnan'), std(pos_h,[],1,'omitnan')/sqrt(size(pos_h,1)), 'm', 'alpha');
    [h2,p2] = boundedline(C.time_vis, mean(pos_l,1,'omitnan'), std(pos_l,[],1,'omitnan')/sqrt(size(pos_l,1)), 'c--', 'alpha');
    p1.HandleVisibility='off'; p2.HandleVisibility='off'; ylabel('CS rate (Hz)');
    yyaxis right; area(C.time_vis, mean(pos_vm,1,'omitnan'), 'FaceColor',C.cfg.vm_color,'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off'); ylim(C.cfg.vm_lim);
    xline(0,'--k'); xlabel('Time from target onset (ms)'); title(sprintf('Reward-Positive Cells (n=%d)', size(pos_h,1)));
    legend([h1 h2], {'High','Low'}, 'Location','northeast','Box','off');

    subplot(1,2,2); hold on; colororder({'k','k'}); yyaxis left
    [h3,p3] = boundedline(C.time_vis, mean(neg_h,1,'omitnan'), std(neg_h,[],1,'omitnan')/sqrt(size(neg_h,1)), 'm', 'alpha');
    [h4,p4] = boundedline(C.time_vis, mean(neg_l,1,'omitnan'), std(neg_l,[],1,'omitnan')/sqrt(size(neg_l,1)), 'c--', 'alpha');
    p3.HandleVisibility='off'; p4.HandleVisibility='off';
    yyaxis right; area(C.time_vis, mean(neg_vm,1,'omitnan'), 'FaceColor',C.cfg.vm_color,'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off'); ylim(C.cfg.vm_lim); ylabel('velocity (deg/s)');
    xline(0,'--k'); xlabel('Time from target onset (ms)'); title(sprintf('Reward-Negative Cells (n=%d)', size(neg_h,1)));
    legend([h3 h4], {'High','Low'}, 'Location','northeast','Box','off');

    ESN_Beautify_Plot(fig, C.cfg.print_size);
    print(fig, fullfile(data_path,'population_figs','ephys','Fig1-C7.pdf'), '-dpdf', '-bestfit');
end

function make_figC8_pos_group_saccade(C, data_path)
    % Reward-positive vs negative (saccade-aligned, CS-on)
    pos_h = squeeze(C.cs_rate_h(C.idx_pos, C.cfg.t_ind_sac, C.cfg.dir_on, 2));
    pos_l = squeeze(C.cs_rate_l(C.idx_pos, C.cfg.t_ind_sac, C.cfg.dir_on, 2));
    neg_h = squeeze(C.cs_rate_h(C.idx_neg, C.cfg.t_ind_sac, C.cfg.dir_on, 2));
    neg_l = squeeze(C.cs_rate_l(C.idx_neg, C.cfg.t_ind_sac, C.cfg.dir_on, 2));

    pos_vm = squeeze(C.cs_vm_h(C.idx_pos, C.cfg.t_ind_sac, C.cfg.dir_on, 2));
    neg_vm = squeeze(C.cs_vm_h(C.idx_neg, C.cfg.t_ind_sac, C.cfg.dir_on, 2));

    fig = figure;
    subplot(1,2,1); hold on; colororder({'k','k'}); yyaxis left
    [h1,p1] = boundedline(C.time_sac, mean(pos_h,1,'omitnan'), std(pos_h,[],1,'omitnan')/sqrt(size(pos_h,1)), 'm', 'alpha');
    [h2,p2] = boundedline(C.time_sac, mean(pos_l,1,'omitnan'), std(pos_l,[],1,'omitnan')/sqrt(size(pos_l,1)), 'c--', 'alpha');
    p1.HandleVisibility='off'; p2.HandleVisibility='off'; ylabel('CS rate (Hz)');
    yyaxis right; area(C.time_sac, mean(pos_vm,1,'omitnan'), 'FaceColor',C.cfg.vm_color,'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off'); ylim(C.cfg.vm_lim);
    xline(0,'--k'); xlabel('Time from saccade onset (ms)'); title(sprintf('Reward-Positive Cells (n=%d)', size(pos_h,1)));
    legend([h1 h2], {'High','Low'}, 'Location','northeast','Box','off');

    subplot(1,2,2); hold on; colororder({'k','k'}); yyaxis left
    [h3,p3] = boundedline(C.time_sac, mean(neg_h,1,'omitnan'), std(neg_h,[],1,'omitnan')/sqrt(size(neg_h,1)), 'm', 'alpha');
    [h4,p4] = boundedline(C.time_sac, mean(neg_l,1,'omitnan'), std(neg_l,[],1,'omitnan')/sqrt(size(neg_l,1)), 'c--', 'alpha');
    p3.HandleVisibility='off'; p4.HandleVisibility='off';
    yyaxis right; area(C.time_sac, mean(neg_vm,1,'omitnan'), 'FaceColor',C.cfg.vm_color,'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off'); ylim(C.cfg.vm_lim); ylabel('velocity (deg/s)');
    xline(0,'--k'); xlabel('Time from saccade onset (ms)'); title(sprintf('Reward-Negative Cells (n=%d)', size(neg_h,1)));
    legend([h3 h4], {'High','Low'}, 'Location','northeast','Box','off');

    ESN_Beautify_Plot(fig, C.cfg.print_size);
    print(fig, fullfile(data_path,'population_figs','ephys','Fig1-C8.pdf'), '-dpdf', '-bestfit');
end

function make_figC9_sorted_heatmap_target(C, data_path)
    cs_h = squeeze(C.cs_rate_h(:, C.cfg.t_ind_vis, C.cfg.dir_on, 1));
    cs_l = squeeze(C.cs_rate_l(:, C.cfg.t_ind_vis, C.cfg.dir_on, 1));
    cs_h_z = (cs_h - mean(cs_h,2,'omitnan'))./std(cs_h,0,2,'omitnan');
    cs_l_z = (cs_l - mean(cs_l,2,'omitnan'))./std(cs_l,0,2,'omitnan');

    H = cs_h_z(C.sort_idx_target, :);
    L = cs_l_z(C.sort_idx_target, :);

    fig = figure;
    subplot(1,2,1); imagesc(C.time_vis, 1:size(L,1), L);
    xlabel('Time from target onset (ms)'); ylabel('Neurons (sorted by reward modulation)'); title('Low Reward — CS-on'); xline(0,'--w','LineWidth',1.2); colormap(C.cfg.rwb_map); caxis([0 4]); colorbar;
    subplot(1,2,2); imagesc(C.time_vis, 1:size(H,1), H);
    xlabel('Time from target onset (ms)'); ylabel('Neurons (same order)'); title('High Reward — CS-on');       xline(0,'--w','LineWidth',1.2); colormap(C.cfg.rwb_map); caxis([0 4]); colorbar;

    ESN_Beautify_Plot(fig, C.cfg.print_size);
    print(fig, fullfile(data_path,'population_figs','ephys','Fig1-C9.pdf'), '-dpdf', '-bestfit');
end

function make_figC10_sorted_heatmap_saccade(C, data_path)
    cs_h = squeeze(C.cs_rate_h(:, C.cfg.t_ind_sac, C.cfg.dir_on, 2));
    cs_l = squeeze(C.cs_rate_l(:, C.cfg.t_ind_sac, C.cfg.dir_on, 2));
    cs_h_z = (cs_h - mean(cs_h,2,'omitnan'))./std(cs_h,0,2,'omitnan');
    cs_l_z = (cs_l - mean(cs_l,2,'omitnan'))./std(cs_l,0,2,'omitnan');

    H = cs_h_z(C.sort_idx_saccade, :);
    L = cs_l_z(C.sort_idx_saccade, :);

    fig = figure;
    subplot(1,2,1); imagesc(C.time_sac, 1:size(L,1), L);
    xlabel('Time from saccade onset (ms)'); ylabel('Neurons (sorted by reward modulation)'); title('Low Reward — CS-on (saccade-aligned)'); xline(0,'--w','LineWidth',1.2); colormap(C.cfg.rwb_map); caxis([0 4]); colorbar;
    subplot(1,2,2); imagesc(C.time_sac, 1:size(H,1), H);
    xlabel('Time from saccade onset (ms)'); ylabel('Neurons (same order)'); title('High Reward — CS-on (saccade-aligned)'); xline(0,'--w','LineWidth',1.2); colormap(C.cfg.rwb_map); caxis([0 4]); colorbar;

    ESN_Beautify_Plot(fig, C.cfg.print_size);
    print(fig, fullfile(data_path,'population_figs','ephys','Fig1-C10.pdf'), '-dpdf', '-bestfit');
end
