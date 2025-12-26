%% ==================== NEW: Projected-space plots per amplitude ====================
function plot_SS_rate_projected_ampbins(data_path)
    JDM_params_funcs
    save_path = fullfile(data_path,'population_data');

    S        = load(fullfile(save_path,'SS_population_clique_combined_sac.mat'),'data','meta');
    data_ss  = S.data;
    meta     = S.meta;
    idx      = load(fullfile(save_path,'purkinje_cell_ids.mat'),'ind_p','ind_b');

    ind_b = idx.ind_b;
    ind_p = idx.ind_p;
    ind_all = ind_b | ind_p;

    % Shapes: [cell x time x dir x reward x amp]
    ss_rate  = data_ss.rate_tot_sac_amp;
    vm_trace = data_ss.vm_tot_sac_amp;

    ss_cs_rho = data_ss.cs_on_rho_tot;       % [cell x 1] weights

    time_ind  = (-100:150)+250;
    time_axis = time_ind - 250;

    high_inds = [1 2 5 6];
    low_inds  = [3 4 7 8];

    amp_edges  = meta.amp_edges(:)';
    n_amp_bins = numel(amp_edges)-1;

    colors = {[0.85 0.20 0.20], [0.20 0.40 0.90], [0 0 0]}; % burster, pauser, all
    labels = {'Bursters','Pausers','All SS'};

    for b = 1:n_amp_bins
        amp_label = sprintf('%g–%g°', amp_edges(b), amp_edges(b+1));
        fig = figure('Position',[80 80 1200 650]);
        colororder({'k','k'});
        sgtitle(sprintf('Projected SS (CS-on axis)  |  Amplitude %s', amp_label),'FontSize',14)

        % ---- HIGH ----
        subplot(1,2,1); hold on; title('High reward'); xlabel('Time (ms)'); ylabel('Projected rate (Hz)')
        plot_projected_panel(ss_rate, ss_cs_rho, {ind_b, ind_p, ind_all}, labels, colors, ...
                             time_ind, time_axis, high_inds, b);
        add_velocity_fill(vm_trace, time_ind, time_axis, high_inds, b);
        legend('Location','northeast','Box','off')

        % ---- LOW ----
        subplot(1,2,2); hold on; title('Low reward'); xlabel('Time (ms)'); ylabel('Projected rate (Hz)')
        plot_projected_panel(ss_rate, ss_cs_rho, {ind_b, ind_p, ind_all}, labels, colors, ...
                             time_ind, time_axis, low_inds, b);
        add_velocity_fill(vm_trace, time_ind, time_axis, low_inds, b);

        ESN_Beautify_Plot(fig,[10,6])
        % saveas(fig, fullfile(data_path,'population_figs','ephys', ...
        %     sprintf('SS_rate_projected_amp_%g-%g.png', amp_edges(b), amp_edges(b+1))));
    end
end

%% helper: plot the 3 projected traces (bursters/pausers/all) with ρ-weighted mean ± SEM
function plot_projected_panel(ss_rate, ss_cs_rho, masks, labels, colors, time_ind, time_axis, reward_sel, amp_bin)
    proj_w = -cosd(0:45:315);                         % 8-dir projection weights
    Wdir   = reshape(proj_w, 1, 1, 8);                % [1 x 1 x dir] for implicit expansion

    for k = 1:numel(masks)
        cell_mask = masks{k};

        % ss_rate: [cell x time x dir x reward x amp]
        tmp = ss_rate(cell_mask, time_ind, :, reward_sel, amp_bin);              % [Nc x T x D x R]
        tmp = mean(tmp, 4, 'omitnan');                                           % avg over reward -> [Nc x T x D]
        proj_mat = sum(tmp .* Wdir, 3, 'omitnan');                               % weighted sum over dir -> [Nc x T]

        w = ss_cs_rho(cell_mask);                                               % [Nc x 1] weights
        [mu, sem] = weighted_mean_sem(proj_mat, w);                              % [1 x T] each

        [hl, hp] = boundedline(time_axis, mu, sem, 'color', colors{k}, 'alpha');
        hl.DisplayName = labels{k};
        hp.HandleVisibility = 'off';
        xline(0,'--','HandleVisibility','off')
        yline(0,'--','HandleVisibility','off')
        ylim([-45,50])
    end
end

%% helper: velocity fill (averaged reward -> dir -> cells), safe 1×T vector
function add_velocity_fill(vm_trace, time_ind, time_axis, reward_sel, amp_bin)
    vm_color = [204,174,98]/256;
    vtmp = vm_trace(:, time_ind, :, reward_sel, amp_bin);   % [cells x T x dir x R]
    vtmp = mean(vtmp, 4, 'omitnan');                        % over reward
    vtmp = mean(vtmp, 3, 'omitnan');                        % over dir
    vtmp = mean(vtmp, 1, 'omitnan');                        % over cells
    sig_v = reshape(squeeze(vtmp), 1, []);
    yyaxis right
    area(time_axis, sig_v, 'FaceColor', vm_color, 'FaceAlpha', .30, 'EdgeColor','none','HandleVisibility','off');
    ylabel('Velocity (deg/s)')
    ylim([0, max(600, ceil(max(sig_v)*1.1))]);
    yyaxis left
end

%% helper: ρ-weighted mean & SEM across cells (handles NaNs)
function [mu, sem] = weighted_mean_sem(X, w)
    % X: [Nc x T], w: [Nc x 1] (nonnegative)
    w = w(:);
    Nc = size(X,1);
    T  = size(X,2);
    W  = repmat(w, 1, T);
    M  = isfinite(X);
    W  = W .* M;                                  % zero weight where X is NaN

    wsum = sum(W, 1);                              % [1 x T]
    wsum_safe = max(wsum, eps);

    Wn = W ./ wsum_safe;                           % normalized weights per time
    mu = nansum(Wn .* X, 1);                       % [1 x T]

    diff = X - repmat(mu, Nc, 1);
    varw = nansum(Wn .* (diff.^2), 1);             % weighted variance

    Neff = (wsum_safe.^2) ./ nansum(W.^2, 1 + 0*(1)); % effective N per time
    sem  = sqrt(varw) ./ sqrt(max(Neff, 1));
end
