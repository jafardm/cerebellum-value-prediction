function plot_SS_rate_pime_corr_combined(data_path)
    % === Load ===
    JDM_params_funcs
    img_save_path = fullfile(data_path,'population_figs','ephys');
    save_path     = fullfile(data_path,'population_data');

    S        = load(fullfile(save_path,'SS_population_clique_combined_sac.mat'),'data','meta');
    data_ss  = S.data;
    meta     = S.meta;                    % to get Amp_edges for labels
    idx      = load(fullfile(save_path,'purkinje_cell_ids.mat'),'ind_p','ind_b');

    ind_b = idx.ind_b;
    ind_p = idx.ind_p;

    % Shapes: [cell x time x dir x reward x amp]
    ss_rate  = data_ss.rate_tot_sac_amp;
    vm_trace = data_ss.vm_tot_sac_amp;

    % --- plotting params ---
    vm_color = [204,174,98]/256;
    time_ind  = (-100:150)+250;
    time_axis = time_ind - 250;

    high_inds = [1 2 5 6];   % forced high
    low_inds  = [3 4 7 8];   % forced low

    group_names = {'Bursters','Pausers'};

    % Amp range labels
    amp_edges  = meta.amp_edges(:)';         % e.g., [2 3 4 5 100]
    n_amp_bins = numel(amp_edges) - 1;

    for b = 1:n_amp_bins
        amp_label = sprintf('%g–%g°', amp_edges(b), amp_edges(b+1));

        fig = figure('Position',[80 80 1400 900]);
        colororder({'k','k'});
        sgtitle(sprintf('SS rate — 8 directions  |  Amplitude %s', amp_label),'FontSize',14)
        for g = 1:2
            make_dir_plots_one_ampbin( ...
                ss_rate, vm_trace, ind_b, ind_p, group_names{g}, ...
                time_axis, time_ind, high_inds, low_inds, vm_color, g, b);
        end

        ESN_Beautify_Plot(fig,[12,10])
        % Uncomment to save:
        % saveas(fig, fullfile(img_save_path, sprintf('SS_rate_8dirs_amp_%g-%g.png', amp_edges(b), amp_edges(b+1))));
    end
end


%% --- Subfunction: 8-direction plots for a given amp bin ---
function make_dir_plots_one_ampbin(ss_rate, vm_trace, ind_b, ind_p, group_name, ...
                                   time_axis, time_ind, high_inds, low_inds, vm_color, gid, amp_bin)

    colors_ = lines(8);
    leg_ = {'CS-on','+45','+90','+135','+180','-135','-90','-45'};

    % --- Decide which cell mask to use for top panels ---
    if gid == 1
        cell_mask = ind_b;
    else
        cell_mask = ind_p;
    end

    %% ---------- HIGH ----------
    subplot(4,2,(gid-1)*2+1); hold on;
    title([group_name ' — High reward'])
    xlabel('Time (ms)'); ylabel('Rate (Hz)');

    % Left y-axis: 8 directions
    yyaxis left
    for d = 1:8
        tmp = ss_rate(cell_mask, time_ind, d, high_inds, amp_bin);
        tmp = squeeze(mean(tmp, 4, 'omitnan'));          % [cells x time]
        sig = mean(tmp, 1, 'omitnan');                   
        sem = std(tmp, 0, 1, 'omitnan') ./ sqrt(sum(cell_mask));

        [h1,h2] = boundedline(time_axis, sig, sem, 'color', colors_(d,:), 'alpha');
        h1.DisplayName = leg_{d};
        h2.HandleVisibility = 'off';
    end
    xline(0,'--','HandleVisibility','off');

    % Right y-axis: velocity
    yyaxis right
    vtmp = vm_trace(cell_mask, time_ind, :, high_inds, amp_bin);
    vtmp = mean(vtmp,[3 4],'omitnan');  % average over dir+reward
    vtmp = mean(vtmp,1,'omitnan');      % average over cells
    sig_v = squeeze(vtmp);
    area(time_axis, sig_v, 'FaceColor', vm_color, 'FaceAlpha', .30, 'EdgeColor', 'none');
    ylabel('Velocity (deg/s)');

    %% ---------- LOW ----------
    subplot(4,2,(gid-1)*2+2); hold on;
    title([group_name ' — Low reward'])
    xlabel('Time (ms)'); ylabel('Rate (Hz)');

    yyaxis left
    for d = 1:8
        tmp = ss_rate(cell_mask, time_ind, d, low_inds, amp_bin);
        tmp = squeeze(mean(tmp, 4, 'omitnan'));
        sig = mean(tmp, 1, 'omitnan');
        sem = std(tmp, 0, 1, 'omitnan') ./ sqrt(sum(cell_mask));

        [h1,h2] = boundedline(time_axis, sig, sem, 'color', colors_(d,:), 'alpha');
        h1.DisplayName = leg_{d};
        h2.HandleVisibility = 'off';
    end
    xline(0,'--','HandleVisibility','off');

    yyaxis right
    vtmp = vm_trace(cell_mask, time_ind, :, low_inds, amp_bin);
    vtmp = mean(vtmp,[3 4],'omitnan');
    vtmp = mean(vtmp,1,'omitnan');
    sig_v = squeeze(vtmp);
    area(time_axis, sig_v, 'FaceColor', vm_color, 'FaceAlpha', .30, ...
         'EdgeColor', 'none','HandleVisibility','off');
    ylabel('Velocity (deg/s)');
    legend('Location','northeast','Box','off')

%% ---------- NEW SUBPLOTS: CS-on − CS-off grouped by reward ----------

if gid == 2
    % === BURSTERS subplot ===
    subplot(4,2,5); hold on;
    title('Bursters: CS-on − CS-off')
    xlabel('Time (ms)'); ylabel('Δ Rate (Hz)');

    % High reward
    tmpON_b   = squeeze(mean(ss_rate(ind_b, time_ind, 1, high_inds, amp_bin),4,'omitnan'));
    tmpOFF_b  = squeeze(mean(ss_rate(ind_b, time_ind, 5, high_inds, amp_bin),4,'omitnan'));
    diff_b_hi = mean(tmpON_b - tmpOFF_b,1,'omitnan');
    sem_b_hi  = std(tmpON_b - tmpOFF_b,0,1,'omitnan') ./ sqrt(sum(ind_b));

    [h1,h2] = boundedline(time_axis, diff_b_hi, sem_b_hi, 'r','alpha');
    h1.DisplayName = 'High reward'; h2.HandleVisibility = 'off';

    % Low reward
    tmpON_b   = squeeze(mean(ss_rate(ind_b, time_ind, 1, low_inds, amp_bin),4,'omitnan'));
    tmpOFF_b  = squeeze(mean(ss_rate(ind_b, time_ind, 5, low_inds, amp_bin),4,'omitnan'));
    diff_b_lo = mean(tmpON_b - tmpOFF_b,1,'omitnan');
    sem_b_lo  = std(tmpON_b - tmpOFF_b,0,1,'omitnan') ./ sqrt(sum(ind_b));

    [h3,h4] = boundedline(time_axis, diff_b_lo, sem_b_lo, 'b','alpha');
    h3.DisplayName = 'Low reward'; h4.HandleVisibility = 'off';

    xline(0,'--','HandleVisibility','off');
    legend([h1 h3], {'High','Low'}, 'Location','northeast','Box','off');

    % === PAUSERS subplot ===
    subplot(4,2,6); hold on;
    title('Pausers: CS-on − CS-off')
    xlabel('Time (ms)'); ylabel('Δ Rate (Hz)');

    % High reward
    tmpON_p   = squeeze(mean(ss_rate(ind_p, time_ind, 1, high_inds, amp_bin),4,'omitnan'));
    tmpOFF_p  = squeeze(mean(ss_rate(ind_p, time_ind, 5, high_inds, amp_bin),4,'omitnan'));
    diff_p_hi = mean(tmpON_p - tmpOFF_p,1,'omitnan');
    sem_p_hi  = std(tmpON_p - tmpOFF_p,0,1,'omitnan') ./ sqrt(sum(ind_p));

    [h1,h2] = boundedline(time_axis, diff_p_hi, sem_p_hi, 'r','alpha');
    h1.DisplayName = 'High reward'; h2.HandleVisibility = 'off';

    % Low reward
    tmpON_p   = squeeze(mean(ss_rate(ind_p, time_ind, 1, low_inds, amp_bin),4,'omitnan'));
    tmpOFF_p  = squeeze(mean(ss_rate(ind_p, time_ind, 5, low_inds, amp_bin),4,'omitnan'));
    diff_p_lo = mean(tmpON_p - tmpOFF_p,1,'omitnan');
    sem_p_lo  = std(tmpON_p - tmpOFF_p,0,1,'omitnan') ./ sqrt(sum(ind_p));

    [h3,h4] = boundedline(time_axis, diff_p_lo, sem_p_lo, 'b','alpha');
    h3.DisplayName = 'Low reward'; h4.HandleVisibility = 'off';

    xline(0,'--','HandleVisibility','off');
    legend([h1 h3], {'High','Low'}, 'Location','northeast','Box','off');
end

%% ---------- High − Low at CS-on and CS+180 ----------
if gid == 2
    % ---------- BURSTERS ----------
    subplot(4,2,7); hold on;
    title('Bursters: High − Low  |  CS-on vs CS+180')
    xlabel('Time (ms)'); ylabel('Δ Rate (Hz)');

    % CS-on
    tmp_hi_on = squeeze(mean(ss_rate(ind_b, time_ind, 1, high_inds, amp_bin), 4, 'omitnan'));
    tmp_lo_on = squeeze(mean(ss_rate(ind_b, time_ind, 1, low_inds,  amp_bin), 4, 'omitnan'));
    diff_on_b = mean(tmp_hi_on - tmp_lo_on, 1, 'omitnan');
    sem_on_b  = std(tmp_hi_on - tmp_lo_on, 0, 1, 'omitnan') ./ sqrt(sum(ind_b));

    [h1, h2] = boundedline(time_axis, diff_on_b, sem_on_b, 'r', 'alpha');
    h1.DisplayName = 'CS-on'; h2.HandleVisibility = 'off';

    % CS+180
    tmp_hi_180 = squeeze(mean(ss_rate(ind_b, time_ind, 5, high_inds, amp_bin), 4, 'omitnan'));
    tmp_lo_180 = squeeze(mean(ss_rate(ind_b, time_ind, 5, low_inds,  amp_bin), 4, 'omitnan'));
    diff_180_b = mean(tmp_hi_180 - tmp_lo_180, 1, 'omitnan');
    sem_180_b  = std(tmp_hi_180 - tmp_lo_180, 0, 1, 'omitnan') ./ sqrt(sum(ind_b));

    [h3, h4] = boundedline(time_axis, diff_180_b, sem_180_b, 'b', 'alpha');
    h3.DisplayName = 'CS+180'; h4.HandleVisibility = 'off';

    xline(0, '--', 'HandleVisibility', 'off');
    legend([h1 h3], {'CS-on','CS+180'}, 'Location', 'northeast', 'Box', 'off');

    % ---------- PAUSERS ----------
    subplot(4,2,8); hold on;
    title('Pausers: High − Low  |  CS-on vs CS+180')
    xlabel('Time (ms)'); ylabel('Δ Rate (Hz)');

    % CS-on
    tmp_hi_on = squeeze(mean(ss_rate(ind_p, time_ind, 1, high_inds, amp_bin), 4, 'omitnan'));
    tmp_lo_on = squeeze(mean(ss_rate(ind_p, time_ind, 1, low_inds,  amp_bin), 4, 'omitnan'));
    diff_on_p = mean(tmp_hi_on - tmp_lo_on, 1, 'omitnan');
    sem_on_p  = std(tmp_hi_on - tmp_lo_on, 0, 1, 'omitnan') ./ sqrt(sum(ind_p));

    [h1, h2] = boundedline(time_axis, diff_on_p, sem_on_p, 'r', 'alpha');
    h1.DisplayName = 'CS-on'; h2.HandleVisibility = 'off';

    % CS+180
    tmp_hi_180 = squeeze(mean(ss_rate(ind_p, time_ind, 5, high_inds, amp_bin), 4, 'omitnan'));
    tmp_lo_180 = squeeze(mean(ss_rate(ind_p, time_ind, 5, low_inds,  amp_bin), 4, 'omitnan'));
    diff_180_p = mean(tmp_hi_180 - tmp_lo_180, 1, 'omitnan');
    sem_180_p  = std(tmp_hi_180 - tmp_lo_180, 0, 1, 'omitnan') ./ sqrt(sum(ind_p));

    [h3, h4] = boundedline(time_axis, diff_180_p, sem_180_p, 'b', 'alpha');
    h3.DisplayName = 'CS+180'; h4.HandleVisibility = 'off';

    xline(0, '--', 'HandleVisibility', 'off');
    legend([h1 h3], {'CS-on','CS+180'}, 'Location', 'northeast', 'Box', 'off');
end

end
