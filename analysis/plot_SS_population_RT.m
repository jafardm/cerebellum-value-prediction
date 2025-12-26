function plot_SS_population_RT(data_path)
% Plot SS population rates binned by reaction time (RT).
% For each reward (High/Low), make ONE big figure:
%   Rows: 1 = bursters, 2 = pausers
%   Columns: RT bins
% Each subplot: 5 direction groups, plus eye velocity trace (vm_tot_sac_RT).

    img_save_path = fullfile(data_path,'population_figs','ephys');
    save_path     = fullfile(data_path,'population_data');
     
    if ~exist(img_save_path,'dir')
        mkdir(img_save_path);
    end

    JDM_params_funcs

    % ----- Load SS data (RT-binned) -----
    tmp     = load(fullfile(save_path,'SS_population_clique_prim_sac_RT.mat'), ...
                   'data','meta');
    data_ss = tmp.data;
    meta    = tmp.meta;

    % RT edges (for labeling)
    if isfield(meta,'RT_edges')
        RT_edges  = meta.RT_edges(:)';      % row vector
        n_RT_bins = numel(RT_edges)-1;
    else
        % fallback from data dimension
        [~,~,~,~,n_RT_bins] = size(data_ss.rate_tot_sac_RT);
        RT_edges            = 1:n_RT_bins+1;
    end

    % burster / pauser indices
    load(fullfile(save_path,'purkinje_cell_ids.mat'),'ind_p','ind_b')
    num_b = sum(ind_b);
    num_p = sum(ind_p);

    fprintf('\n================  number of bursters %d \n',num_b)
    fprintf('\n================  number of pausers %d \n',num_p)

    % ----- Plotting params -----
    time_ind = (-50:50) + 250;   % indices around vmax
    t_ms     = time_ind - 250;     % x-axis in ms

    rate_lim = [30, 95];

    vm_color = [204, 174, 98]/256;
    vm_lim   = [0, 600];

    leg_ = {'cs±180','cs±135','cs±90','cs±45','cs-on'};
    colors_ = cool(5);
    colors_ = colors_(end:-1:1,:);   % flip so cs-on is last / warm

    % Sanity: reward / RT dims
    [~, ~, ~, n_rew, n_rt_bins_data] = size(data_ss.rate_tot_sac_RT);
    if n_rt_bins_data ~= n_RT_bins
        warning('RT bin count mismatch between meta and data. Using data dimension.');
        n_RT_bins = n_rt_bins_data;
    end
    if n_rew ~= 2
        warning('Expected 2 reward conditions, found %d.', n_rew);
    end

    % =====================================================================
    % Loop over reward: ONE figure for High, ONE for Low
    % =====================================================================
    for ii = 1:n_rew   % 1=High, 2=Low

        if ii == 1
            rew_label = 'High reward';
        else
            rew_label = 'Low reward';
        end

        fig = figure('Name',sprintf('SS RT %s', rew_label), 'Color','w');
        sgtitle(sprintf('SS population by RT bin — %s', rew_label), 'FontSize', 12)

        % Columns are RT bins
        for iRT = 1:n_RT_bins
            rt_str = sprintf('RT %.0f–%.0f ms', RT_edges(iRT), RT_edges(iRT+1));

            % ---------- Row 1: BURSTERS ----------
            subplot(6, 5, iRT);   % row 1, column iRT
            hold on;

            for counter_dir = 5:-1:1
    
                current_dir = counter_dir;
                if current_dir == 1 || current_dir == 5
                    current_dir_ = current_dir;
                else
                    current_dir_ = [current_dir, 10-current_dir]; % symmetric dirs
                end

                % Extract: [cells x time] for bursters, this dir group, reward ii, RT bin iRT
                current_ss_rate_sac = squeeze( ...
                    mean(data_ss.rate_tot_sac_RT(ind_b, time_ind, current_dir_, ii, iRT), 3)); 
                % current_ss_rate_sac: [num_b x nTime]

                if isempty(current_ss_rate_sac)
                    continue;
                end

                % nan-robust mean and SEM across cells
                sig     = squeeze(mean(current_ss_rate_sac, 1, 'omitnan'));              % [1 x nTime]
                sig_std = squeeze(std(current_ss_rate_sac, 0, 1, 'omitnan'));            % [1 x nTime]
                n_eff   = sum(~isnan(current_ss_rate_sac), 1);                            % [1 x nTime]
                n_eff(n_eff==0) = 1;
                sig_se  = sig_std ./ sqrt(n_eff);

                [hl1, hl2] = boundedline(t_ms, sig, sig_se, 'alpha');
                set(hl1, 'Color', colors_(counter_dir,:));
                set(hl2, 'FaceColor', colors_(counter_dir,:), 'HandleVisibility','off');
                % bursters: no legend labels to avoid clutter
            end

            % Eye velocity (all cells, all directions) for this reward & RT bin
            yyaxis right
            vm_all = squeeze(mean( ...
                        mean(data_ss.vm_tot_sac_RT(:, time_ind, :, ii, iRT), 1, 'omitnan'), ...
                        3, 'omitnan'));   % [1 x nTime]
            if ~all(isnan(vm_all))
                area(t_ms, vm_all, 'FaceColor', vm_color, ...
                     'FaceAlpha', .3, 'EdgeColor','none');
                ylim(vm_lim);
            end
            yyaxis left
            xline(0,'--','HandleVisibility','off');
            ylim(rate_lim);
            if iRT == 1
                ylabel('Rate (Hz) — bursters')
            end
            title(rt_str, 'FontSize', 8);
            if ii == n_rew
                xlabel('Time from max vel. (ms)')
            end

            % ---------- Row 2: PAUSERS ----------
            subplot(6, 5, n_RT_bins + iRT);  % row 2, column iRT
            hold on;

            for counter_dir = 5:-1:1
    
                current_dir = counter_dir;
                if current_dir == 1 || current_dir == 5
                    current_dir_ = current_dir;
                else
                    current_dir_ = [current_dir, 10-current_dir];
                end

                current_ss_rate_sac = squeeze( ...
                    mean(data_ss.rate_tot_sac_RT(ind_p, time_ind, current_dir_, ii, iRT), 3));
                if isempty(current_ss_rate_sac)
                    continue;
                end

                sig     = squeeze(mean(current_ss_rate_sac, 1, 'omitnan'));
                sig_std = squeeze(std(current_ss_rate_sac, 0, 1, 'omitnan'));
                n_eff   = sum(~isnan(current_ss_rate_sac), 1);
                n_eff(n_eff==0) = 1;
                sig_se  = sig_std ./ sqrt(n_eff);
        
                [hl1, hl2] = boundedline(t_ms, sig, sig_se, 'alpha');
                set(hl1, 'Color', colors_(counter_dir,:));
                set(hl2, 'FaceColor', colors_(counter_dir,:), 'HandleVisibility','off');
                hl1.DisplayName = leg_{counter_dir};
            end

            % Eye velocity (same as above, for consistency)
            yyaxis right
            vm_all = squeeze(mean( ...
                        mean(data_ss.vm_tot_sac_RT(:, time_ind, :, ii, iRT), 1, 'omitnan'), ...
                        3, 'omitnan'));   % [1 x nTime]
            if ~all(isnan(vm_all))
                area(t_ms, vm_all, 'FaceColor', vm_color, ...
                     'FaceAlpha', .3, 'EdgeColor','none');
                ylim(vm_lim);
            end
            yyaxis left
            xline(0,'--','HandleVisibility','off');
            ylim(rate_lim);
            if iRT == 1
                ylabel('Rate (Hz) — pausers')
            end
            if ii == n_rew
                xlabel('Time from max vel. (ms)')
            end

            if iRT == n_RT_bins
                legend('Location','bestoutside');
            end
        end

        % save each reward-specific figure
%         fname = sprintf('SS_pop_RT_allbins_%s.png', strrep(rew_label,' ',''));
%         saveas(fig, fullfile(img_save_path, fname));
    end
end
