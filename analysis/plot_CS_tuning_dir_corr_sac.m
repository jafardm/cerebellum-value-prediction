function plot_CS_tuning_dir_corr_sac(data_path)
% Corrective saccades: CS directional tuning (all cells).
% Produces two PDFs: one for Visual epoch and one for Saccade epoch.
% Left subplot: LH & HH.  Right subplot: LL & HL.
% Each trace is mean ± SEM across cells, aligned to each neuron's CS-on.

JDM_params_funcs;

img_save_path = fullfile(data_path,'population_figs','ephys');
if ~exist(img_save_path,'dir'), mkdir(img_save_path); end

load(fullfile(data_path,'population_data','CS_rate_dir_corr_sac_50-150.mat'));

% ------------------------------
% Prepare directional alignment
% ------------------------------
cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_data, 'UniformOutput', false));
cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_data, 'UniformOutput', false));
cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, cs_on_data, 'UniformOutput', false));
cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, cs_on_data, 'UniformOutput', false));

num_cs = numel(cs_on_data);
use_sac = cs_on_rho_sac > cs_on_rho_vis;

cs_on_ang = cs_on_ang_vis / 180 * pi;
cs_on_ang(use_sac) = cs_on_ang_sac(use_sac) / 180 * pi;

ang_bins = params.sac.ang_values / 180 * pi;  % 8 bins in radians
[~, cs_on_bin] = min(abs(angdiff(repmat(ang_bins, num_cs, 1), cs_on_ang * ones(size(ang_bins)))), [], 2);

order_ang = nan(num_cs, 8);
for i = 1:num_cs
    order_ang(i, :) = circshift(1:8, 5 - cs_on_bin(i)); % CS-on at position 5
end
order_ang = [order_ang, order_ang(:,1)]; % wrap to length 9

% Rotate per-cell curves: cs_rate is [cell x 8dir x 4cond x 2epoch]
cs_rate_rot = nan(num_cs, 9, 4, 2);
for i = 1:num_cs
    cs_rate_rot(i,:,:,:) = cs_rate(i, order_ang(i,:), :, :);
end

% -------- Visual (epoch = 1) --------
make_epoch_figure_CS_corr('Visual','VIS',1, cs_rate_rot, num_cs, img_save_path);

% -------- Saccade (epoch = 2) --------
make_epoch_figure_CS_corr('Saccade','SAC',2, cs_rate_rot, num_cs, img_save_path);

end  % <- end of main function



% ======================================================================
% Subfunction (not nested) to avoid "mismatch with end of functions" error
% ======================================================================
function make_epoch_figure_CS_corr(epoch_name, tag, ep_idx, cs_rate_rot, num_cs, img_save_path)

xaxis   = 1:9;
xlabels = {'cs-180','cs-135','cs-90','cs-45','cs-on','cs+45','cs+90','cs+135','cs+180'};

% Colors
hh_color = [0.8392, 0.1529, 0.1569]; % red
hl_color = [1.0000, 0.4980, 0.0549]; % orange
lh_color = [0.1725, 0.6275, 0.1725]; % green
ll_color = [0.1216, 0.4667, 0.7059]; % blue

% Baseline for this epoch (across dirs & conditions)
tmp = squeeze(cs_rate_rot(:,:,:,ep_idx));  % [cell x dir9 x cond4]
tot_cs_rate   = squeeze(mean(mean(tmp, 3, 'omitnan'), 2, 'omitnan')); % [cell x 1]
mean_baseline = mean(tot_cs_rate, 1, 'omitnan');
sem_baseline  = std(tot_cs_rate, 0, 1, 'omitnan') / sqrt(num_cs);

% Aggregate across ALL cells for this epoch
cs_LL = squeeze(cs_rate_rot(:, :, 3, ep_idx));  % [cell x 9]
cs_LH = squeeze(cs_rate_rot(:, :, 4, ep_idx));
cs_HL = squeeze(cs_rate_rot(:, :, 2, ep_idx));
cs_HH = squeeze(cs_rate_rot(:, :, 1, ep_idx));

n_LL = sum(any(~isnan(cs_LL),2));
n_LH = sum(any(~isnan(cs_LH),2));
n_HL = sum(any(~isnan(cs_HL),2));
n_HH = sum(any(~isnan(cs_HH),2));

mu_LL = mean(cs_LL, 1, 'omitnan'); se_LL = std(cs_LL, 0, 1, 'omitnan') / max(1,sqrt(n_LL));
mu_LH = mean(cs_LH, 1, 'omitnan'); se_LH = std(cs_LH, 0, 1, 'omitnan') / max(1,sqrt(n_LH));
mu_HL = mean(cs_HL, 1, 'omitnan'); se_HL = std(cs_HL, 0, 1, 'omitnan') / max(1,sqrt(n_HL));
mu_HH = mean(cs_HH, 1, 'omitnan'); se_HH = std(cs_HH, 0, 1, 'omitnan') / max(1,sqrt(n_HH));

fig = figure;
set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 10 4], 'PaperSize', [10 4]);

% --- Subplot 1: LH & HH
subplot(1,2,1); hold on; box on
[h1, p1] = boundedline(xaxis, mu_LH, se_LH, 'cmap', lh_color); h1.DisplayName = sprintf('LH (n=%d)', n_LH); p1.HandleVisibility='off';
[h2, p2] = boundedline(xaxis, mu_HH, se_HH, 'cmap', hh_color); h2.DisplayName = sprintf('HH (n=%d)', n_HH); p2.HandleVisibility='off';
[hb1, hb1p] = boundedline(xaxis, repmat(mean_baseline, [1 9]), repmat(sem_baseline, [1 9]), '--', 'cmap', [0 0 0]);
hb1.DisplayName = 'Baseline'; hb1p.HandleVisibility='off';
xticks(xaxis); xticklabels(xlabels);
xlabel('CS-on aligned direction'); ylabel('CS rate (Hz)');
title([epoch_name ' — LH & HH']); legend('Box','off','Location','northwest');

% --- Subplot 2: LL & HL
subplot(1,2,2); hold on; box on
[h3, p3] = boundedline(xaxis, mu_LL, se_LL, 'cmap', ll_color); h3.DisplayName = sprintf('LL (n=%d)', n_LL); p3.HandleVisibility='off';
[h4, p4] = boundedline(xaxis, mu_HL, se_HL, 'cmap', hl_color); h4.DisplayName = sprintf('HL (n=%d)', n_HL); p4.HandleVisibility='off';
[hb2, hb2p] = boundedline(xaxis, repmat(mean_baseline, [1 9]), repmat(sem_baseline, [1 9]), '--', 'cmap', [0 0 0]);
hb2.DisplayName = 'Baseline'; hb2p.HandleVisibility='off';
xticks(xaxis); xticklabels(xlabels);
xlabel('CS-on aligned direction'); ylabel('CS rate (Hz)');
title([epoch_name ' — LL & HL']); legend('Box','off','Location','northwest');

sgtitle(sprintf('%s epoch — CS cells (n = %d)', epoch_name, num_cs), 'FontSize', 10, 'FontWeight', 'bold');

ESN_Beautify_Plot(fig, [10, 4]);
print(fig, fullfile(img_save_path, sprintf('cs-rate-corr-sac_%s_LH-HH_vs_LL-HL.pdf', tag)), '-dpdf', '-bestfit');

end  % <- end of subfunction
