function plot_CS_tuning_dir_prim_sac(data_path)
% Plots CS directional tuning (primary saccade epoch) for the entire cell population.
% Shows Low vs High reward traces aligned to each neuron's CS-on and wrapped
% from cs-180 ... cs+180, plus a population baseline band.

JDM_params_funcs;

img_save_path = fullfile(data_path,'population_figs','ephys');
if ~exist(img_save_path,'dir'), mkdir(img_save_path); end

load(fullfile(data_path,'population_data','CS_rate_dir_prim_sac.mat'));

% --- Choose each neuron's dominant tuning (visual vs saccade)
cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_data, 'UniformOutput', false));
cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_data, 'UniformOutput', false));

cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, cs_on_data, 'UniformOutput', false));
cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, cs_on_data, 'UniformOutput', false));

num_cs = numel(cs_on_data);
use_sac = cs_on_rho_sac > cs_on_rho_vis;

cs_on_ang = cs_on_ang_vis / 180 * pi;
cs_on_ang(use_sac) = cs_on_ang_sac(use_sac) / 180 * pi;

% --- Build rotation to put CS-on at the center
ang_bins = params.sac.ang_values / 180 * pi;           % 8 bins in radians
[~, cs_on_bin] = min(abs(angdiff(repmat(ang_bins, num_cs, 1), cs_on_ang * ones(size(ang_bins)))), [], 2);

order_ang = nan(num_cs, 8);
for c = 1:num_cs
    order_ang(c, :) = circshift(1:8, 5 - cs_on_bin(c));  % put CS-on at position 5
end
order_ang = [order_ang, order_ang(:,1)];                 % wrap to length 9

% --- Rotate per-cell tuning curves
cs_on_rate_sac = cell2mat(cellfun(@(x) x.sac.fr_avg, cs_on_data, 'UniformOutput', false)); % [cell x 8]
cs_on_rate_sac_rot = nan(num_cs, 9);
cs_rate_rot = nan(num_cs, 9, 2, 2);  % [cell x dir9 x reward2 x epoch2]

for c = 1:num_cs
    cs_on_rate_sac_rot(c, :) = cs_on_rate_sac(c, order_ang(c, :));
    cs_rate_rot(c, :, :, :)  = cs_rate(c, order_ang(c, :), :, :);
end

% --- Baseline across directions (per cell), then pop mean±SEM
tot_cs_rate       = mean(cs_on_rate_sac_rot, 2);                 % [cell x 1]
mean_sac_baseline = mean(tot_cs_rate, 1, 'omitnan');
std_sac_baseline  = std(tot_cs_rate, 0, 1, 'omitnan') / sqrt(num_cs);

% --- X axis / labels
xaxis   = 1:9;
xlabels = {'\theta-180','\theta-135','\theta-90','\theta-45',...
    '\theta','\theta+45','\theta+90','\theta+135','\theta+180'};

% --- Colors (Low, High)
colors_ = cool(2);  % colors_(1,:) = Low, colors_(2,:) = High

% --- Aggregate across ALL cells (epoch 1 = visual)
cs_rate_high_vis = squeeze(cs_rate_rot(:, :, 1, 1));   % [cell x 9]
cs_rate_low_vis  = squeeze(cs_rate_rot(:, :, 2, 1));  

mean_high_vis = mean(cs_rate_high_vis, 1, 'omitnan');     
sem_high_vis  = std(cs_rate_high_vis, 0, 1, 'omitnan') / sqrt(num_cs);

mean_low_vis  = mean(cs_rate_low_vis, 1, 'omitnan');
sem_low_vis  = std(cs_rate_low_vis, 0, 1, 'omitnan') / sqrt(num_cs);

% ---  (epoch 2 = saccade panel)
cs_rate_high_sac = squeeze(cs_rate_rot(:, :, 1, 2));   % [cell x 9]
cs_rate_low_sac  = squeeze(cs_rate_rot(:, :, 2, 2));   

mean_high_sac = mean(cs_rate_high_sac, 1, 'omitnan');     
sem_high_sac  = std(cs_rate_high_sac, 0, 1, 'omitnan') / sqrt(num_cs);

mean_low_sac  = mean(cs_rate_low_sac, 1, 'omitnan');
sem_low_sac  = std(cs_rate_low_sac, 0, 1, 'omitnan') / sqrt(num_cs);

% --- Plot
fig = figure;
set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6.5 4], 'PaperSize', [6.5 4]);
subplot(1,2,1)
hold on; box on
[hl1, hp1] = boundedline(xaxis, mean_low_vis,  sem_low_vis,  'cmap', colors_(1,:)); hl1.DisplayName = 'Low';  hp1.HandleVisibility = 'off';
[hl2, hp2] = boundedline(xaxis, mean_high_vis, sem_high_vis, 'cmap', colors_(2,:)); hl2.DisplayName = 'High'; hp2.HandleVisibility = 'off';

% Baseline band
[hb, hbpatch] = boundedline(xaxis, repmat(mean_sac_baseline, [1 9]), ...
                            repmat(std_sac_baseline, [1 9]), '--', 'cmap', [0 0 0]);
hb.DisplayName = 'Baseline'; hbpatch.HandleVisibility = 'off';

xlabel('CS angle (deg)')
ylabel('CS rate (Hz)')
xticks(xaxis); xticklabels(xlabels)
legend('Box','off','Location','northwest')
title(sprintf('primary sac(visual) (50 to 150 ms), (n = %d)', num_cs))


subplot(1,2,2)
hold on; box on
[hl1, hp1] = boundedline(xaxis, mean_low_sac,  sem_low_sac,  'cmap', colors_(1,:)); hl1.DisplayName = 'Low';  hp1.HandleVisibility = 'off';
[hl2, hp2] = boundedline(xaxis, mean_high_sac, sem_high_sac, 'cmap', colors_(2,:)); hl2.DisplayName = 'High'; hp2.HandleVisibility = 'off';

% Baseline band
[hb, hbpatch] = boundedline(xaxis, repmat(mean_sac_baseline, [1 9]), ...
                            repmat(std_sac_baseline, [1 9]), '--', 'cmap', [0 0 0]);
hb.DisplayName = 'Baseline'; hbpatch.HandleVisibility = 'off';

xlabel('CS angle (deg)')
ylabel('CS rate (Hz)')
xticks(xaxis); xticklabels(xlabels)
legend('Box','off','Location','northwest')
title(sprintf('primary sac(saccade)(-50 to 50 ms), (n = %d)', num_cs))

ESN_Beautify_Plot(fig, [8, 4])
print(fig, fullfile(img_save_path, 'Fig1F_cs-dir-tuning.pdf'), '-dpdf', '-bestfit')

end
