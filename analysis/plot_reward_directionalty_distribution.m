function plot_reward_directionalty_distribution(data_path)
% 2D distribution per cell:
%   x = CS_on - CS+180  (directional contrast, time-windowed, averaged across ALL conditions)
%   y = High - Low at CS_on (reward contrast, time-windowed)
% Two subplots: Visual window and Saccade window.

JDM_params_funcs;

img_save_path = fullfile(data_path,'population_figs','ephys');
if ~exist(img_save_path,'dir'), mkdir(img_save_path); end

load(fullfile(data_path,'population_data','cs_on_rate_prim_sac'));

% ------------------------------
% CS-on alignment
% ------------------------------
cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_data, 'UniformOutput', false));
cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_data, 'UniformOutput', false));
cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, cs_on_data, 'UniformOutput', false));
cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, cs_on_data, 'UniformOutput', false));

num_cs = numel(cs_on_data);
use_sac = cs_on_rho_sac > cs_on_rho_vis;

cs_on_ang = cs_on_ang_vis/180*pi;
cs_on_ang(use_sac) = cs_on_ang_sac(use_sac)/180*pi;

ang_bins = params.sac.ang_values/180*pi;  % 8 bins
[~, cs_on_bin] = min(abs(angdiff(repmat(ang_bins, num_cs, 1), cs_on_ang * ones(size(ang_bins)))), [], 2);

order_ang = nan(num_cs,8);
for c = 1:num_cs
    order_ang(c,:) = circshift(1:8, 5 - cs_on_bin(c));  % CS-on -> index 5
end
order_ang = [order_ang, order_ang(:,1)];                % wrap to 9 (repeat first)

% ------------------------------
% Rotate cs_rate to CS-on frame
% cs_rate dims expected: [cell x time x 8dir x cond x epoch]
% ------------------------------
[~, num_time, ~, num_cond, num_epoch] = size(cs_rate);
num_epoch = min(num_epoch, 2); % keep [visual,saccade] if present

cs_rate_rot = nan(num_cs, num_time, 9, num_cond, num_epoch);
for c = 1:num_cs
    % first 8 directions
    cs_rate_rot(c,:,1:8,:,:) = cs_rate(c, :, order_ang(c,1:8), :, 1:num_epoch);
    % 9th is wrap (+180)
    cs_rate_rot(c,:,9,:,:)   = cs_rate(c, :, order_ang(c,9),   :, 1:num_epoch);
end

% Time base (assumes zero ~ index 250; adjust if needed)
t = (1:num_time) - 250;

% Windows and condition groups
time_window_vis = [50 150];
time_window_sac = [-50 50];
mask_vis = (t >= time_window_vis(1)) & (t <= time_window_vis(2));
mask_sac = (t >= time_window_sac(1)) & (t <= time_window_sac(2));

high_conds = 1;
low_conds  = 2;
high_conds = high_conds(high_conds <= num_cond);
low_conds  = low_conds( low_conds  <= num_cond);

% ------------------------------
% Compute X and Y per cell (VISUAL)
% ------------------------------
vis_allcond = squeeze(mean(cs_rate_rot(:,:,:,:,1), 4, 'omitnan')); % [cell x time x 9]
x_vis = squeeze(mean(vis_allcond(:, mask_vis, 5) - vis_allcond(:, mask_vis, 9), 2, 'omitnan')); % [cell x 1]

vis_high = squeeze(mean(cs_rate_rot(:,:,5, high_conds, 1), 4, 'omitnan')); % [cell x time]
vis_low  = squeeze(mean(cs_rate_rot(:,:,5, low_conds, 1),  4, 'omitnan')); % [cell x time]
y_vis = squeeze(mean(vis_high(:, mask_vis) - vis_low(:, mask_vis), 2, 'omitnan'));

% ------------------------------
% Compute X and Y per cell (SACCADE)
% ------------------------------
sac_allcond = squeeze(mean(cs_rate_rot(:,:,:,:,2), 4, 'omitnan')); % [cell x time x 9]
x_sac = squeeze(mean(sac_allcond(:, mask_sac, 5) - sac_allcond(:, mask_sac, 9), 2, 'omitnan'));

sac_high = squeeze(mean(cs_rate_rot(:,:,5, high_conds, 2), 4, 'omitnan'));
sac_low  = squeeze(mean(cs_rate_rot(:,:,5, low_conds, 2),  4, 'omitnan'));
y_sac = squeeze(mean(sac_high(:, mask_sac) - sac_low(:, mask_sac), 2, 'omitnan'));

% ------------------------------
% Plot
% ------------------------------
fig = figure;
set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 10 4], 'PaperSize', [10 4]);

% Visual
subplot(1,2,1); hold on; box on
scatter(x_vis, y_vis, 18, 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0.2);
xline(0,'--k'); yline(0,'--k');
xlabel(sprintf('CS_{on} - CS_{+180}  (Hz)\n[%d..%d ms]', time_window_vis(1), time_window_vis(2)));
ylabel('High - Low at CS_{on}  (Hz)');
title(sprintf('Visual (n = %d)', num_cs));
axis tight; grid on

% Saccade
subplot(1,2,2); hold on; box on
scatter(x_sac, y_sac, 18, 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0.2);
xline(0,'--k'); yline(0,'--k');
xlabel(sprintf('CS_{on} - CS_{+180}  (Hz)\n[%d..%d ms]', time_window_sac(1), time_window_sac(2)));
ylabel('High - Low at CS_{on}  (Hz)');
title(sprintf('Saccade (n = %d)', num_cs));
axis tight; grid on
sgtitle('2D distribution per cell: Directionality vs Reward contrast');

ESN_Beautify_Plot(fig, [10, 4])
print(fig, fullfile(img_save_path, 'reward_directionality_2D_vis_sac.pdf'), '-dpdf', '-bestfit')

end
