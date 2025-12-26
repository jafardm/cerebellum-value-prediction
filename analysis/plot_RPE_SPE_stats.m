function plot_RPE_SPE_stats(data_path)
JDM_params_funcs 

img_save_path = fullfile(data_path,'population_figs\ephys');
load(fullfile(data_path,'population_data','CS_rate_dir_corr_sac_50-150'));

% --- Extract visual and motor tuning ---
cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_data, 'UniformOutput', false));
cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_data, 'UniformOutput', false));

cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, cs_on_data, 'UniformOutput', false));
cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, cs_on_data, 'UniformOutput', false));

num_cs = numel(cs_on_data);

% --- Select preferred angle (visual or saccade) ---
ind_m = cs_on_rho_sac > cs_on_rho_vis;

cs_on_ang = cs_on_ang_vis / 180 * pi;
cs_on_ang(ind_m) = cs_on_ang_sac(ind_m) / 180 * pi;

cs_on_rho = cs_on_rho_vis;
cs_on_rho(ind_m) = cs_on_rho_sac(ind_m);

ang_bins = params.sac.ang_values / 180 * pi;

[~, cs_on_bin] = min(abs(angdiff(repmat(ang_bins, num_cs, 1), cs_on_ang * ones(size(ang_bins)))), [], 2);

order_ang = nan(num_cs, 8);
for counter_cs = 1:num_cs
    order_ang(counter_cs, :) = circshift(1:8, 5 - cs_on_bin(counter_cs));
end

% --- Add wraparound index ---
order_ang = [order_ang, order_ang(:,1)];  % now 9 angles

% --- Rotate visual/saccade firing rates ---
cs_on_rate_vis = cell2mat(cellfun(@(x) x.vis.fr_avg, cs_on_data, 'UniformOutput', false));
cs_on_rate_sac = cell2mat(cellfun(@(x) x.sac.fr_avg, cs_on_data, 'UniformOutput', false));

cs_on_rate_vis_rot = nan(num_cs, 9);
cs_on_rate_sac_rot = nan(num_cs, 9);
cs_rate_rot        = nan(num_cs, 9, 4, 2);  % 9 directions now

for counter_cs = 1:num_cs
    cs_on_rate_vis_rot(counter_cs, :) = cs_on_rate_vis(counter_cs, order_ang(counter_cs, :));
    cs_on_rate_sac_rot(counter_cs, :) = cs_on_rate_sac(counter_cs, order_ang(counter_cs, :));
    cs_rate_rot(counter_cs, :, :, :)  = cs_rate(counter_cs, order_ang(counter_cs, :), :, :);
end

% --------------------
% Plot
% --------------------

cs_rate_HL  = mean(cs_rate_rot(:, :, 2, 1),2);  % reward = high, event = sac
cs_rate_LL  = mean(cs_rate_rot(:, :, 3, 1),2);  % reward = low,  event = sac

cs_rate_HH  = mean(cs_rate_rot(:, :, 1, 1),2);  % reward = high, event = sac
cs_rate_LH  = mean(cs_rate_rot(:, :, 4, 1),2);  % reward = low,  event = sac

LH_HH_diff = cs_rate_LH - cs_rate_HH;
LL_HL_diff = cs_rate_LL - cs_rate_HL;


% ----------------------------
% Violin plot of CS rate differences
% ----------------------------

group_labels = [repmat({'LH - HH'}, numel(LH_HH_diff), 1); ...
                repmat({'LL - HL'}, numel(LL_HL_diff), 1)];
combined_data = [LH_HH_diff(:); LL_HL_diff(:)];

% Create figure
fig=figure;
violinplot(combined_data, group_labels, 'ShowMean', true, 'ViolinColor', [1 0.4 0.4;0.2 0.6 1]);

ylabel('\Delta CS rate (sp/s)');
title('CS rate diff RPE+- on SPE, N = 900');
box off
ESN_Beautify_Plot(fig, [10, 4])
print(fig,fullfile(img_save_path,'cs-rate-diff_stats-RPE-SPE.pdf'),'-dpdf', '-bestfit') 

% ----------------------------
% t-test: 
% ----------------------------
[h, p, ci, stats] = ttest(LL_HL_diff);
fprintf('T-test LL - HL: t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p);

[h, p, ci, stats] = ttest(LH_HH_diff);
fprintf('T-test LH - HH: t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p);
end