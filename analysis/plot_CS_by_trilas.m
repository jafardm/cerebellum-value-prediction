function plot_CS_by_trilas(data_path)

MAF_params_funcs 
img_save_path = fullfile(data_path,'population_figs\ephys');
load(fullfile(data_path,'population_data','CS_rate_by_trial'));

colors_1 = cool(4);
t_ind_sac = (-100:250) + 250;

% conditions orders
% HH→HX, HL→HX, LL→LX, LH→LX
num_cs = size(cs_cell_ids,1);
% visual
HH_HX = cs_rate(:,t_ind_sac,1,1);
HL_HX = cs_rate(:,t_ind_sac,2,1);
LL_LX = cs_rate(:,t_ind_sac,3,1);
LH_LX = cs_rate(:,t_ind_sac,4,1);

avg_HH_HX = mean(HH_HX,1,'omitnan');
std_HH_HX = std(HH_HX,0,1,'omitnan')./sqrt(num_cs);

avg_HL_HX = mean(HL_HX,1,'omitnan');
std_HL_HX = std(HL_HX,0,1,'omitnan')./sqrt(num_cs);

avg_LL_LX = mean(LL_LX,1,'omitnan');
std_LL_LX = std(LL_LX,0,1,'omitnan')./sqrt(num_cs);

avg_LH_LX = mean(LH_LX,1,'omitnan');
std_LH_LX = std(LH_LX,0,1,'omitnan')./sqrt(num_cs);

fig = figure;
subplot(1,2,1);
hold on
colororder({'k','k'})
% yyaxis left
[hl1,hl2] = boundedline(t_ind_sac - 250, avg_HH_HX, std_HH_HX, 'Color', colors_1(4,:));
hl1.DisplayName = 'HH-HX';
hl2.HandleVisibility = 'off';

[hl1,hl2] = boundedline(t_ind_sac - 250, avg_HL_HX, std_HL_HX, 'Color', colors_1(3,:));
hl1.DisplayName = 'HL-HX';
hl2.HandleVisibility = 'off';

diff_HH_HL = avg_HH_HX - avg_HL_HX;
h = plot(t_ind_sac- 250,diff_HH_HL,'k');
h.DisplayName = 'HH-HL';
h.HandleVisibility = 'off';

hold off;
xlabel('time from cue onset (ms)');
ylabel(' CS rate (Hz)');
legend('Box','off');
xline(0, '--', 'HandleVisibility', 'off');

% yyaxis right
% area(t_ind_sac - 250, current_cs_vm, ...
%     'FaceColor', vm_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
% ylim(vm_lim);

subplot(1,2,2);
hold on
colororder({'k','k'})
% yyaxis left
[hl1,hl2] = boundedline(t_ind_sac - 250, avg_LL_LX, std_LL_LX, 'Color', colors_1(2,:));
hl1.DisplayName = 'LL-LX';
hl2.HandleVisibility = 'off';

[hl1,hl2] = boundedline(t_ind_sac - 250, avg_LH_LX, std_LH_LX, 'Color', colors_1(1,:));
hl1.DisplayName = 'LH-LX';
hl2.HandleVisibility = 'off';

diff_HH_HL = avg_LL_LX - avg_LH_LX;
h = plot(t_ind_sac- 250,diff_HH_HL,'k');
h.DisplayName = 'LL-LH';
h.HandleVisibility = 'off';

xlabel('time from cue onset (ms)');
legend('Box','off');
xline(0, '--', 'HandleVisibility', 'off');

legend('Box','off','Location','northeast')
hold off
ESN_Beautify_Plot(fig, [10, 4])
sgtitle('primary saccade', 'FontSize', 12, 'FontWeight', 'bold'); 
print(fig,fullfile(img_save_path,'cs-rate-f-by_trial-vis.pdf'),'-dpdf', '-bestfit') 

% sac
HH_HX = cs_rate(:,t_ind_sac,1,2);
HL_HX = cs_rate(:,t_ind_sac,2,2);
LL_LX = cs_rate(:,t_ind_sac,3,2);
LH_LX = cs_rate(:,t_ind_sac,4,2);

avg_HH_HX = mean(HH_HX,1,'omitnan');
std_HH_HX = std(HH_HX,0,1,'omitnan')./sqrt(num_cs);

avg_HL_HX = mean(HL_HX,1,'omitnan');
std_HL_HX = std(HL_HX,0,1,'omitnan')./sqrt(num_cs);

avg_LL_LX = mean(LL_LX,1,'omitnan');
std_LL_LX = std(LL_LX,0,1,'omitnan')./sqrt(num_cs);

avg_LH_LX = mean(LH_LX,1,'omitnan');
std_LH_LX = std(LH_LX,0,1,'omitnan')./sqrt(num_cs);

fig = figure;
subplot(1,2,1);
hold on
colororder({'k','k'})
% yyaxis left
[hl1,hl2] = boundedline(t_ind_sac - 250, avg_HH_HX, std_HH_HX, 'Color', colors_1(4,:));
hl1.DisplayName = 'HH-HX';
hl2.HandleVisibility = 'off';

[hl1,hl2] = boundedline(t_ind_sac - 250, avg_HL_HX, std_HL_HX, 'Color', colors_1(3,:));
hl1.DisplayName = 'HL-HX';
hl2.HandleVisibility = 'off';

diff_HH_HL = avg_HH_HX - avg_HL_HX;
h = plot(t_ind_sac- 250,diff_HH_HL,'k');
h.DisplayName = 'HH-HL';
h.HandleVisibility = 'off';

xlabel('time from sac onset (ms)');
ylabel(' CS rate (Hz)');
legend('Box','off');
xline(0, '--', 'HandleVisibility', 'off');

% yyaxis right
% area(t_ind_sac - 250, current_cs_vm, ...
%     'FaceColor', vm_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
% ylim(vm_lim);

subplot(1,2,2);
hold on
colororder({'k','k'})
% yyaxis left
[hl1,hl2] = boundedline(t_ind_sac - 250, avg_LL_LX, std_LL_LX, 'Color', colors_1(2,:));
hl1.DisplayName = 'LL-LX';
hl2.HandleVisibility = 'off';

[hl1,hl2] = boundedline(t_ind_sac - 250, avg_LH_LX, std_LH_LX, 'Color', colors_1(1,:));
hl1.DisplayName = 'LH-LX';
hl2.HandleVisibility = 'off';

diff_HH_HL = avg_LL_LX - avg_LH_LX;
h = plot(t_ind_sac- 250,diff_HH_HL,'k');
h.DisplayName = 'LL-LH';
h.HandleVisibility = 'off';

hold on;
xlabel('time from sac onset (ms)');
legend('Box','off');
xline(0, '--', 'HandleVisibility', 'off');

legend('Box','off','Location','northeast')
hold off
ESN_Beautify_Plot(fig, [10, 4])
sgtitle('primary saccade', 'FontSize', 12, 'FontWeight', 'bold'); 
print(fig,fullfile(img_save_path,'cs-rate-f-by_trial-sac.pdf'),'-dpdf', '-bestfit') 

end