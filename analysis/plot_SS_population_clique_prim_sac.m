function plot_SS_population_clique_prim_sac(data_path)

 img_save_path = fullfile(data_path,'population_figs\ephys');
 save_path = fullfile(data_path,'population_data');
 
JDM_params_funcs
data_ss  = load(fullfile(save_path,'SS_population_clique_prim_sac.mat'), 'data');
data_ml1 = load(fullfile(save_path,'MLI_population_clique_prim_sac.mat'), 'data');
data_ml2 = load(fullfile(save_path,'MLI2_population_clique_prim_sac.mat'), 'data');
load(fullfile(save_path,'purkinje_cell_ids.mat'),'ind_p','ind_b')
 
%data_cs = load(fullfile(save_path,'CS_population_clique_prim_sac'), 'data');

data_ss  = data_ss.data;
data_ml1 = data_ml1.data;
data_ml2 = data_ml2.data;

sess_SS      = cell2mat(cellfun(@(x) str2double(x(1:6)),...
    data_ss.cell_ids_tot,'UniformOutput',false));
sess_MLI     = cell2mat(cellfun(@(x) str2double(x(1:6)),...
    data_ml1.cell_ids_tot,'UniformOutput',false));
sess_MLI2    = cell2mat(cellfun(@(x) str2double(x(1:6)),...
    data_ml2.cell_ids_tot,'UniformOutput',false));

cliques_SS   = sess_SS.*10+data_ss.cliques;
cliques_MLI  = sess_MLI.*10+data_ml1.cliques;
cliques_MLI2 = sess_MLI2.*10+data_ml2.cliques;

num_SS   = numel(cliques_SS);
num_MLI  = numel(cliques_MLI);
num_MLI2 = numel(cliques_MLI2);
num_b = sum(ind_b);
num_p = sum(ind_p);

SS_cs_on_angle = data_ss.cs_on_ang_tot;
ss_cs_rho      = data_ss.cs_on_rho_tot;
ml1_cs_rho     = data_ml1.cs_on_rho_tot;
ml2_cs_rho     = data_ml2.cs_on_rho_tot;
%% plot SS MLI1 and MLI2 plots
vm_color = [204, 174, 98]/256;
vm_lim   = [0,600];
time_ind = (-50:150)+250;
rate_lim = [10,110];
proj_lim = [-30,60];

cond_labels = {'High', 'Low'};
cs_order_labels = {'cs±180','cs±135','cs±90','cs±45','cs-on'};
colors_ = cool(5); 

% SS rate and vm of high
data_ss_high  = mean(data_ss.rate_tot_sac(:,:,:, [1 2]), 4, 'omitnan');
data_vm_high = mean(data_ss.vm_tot_sac(:,:,:, [1 2]), 4, 'omitnan');
% SS rate and vm of low
data_ss_low =mean(data_ss.rate_tot_sac(:,:,:, [3 4]), 4, 'omitnan');
data_vm_low = mean(data_ss.vm_tot_sac(:,:,:, [3 4]), 4, 'omitnan');
% MLI1 rate and vm
data_ml1_high  = mean(data_ml1.rate_tot_sac(:,:,:, [1 2]), 4, 'omitnan');
data_ml1_low  = mean(data_ml1.rate_tot_sac(:,:,:, [3 4]), 4, 'omitnan');
% MLI2 rate and vm
data_ml2_high  = mean(data_ml2.rate_tot_sac(:,:,:, [1 2]), 4, 'omitnan');
data_ml2_low  = mean(data_ml2.rate_tot_sac(:,:,:, [3 4]), 4, 'omitnan');


% === SS Plot ===
fig = figure;
tiledlayout(4,2,'TileSpacing','compact','Padding','compact')
colororder({'k','k'})
%=== Bursterss Plot ===
for i_group = 1:2
    nexttile
    hold on;
    title([cond_labels{i_group} '-Bursters ' 'N = ' num2str(num_b)])

    for counter_dir = 5:-1:1
        current_dir = counter_dir;
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            current_dir_ = [current_dir, 10-current_dir];
        end

        if i_group == 1
            current_data = data_ss_high;
            current_vm = data_vm_high;
        else
            current_data = data_ss_low;
            current_vm = data_vm_low;
        end

        current_ss_rate_sac_b = squeeze(mean(current_data(ind_b, time_ind, current_dir_), 3, 'omitnan'));
        sig_b = squeeze(mean(current_ss_rate_sac_b,'omitnan'));
        sig_se_b = squeeze(std(current_ss_rate_sac_b,'omitnan'))./sqrt(num_b);

        yyaxis left
        [hl1_b, hl2_b] = boundedline(time_ind-250, sig_b, sig_se_b, 'Color', colors_(6-counter_dir,:), 'alpha');
        hl2_b.HandleVisibility = 'off';
        hl1_b.DisplayName = cs_order_labels{counter_dir};
        ylim([40 90])

        if counter_dir == 1 && i_group == 1
            ylabel('rate(Hz)')
        end

        if counter_dir == 1
            yyaxis right
            sig = mean(mean(current_vm(ind_b, time_ind, :), 3, 'omitnan'));
            h_area = area(time_ind-250,sig,'FaceColor',vm_color,'FaceAlpha',.3,'EdgeColor','none');
            h_area.HandleVisibility = 'off';
            xline(0,'--','HandleVisibility','off');
            ylim(vm_lim);
        end

        if counter_dir == 1 && i_group == 2
            ylabel('Velocity (deg/sec)')
        end
    end
end

%=== Pausers Plot ===
for i_group = 1:2
    nexttile
    title([cond_labels{i_group} '-Pausers ' 'N = ' num2str(num_p)])
    hold on;
    for counter_dir = 5:-1:1
        current_dir = counter_dir;
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            current_dir_ = [current_dir, 10-current_dir];
        end

        if i_group == 1
            current_data = data_ss_high;
            current_vm = data_vm_high;
        else
            current_data = data_ss_low;
            current_vm = data_vm_low;
        end

        current_ss_rate_sac_p = squeeze(mean(current_data(ind_p, time_ind, current_dir_), 3, 'omitnan'));
        sig_b = squeeze(mean(current_ss_rate_sac_p,'omitnan'));
        sig_se_b = squeeze(std(current_ss_rate_sac_p,'omitnan'))./sqrt(num_p);

        yyaxis left
        [hl1_b, hl2_b] = boundedline(time_ind-250, sig_b, sig_se_b, 'Color', colors_(6-counter_dir,:), 'alpha');
        hl2_b.HandleVisibility = 'off';
        hl1_b.DisplayName =  cs_order_labels{counter_dir};
        ylim([35 70])

        if counter_dir == 1 && i_group == 1
            ylabel('rate(Hz)')
        end

        if counter_dir == 1
            yyaxis right
            sig = mean(squeeze(mean(current_vm(ind_p, time_ind, :), 3, 'omitnan')));
            h_area = area(time_ind-250,sig,'FaceColor',vm_color,'FaceAlpha',.3,'EdgeColor','none');
            h_area.HandleVisibility = 'off';
            xline(0,'--','HandleVisibility','off');
            ylim(vm_lim);
        end

        if counter_dir == 1 && i_group == 2
            ylabel('Velocity (deg/sec)')
        end
    end
end

%=== MLI1 Plot ===
colororder({'k','k'})
for i_group = 1:2
    nexttile
    title([cond_labels{i_group} '-MLI1 ' 'N = ' num2str(num_MLI)])
    hold on;
    for counter_dir = 5:-1:1
        current_dir = counter_dir;
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            current_dir_ = [current_dir, 10-current_dir];
        end

        if i_group == 1
            current_data = data_ml1_high;
        else
            current_data = data_ml1_low;
        end

        current_rate = squeeze(mean(current_data(:, time_ind, current_dir_), 3, 'omitnan'));
        sig_b = squeeze(mean(current_rate,'omitnan'));
        sig_se_b = squeeze(std(current_rate,'omitnan'))./sqrt(num_MLI);

        yyaxis left
        [hl1_b, hl2_b] = boundedline(time_ind-250, sig_b, sig_se_b, 'Color', colors_(6-counter_dir,:), 'alpha');
        hl2_b.HandleVisibility = 'off';
        hl1_b.DisplayName =  cs_order_labels{counter_dir};
        ylim([20 110])

        if counter_dir == 1 && i_group == 1
            ylabel('rate(Hz)')
        end

        if counter_dir == 1
            yyaxis right
            sig = mean(mean(current_vm(:, time_ind, :), 3, 'omitnan'));
            h_area = area(time_ind-250,sig,'FaceColor',vm_color,'FaceAlpha',.3,'EdgeColor','none');
            h_area.HandleVisibility = 'off';
            xline(0,'--','HandleVisibility','off');
            ylim(vm_lim);
        end

        if counter_dir == 1 && i_group == 2
            ylabel('Velocity (deg/sec)')
        end
    end
end

%=== MLI2 Plot ===
colororder({'k','k'})
for i_group = 1:2
    nexttile
    title([cond_labels{i_group} '-MLI2 ' 'N = ' num2str(num_MLI2)])
    hold on;
    for counter_dir = 5:-1:1
        current_dir = counter_dir;
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            current_dir_ = [current_dir, 10-current_dir];
        end

        if i_group == 1
            current_data = data_ml2_high;
        else
            current_data = data_ml2_low;
        end

        current_rate = squeeze(mean(current_data(:, time_ind, current_dir_), 3, 'omitnan'));
        sig_b = squeeze(mean(current_rate,'omitnan'));
        sig_se_b = squeeze(std(current_rate,'omitnan'))./sqrt(num_MLI2);

        yyaxis left
        [hl1_b, hl2_b] = boundedline(time_ind-250, sig_b, sig_se_b, 'Color', colors_(6-counter_dir,:), 'alpha');
        hl2_b.HandleVisibility = 'off';
        hl1_b.DisplayName =  cs_order_labels{counter_dir};
        ylim([10 55])
        xlabel('Time (ms)')

        if counter_dir == 1 && i_group == 1
            ylabel('rate(Hz)')
        end

        if counter_dir == 1
            yyaxis right
            sig = mean(mean(current_vm(:, time_ind, :), 3, 'omitnan'));
            h_area = area(time_ind-250,sig,'FaceColor',vm_color,'FaceAlpha',.3,'EdgeColor','none');
            h_area.HandleVisibility = 'off';
            xline(0,'--','HandleVisibility','off');
            ylim(vm_lim);
        end

        if counter_dir == 1 && i_group == 2
            ylabel('Velocity (deg/sec)')
        end
    end
end

lgd = legend(cs_order_labels, 'Position', [0.4 0.01 0.2 0.05], 'Orientation','horizontal');
lgd.Layout.Tile = 'north';
ESN_Beautify_Plot(fig,[10,12])
print(fig, fullfile(img_save_path,'ss-rate-clique-prim-sac.pdf'), '-dpdf', '-bestfit');


% projection to CS-on and CS+90
% === SS Bursters and Pausers Plot ===

 
% %% Axis of symmetry: Null Cancelation
% for i_group = 1:2
%     if i_group == 1
%         current_data = data_ss.rate_tot_sac(:,:,:, [1 2]); % High
%         group_name = 'High';
%     else
%         current_data = data_ss.rate_tot_sac(:,:,:, [3 4]); % Low
%         group_name = 'Low';
%     end
% 
%     current_data = mean(current_data, 4, 'omitnan');  % [neurons x time x 8 directions]
% 
%     fig = figure;
%     n_col = 4;
%     n_row = 2;
% 
%     t_ind = (-100:200)+250;
%     n_t = length(t_ind);
% 
%     % -------- Null Axis: Bursters --------
%     subplot(n_row, n_col, 1)
%     current_rate = current_data(ind_b, t_ind, :);
%     num_b = sum(ind_b);
%     rate_ss = squeeze(mean(current_rate));
%     rate_se = squeeze(std(current_rate)) ./ sqrt(num_b);
%     rate_ss = reshape([rate_ss(:,1:5), rate_ss(:,1), rate_ss(:,8:-1:6), rate_ss(:,5)], 1, []);
%     rate_se = reshape([rate_se(:,1:5), rate_se(:,1), rate_se(:,8:-1:6), rate_se(:,5)], 1, []);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind, rate_ss(ind), rate_se(ind), 'c', 'alpha');
%         hold on;
%         ind2 = ind + n_t*5;
%         boundedline(ind, rate_ss(ind2), rate_se(ind2), 'm', 'alpha');
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     ylim([30,110])
%     ylabel('rate (Hz)')
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     text(10:n_t:5*n_t,45*ones(5,1),{'\theta','\theta+45','\theta+90','\theta+135','\theta+180'},'Color','c');
%     text(10:n_t:5*n_t,35*ones(5,1),{'\theta','\theta-45','\theta-90','\theta-135','\theta-180'},'Color','m');
%     title('PCs burster (null axis)')
% 
%     subplot(n_row, n_col, n_col+1)
%     rate_diff = current_rate(:,:,1:5) - cat(3, current_rate(:,:,1), current_rate(:,:,8:-1:6), current_rate(:,:,5));
%     sig = reshape(mean(rate_diff), 1, []);
%     sig_se = reshape(std(rate_diff), 1, []) ./ sqrt(num_b);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind, sig(ind), sig_se(ind), 'k', 'alpha')
%         hold on;
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     yline(0)
%     ylim([30,110]-70)
%     title('Difference')
% 
%     % -------- Null Axis: Pausers --------
%     subplot(n_row, n_col, 2)
%     current_rate = current_data(ind_p, t_ind, :);
%     num_p = sum(ind_p);
%     rate_ss = squeeze(mean(current_rate));
%     rate_se = squeeze(std(current_rate)) ./ sqrt(num_p);
%     rate_ss = reshape([rate_ss(:,1:5), rate_ss(:,1), rate_ss(:,8:-1:6), rate_ss(:,5)], 1, []);
%     rate_se = reshape([rate_se(:,1:5), rate_se(:,1), rate_se(:,8:-1:6), rate_se(:,5)], 1, []);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind, rate_ss(ind), rate_se(ind), 'c', 'alpha');
%         hold on;
%         ind2 = ind + n_t*5;
%         boundedline(ind, rate_ss(ind2), rate_se(ind2), 'm', 'alpha');
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     ylim([30,110])
%     ylabel('rate (Hz)')
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     text(10:n_t:5*n_t,105*ones(5,1),{'\theta','\theta+45','\theta+90','\theta+135','\theta+180'},'Color','c');
%     text(10:n_t:5*n_t,95*ones(5,1),{'\theta','\theta-45','\theta-90','\theta-135','\theta-180'},'Color','m');
%     title('PCs pauser (null axis)')
% 
%     subplot(n_row, n_col, n_col+2)
%     rate_diff = current_rate(:,:,1:5) - cat(3, current_rate(:,:,1), current_rate(:,:,8:-1:6), current_rate(:,:,5));
%     sig = reshape(mean(rate_diff), 1, []);
%     sig_se = reshape(std(rate_diff), 1, []) ./ sqrt(num_p);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind, sig(ind), sig_se(ind), 'k', 'alpha')
%         hold on;
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     yline(0)
%     ylim([30,110]-70)
%     title('Difference')
% 
%     % -------- Bursters Potent Axis --------
%     subplot(n_row, n_col, 3)
%     current_rate = current_data(ind_b, t_ind, :);
%     rate_ss = squeeze(mean(current_rate));
%     rate_se = squeeze(std(current_rate)) ./ sqrt(num_b);
%     rate_ss = reshape([rate_ss(:,1:5), rate_ss(:,5:8), rate_ss(:,1)], 1, []);
%     rate_se = reshape([rate_se(:,1:5), rate_se(:,5:8), rate_se(:,1)], 1, []);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind, rate_ss(ind), rate_se(ind), 'c', 'alpha');
%         hold on;
%         ind2 = ind + n_t*5;
%         boundedline(ind, rate_ss(ind2), rate_se(ind2), 'm', 'alpha');
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     ylim([30,110])
%     ylabel('rate (Hz)')
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     text(10:n_t:5*n_t,45*ones(5,1),{'\theta','\theta+45','\theta+90','\theta+135','\theta+180'},'Color','c');
%     text(10:n_t:5*n_t,35*ones(5,1),{'\theta+180','\theta+225','\theta+270','\theta+315','\theta'},'Color','m');
%     title('PCs burster (potent axis)')
%     %==  Burster Potent Axis Difference==
%     subplot(n_row, n_col, n_col+3)
%     rate_diff = current_rate(:,:,1:5) - cat(3, current_rate(:,:,5:8), current_rate(:,:,1));
%     sig = reshape(mean(rate_diff), 1, []);
%     sig_se = reshape(std(rate_diff), 1, []) ./ sqrt(num_b);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind, sig(ind), sig_se(ind), 'k', 'alpha')
%         hold on;
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     yline(0)
%     ylim([30,110]-70)
%     title('Difference')
%     
%     % Pausers — Potent Axis
%     subplot(n_row, n_col, 4)
%     current_rate = current_data(ind_p, t_ind, :);
%     rate_ss = squeeze(mean(current_rate));
%     rate_se = squeeze(std(current_rate)) ./ sqrt(num_p);
%     rate_ss = reshape([rate_ss(:,1:5), rate_ss(:,5:8), rate_ss(:,1)], 1, []);
%     rate_se = reshape([rate_se(:,1:5), rate_se(:,5:8), rate_se(:,1)], 1, []);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind, rate_ss(ind), rate_se(ind), 'c', 'alpha');
%         hold on;
%         ind2 = ind + n_t*5;
%         boundedline(ind, rate_ss(ind2), rate_se(ind2), 'm', 'alpha');
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     ylim([30,110])
%     ylabel('rate (Hz)')
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     text(10:n_t:5*n_t,105*ones(5,1),{'\theta','\theta+45','\theta+90','\theta+135','\theta+180'},'Color','c');
%     text(10:n_t:5*n_t,95*ones(5,1),{'\theta+180','\theta+225','\theta+270','\theta+315','\theta'},'Color','m');
%     title('PCs pauser (potent axis)')
%     %==  Pausers Potent Axis Difference==
%     subplot(n_row, n_col, n_col+4)
%     rate_diff = current_rate(:,:,1:5) - cat(3, current_rate(:,:,5:8), current_rate(:,:,1));
%     sig = reshape(mean(rate_diff), 1, []);
%     sig_se = reshape(std(rate_diff), 1, []) ./ sqrt(num_p);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind, sig(ind), sig_se(ind), 'k', 'alpha')
%         hold on;
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     yline(0)
%     ylim([30,110]-70)
%     title('Difference')
% 
%     sgtitle(['Symmetry Analysis — ' group_name])
%     ESN_Beautify_Plot(fig, [16, 5])
%     
% end
% 
% %% Axis of symmetry: Null Cancelation
% 
% for i_group = 1:2
%      if i_group == 1
%         ss_data = data_ss_high;
%         ml1_data = data_ml1_high;
%         ml2_data = data_ml2_high;
%     else
%         ss_data = data_ss_low;
%         ml1_data = data_ml1_low;
%         ml2_data = data_ml2_low;
%      end
% 
%     fig = figure;
%     n_col = 6;
%     n_row = 2;
%     
%     % Axis of symmetry: Null Cancelation
%     % SS
%     subplot(n_row,n_col,1)
%     t_ind = (-100:200)+250;
%     n_t = length(t_ind);
%     current_rate = ss_data(:,t_ind,:);
%     rate_ss = squeeze(mean(current_rate));
%     rate_se = squeeze(std(current_rate))./sqrt(size(current_rate,1));
%     rate_ss = reshape([rate_ss(:,1:5),rate_ss(:,1),rate_ss(:,8:-1:6),rate_ss(:,5)],1,[]);
%     rate_se = reshape([rate_se(:,1:5),rate_se(:,1),rate_se(:,8:-1:6),rate_se(:,5)],1,[]);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,rate_ss(ind),rate_se(ind),'c','alpha')
%         hold on;
%         ind2 = (1:n_t)+(counter-1)*n_t+n_t*5;
%         boundedline(ind,rate_ss(ind2),rate_se(ind2),'m','alpha')
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     ylim([0,110])
%     ylabel('rate (Hz)')
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     text(10:n_t:5*n_t,45*ones(5,1),{'\theta','\theta+45','\theta+90','\theta+135','\theta+180'},'Color','c');
%     text(10:n_t:5*n_t,35*ones(5,1),{'\theta','\theta-45','\theta-90','\theta-135','\theta-180'},'Color','m');
%     title('Pcells')
%     
%     subplot(n_row,n_col,n_col+1)
%     rate_diff = current_rate(:,:,1:5)-cat(3,current_rate(:,:,1),current_rate(:,:,8:-1:6),current_rate(:,:,5));
%     sig = reshape(mean(rate_diff),1,[]);
%     sig_se = reshape(std(rate_diff),1,[])./sqrt(size(current_rate,1));
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,sig(ind),sig_se(ind),'k','alpha')
%         hold on;
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     yline(0)
%     ylim([0,110]-55)
%     
%     % MLI1
%     subplot(n_row,n_col,2)
%     t_ind = (-100:200)+250;
%     n_t = length(t_ind);
%     current_rate = ml1_data(:,t_ind,:);
%     rate_ss = squeeze(mean(current_rate));
%     rate_se = squeeze(std(current_rate))./sqrt(size(current_rate,1));
%     rate_ss = reshape([rate_ss(:,1:5),rate_ss(:,1),rate_ss(:,8:-1:6),rate_ss(:,5)],1,[]);
%     rate_se = reshape([rate_se(:,1:5),rate_se(:,1),rate_se(:,8:-1:6),rate_se(:,5)],1,[]);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,rate_ss(ind),rate_se(ind),'c','alpha')
%         hold on;
%         ind2 = (1:n_t)+(counter-1)*n_t+n_t*5;
%         boundedline(ind,rate_ss(ind2),rate_se(ind2),'m','alpha')
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     ylim([0,110])
%     ylabel('rate (Hz)')
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     text(10:n_t:5*n_t,105*ones(5,1),{'\theta','\theta+45','\theta+90','\theta+135','\theta+180'},'Color','c');
%     text(10:n_t:5*n_t,95*ones(5,1),{'\theta','\theta-45','\theta-90','\theta-135','\theta-180'},'Color','m');
%     
%     title('MLI1')
%     
%     subplot(n_row,n_col,n_col+2)
%     rate_diff = current_rate(:,:,1:5)-cat(3,current_rate(:,:,1),current_rate(:,:,8:-1:6),current_rate(:,:,5));
%     sig = reshape(mean(rate_diff),1,[]);
%     sig_se = reshape(std(rate_diff),1,[])./sqrt(size(current_rate,1));
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,sig(ind),sig_se(ind),'k','alpha')
%         hold on;
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     yline(0)
%     ylim([0,110]-55)
%     
%     % MLI2
%     subplot(n_row,n_col,3)
%     t_ind = (-100:200)+250;
%     n_t = length(t_ind);
%     current_rate = ml2_data(:,t_ind,:);
%     rate_ss = squeeze(mean(current_rate));
%     rate_se = squeeze(std(current_rate))./sqrt(size(current_rate,1));
%     rate_ss = reshape([rate_ss(:,1:5),rate_ss(:,1),rate_ss(:,8:-1:6),rate_ss(:,5)],1,[]);
%     rate_se = reshape([rate_se(:,1:5),rate_se(:,1),rate_se(:,8:-1:6),rate_se(:,5)],1,[]);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,rate_ss(ind),rate_se(ind),'c','alpha')
%         hold on;
%         ind2 = (1:n_t)+(counter-1)*n_t+n_t*5;
%         boundedline(ind,rate_ss(ind2),rate_se(ind2),'m','alpha')
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     ylim([0,110])
%     ylabel('rate (Hz)')
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     text(10:n_t:5*n_t,105*ones(5,1),{'\theta','\theta+45','\theta+90','\theta+135','\theta+180'},'Color','c');
%     text(10:n_t:5*n_t,95*ones(5,1),{'\theta','\theta-45','\theta-90','\theta-135','\theta-180'},'Color','m');
%     
%     title('MLI2')
%     
%     subplot(n_row,n_col,n_col+3)
%     rate_diff = current_rate(:,:,1:5)-cat(3,current_rate(:,:,1),current_rate(:,:,8:-1:6),current_rate(:,:,5));
%     sig = reshape(mean(rate_diff),1,[]);
%     sig_se = reshape(std(rate_diff),1,[])./sqrt(size(current_rate,1));
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,sig(ind),sig_se(ind),'k','alpha')
%         hold on;
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     yline(0)
%     ylim([0,110]-55)
%     
%     % Potent axis
%     % SS
%     subplot(n_row,n_col,4)
%     t_ind = (-100:200)+250;
%     n_t = length(t_ind);
%     current_rate = ss_data(:,t_ind,:);
%     rate_ss = squeeze(mean(current_rate));
%     rate_se = squeeze(std(current_rate))./sqrt(size(current_rate,1));
%     rate_ss = reshape([rate_ss(:,1:5),rate_ss(:,5:8),rate_ss(:,1)],1,[]);
%     rate_se = reshape([rate_se(:,1:5),rate_se(:,5:8),rate_se(:,1)],1,[]);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,rate_ss(ind),rate_se(ind),'c','alpha')
%         hold on;
%         ind2 = (1:n_t)+(counter-1)*n_t+n_t*5;
%         boundedline(ind,rate_ss(ind2),rate_se(ind2),'m','alpha')
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     ylim([0,110])
%     ylabel('rate (Hz)')
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     text(10:n_t:5*n_t,45*ones(5,1),{'\theta','\theta+45','\theta+90','\theta+135','\theta+180'},'Color','c');
%     text(10:n_t:5*n_t,35*ones(5,1),{'\theta+180','\theta+225','\theta+270','\theta+315','\theta'},'Color','m');
%     title('Pcells')
%     
%     subplot(n_row,n_col,n_col+4)
%     rate_diff = current_rate(:,:,1:5)-cat(3,current_rate(:,:,5:8),current_rate(:,:,1));
%     sig = reshape(mean(rate_diff),1,[]);
%     sig_se = reshape(std(rate_diff),1,[])./sqrt(num_b);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,sig(ind),sig_se(ind),'k','alpha')
%         hold on;
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     yline(0)
%     ylim([0,110]-55)
%     
%     % MLI1
%     subplot(n_row,n_col,5)
%     t_ind = (-100:200)+250;
%     n_t = length(t_ind);
%     current_rate = ml1_data(:,t_ind,:);
%     rate_ss = squeeze(mean(current_rate));
%     rate_se = squeeze(std(current_rate))./sqrt(size(current_rate,1));
%     rate_ss = reshape([rate_ss(:,1:5),rate_ss(:,5:8),rate_ss(:,1)],1,[]);
%     rate_se = reshape([rate_se(:,1:5),rate_se(:,5:8),rate_se(:,1)],1,[]);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,rate_ss(ind),rate_se(ind),'c','alpha')
%         hold on;
%         ind2 = (1:n_t)+(counter-1)*n_t+n_t*5;
%         boundedline(ind,rate_ss(ind2),rate_se(ind2),'m','alpha')
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     ylim([0,110])
%     ylabel('rate (Hz)')
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     text(10:n_t:5*n_t,105*ones(5,1),{'\theta','\theta+45','\theta+90','\theta+135','\theta+180'},'Color','c');
%     text(10:n_t:5*n_t,95*ones(5,1),{'\theta+180','\theta+225','\theta+270','\theta+315','\theta'},'Color','m');
%     
%     title('MLI1')
%     
%     subplot(n_row,n_col,n_col+5)
%     rate_diff = current_rate(:,:,1:5)-cat(3,current_rate(:,:,5:8),current_rate(:,:,1));
%     sig = reshape(mean(rate_diff),1,[]);
%     sig_se = reshape(std(rate_diff),1,[])./sqrt(num_b);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,sig(ind),sig_se(ind),'k','alpha')
%         hold on;
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     yline(0)
%     ylim([0,110]-55)
%     
%     % MLI2
%     subplot(n_row,n_col,6)
%     t_ind = (-100:200)+250;
%     n_t = length(t_ind);
%     current_rate = ml2_data(:,t_ind,:);
%     rate_ss = squeeze(mean(current_rate));
%     rate_se = squeeze(std(current_rate))./sqrt(size(current_rate,1));
%     rate_ss = reshape([rate_ss(:,1:5),rate_ss(:,5:8),rate_ss(:,1)],1,[]);
%     rate_se = reshape([rate_se(:,1:5),rate_se(:,5:8),rate_se(:,1)],1,[]);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,rate_ss(ind),rate_se(ind),'c','alpha')
%         hold on;
%         ind2 = (1:n_t)+(counter-1)*n_t+n_t*5;
%         boundedline(ind,rate_ss(ind2),rate_se(ind2),'m','alpha')
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     ylim([0,110])
%     ylabel('rate (Hz)')
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     text(10:n_t:5*n_t,105*ones(5,1),{'\theta','\theta+45','\theta+90','\theta+135','\theta+180'},'Color','c');
%     text(10:n_t:5*n_t,95*ones(5,1),{'\theta+180','\theta+225','\theta+270','\theta+315','\theta'},'Color','m');
%     
%     title('MLI2')
%     
%     subplot(n_row,n_col,n_col+6)
%     rate_diff = current_rate(:,:,1:5)-cat(3,current_rate(:,:,5:8),current_rate(:,:,1));
%     sig = reshape(mean(rate_diff),1,[]);
%     sig_se = reshape(std(rate_diff),1,[])./sqrt(num_b);
%     for counter = 1:5
%         ind = (1:n_t)+(counter-1)*n_t;
%         boundedline(ind,sig(ind),sig_se(ind),'k','alpha')
%         hold on;
%     end
%     xticks(1:100:5*n_t)
%     xticklabels(repmat({'','0','100'},1,5));
%     xline(101:n_t-1:5*n_t,'--')
%     xline(0:n_t-1:5*n_t)
%     yline(0)
%     ylim([0,110]-55)
%     sgtitle(['Symmetry Analysis — ' cond_labels{i_group }])
%     ESN_Beautify_Plot(fig,[20,5])
% end

end