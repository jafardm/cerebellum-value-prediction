function plot_projected_space_tag6(data_path)

 img_save_path = fullfile(data_path,'population_figs\ephys');
 save_path     = fullfile(data_path,'population_data');
 
JDM_params_funcs
data_ss  = load(fullfile(save_path,'SS_population_clique_tag6_sac.mat'), 'data');
data_ml1 = load(fullfile(save_path,'MLI_population_clique_tag6_sac.mat'), 'data');
load(fullfile(save_path,'purkinje_cell_ids.mat'),'ind_p','ind_b');

data_ss  = data_ss.data;
data_ml1 = data_ml1.data;

sess_SS      = cell2mat(cellfun(@(x) str2double(x(1:6)),...
    data_ss.cell_ids_tot,'UniformOutput',false));
sess_MLI     = cell2mat(cellfun(@(x) str2double(x(1:6)),...
    data_ml1.cell_ids_tot,'UniformOutput',false));

cliques_SS   = sess_SS.*10+data_ss.cliques;
cliques_MLI  = sess_MLI.*10+data_ml1.cliques;

num_SS   = numel(cliques_SS);
num_MLI  = numel(cliques_MLI);

num_b = sum(ind_b);
num_p = sum(ind_p);

ss_cs_rho = data_ss.cs_on_rho_tot;
ml1_cs_rho = data_ml1.cs_on_rho_tot;
%% plot SS MLI1 plots
vm_color = [204, 174, 98]/256;
vm_lim   = [0,400];
time_ind = (-100:150)+250;
proj_lim = [-30,55];

%
% ========= One figure with 1x3 tiles for SS / MLI1 (and optional MLI2) =========
fig = figure;
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
title(t,'Back to center sac. projection — SS & MLI1 (cs-on+90° and cs-on°)')


    % =========================
    % (1) p_y(t): projection onto sac dir + 90°
    % =========================
    ax = nexttile(1);
    hold(ax,'on')
    colororder({'k','k'})
    proj_w = -cosd((0:45:360-45)+90);

    % --- Bursters ---
    rate_pr_b = zeros(num_b,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_b = rate_pr_b + data_ss.rate_tot_sac(ind_b,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_b = rate_pr_b * num_b .* ss_cs_rho(ind_b) / sum(ss_cs_rho(ind_b));

    yyaxis left
    colororder({'k','k'})
    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_b,'omitnan'), ...
        std(rate_pr_b,'omitnan')/sqrt(num_b), 'r','alpha');
    hl.DisplayName = ['SS burst (n = ' num2str(num_b) ')'];
    hp.HandleVisibility = 'off';
    ylabel('Firing rate (Hz)')

    % --- Pausers ---
    rate_pr_p = zeros(num_p,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_p = rate_pr_p + data_ss.rate_tot_sac(ind_p,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_p = rate_pr_p * num_p .* ss_cs_rho(ind_p) / sum(ss_cs_rho(ind_p));

    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_p,'omitnan'), ...
        std(rate_pr_p,'omitnan')/sqrt(num_p), 'b','alpha');
    hl.DisplayName = ['SS pause (n = ' num2str(num_p) ')'];
    hp.HandleVisibility = 'off';

    % --- All SS ---
    rate_pr_ss = zeros(num_SS,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_ss = rate_pr_ss + data_ss.rate_tot_sac(:,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_ss = rate_pr_ss * num_SS .* ss_cs_rho / sum(ss_cs_rho);

    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_ss,'omitnan'), ...
        std(rate_pr_ss,'omitnan')/sqrt(num_SS), 'k','alpha');
    hl.DisplayName = ['SS all (n = ' num2str(num_SS) ')'];
    hp.HandleVisibility = 'off';

    % --- MLI1 ---
    rate_pr_ml1 = zeros(num_MLI, numel(time_ind));
    for counter_dir = 1:8
        rate_pr_ml1 = rate_pr_ml1 + data_ml1.rate_tot_sac(:, time_ind, counter_dir) .* proj_w(counter_dir);
    end
    rate_pr_ml1 = rate_pr_ml1 .* ml1_cs_rho / sum(ml1_cs_rho) * num_MLI;
    [hl, hp] = boundedline(time_ind - 250, ...
        mean(rate_pr_ml1,'omitnan'), ...
        std(rate_pr_ml1,'omitnan') / sqrt(num_MLI), 'g', 'alpha');
    hl.DisplayName = ['MLI1 (n = ' num2str(num_MLI) ')'];
    hp.HandleVisibility = 'off';

    xline(0,'--','HandleVisibility','off')
    yline(0,'--','HandleVisibility','off')
    xlim([-100,100]); ylim(proj_lim)
    xlabel('Time from deceleration')

    yyaxis right
    current_vel_sac = mean(mean(data_ml1.vm_tot_sac(:,time_ind),3,'omitnan'),1);
    area(time_ind-250,current_vel_sac,'FaceColor',vm_color, ...
        'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off')
    ylim(vm_lim)
    legend('Box','off','Location','northwest')
    
    title(' p_y(t)')

    % =========================
    % (2) p_x(t): projection onto sac dir (0°)
    % =========================
    ax = nexttile(2);
    hold(ax,'on')

    proj_w = -cosd((0:45:360-45));

    % --- Bursters ---
    rate_pr_b = zeros(num_b,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_b = rate_pr_b + data_ss.rate_tot_sac(ind_b,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_b = rate_pr_b * num_b .* ss_cs_rho(ind_b) / sum(ss_cs_rho(ind_b));

    yyaxis left
    colororder({'k','k'})
    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_b,'omitnan'), ...
        std(rate_pr_b,'omitnan')/sqrt(num_b), 'r','alpha');
    hl.DisplayName = 'SS burst';
    hp.HandleVisibility = 'off';

    % --- Pausers ---
    rate_pr_p = zeros(num_p,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_p = rate_pr_p + data_ss.rate_tot_sac(ind_p,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_p = rate_pr_p * num_p .* ss_cs_rho(ind_p) / sum(ss_cs_rho(ind_p));

    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_p,'omitnan'), ...
        std(rate_pr_p,'omitnan')/sqrt(num_p), 'b','alpha');
    hl.HandleVisibility = 'off';

    % --- All SS ---
    rate_pr_ss = zeros(num_SS,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_ss = rate_pr_ss + data_ss.rate_tot_sac(:,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_ss = rate_pr_ss * num_SS .* ss_cs_rho / sum(ss_cs_rho);

    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_ss,'omitnan'), ...
        std(rate_pr_ss,'omitnan')/sqrt(num_SS), 'k','alpha');
    hl.DisplayName = 'SS all';
    hp.HandleVisibility = 'off';

    % --- MLI1 ---
    rate_pr_ml1 = zeros(num_MLI, numel(time_ind));
    for counter_dir = 1:8
        rate_pr_ml1 = rate_pr_ml1 + data_ml1.rate_tot_sac(:, time_ind, counter_dir) .* proj_w(counter_dir);
    end
    rate_pr_ml1 = rate_pr_ml1 .* ml1_cs_rho / sum(ml1_cs_rho) * num_MLI;
    [hl, hp] = boundedline(time_ind - 250, ...
        mean(rate_pr_ml1,'omitnan'), ...
        std(rate_pr_ml1,'omitnan') / sqrt(num_MLI), 'g', 'alpha');
    hl.DisplayName = ['MLI1 (n = ' num2str(num_MLI) ')'];
    hp.HandleVisibility = 'off';


    xline(0,'--','HandleVisibility','off')
    yline(0,'--','HandleVisibility','off')
    xlim([-100,100]); ylim(proj_lim)
    xlabel('Time from deceleration')

    yyaxis right
    current_vel_sac = mean(mean(data_ml1.vm_tot_sac(:,time_ind),3,'omitnan'),1);
    area(time_ind-250,current_vel_sac,'FaceColor',vm_color, ...
        'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off')
    ylim(vm_lim)
    ylabel('Velocity (deg/s)')
    title('p_x(t)')

    % =========================
    % (3) Phase plot: PC vs MLI1 (and optional MLI2) using p_x(t)
    % =========================
    ax = nexttile(3);
    hold(ax,'on')

    % Use the p_x(t) projections just computed above:
    sig_ml1 = mean(rate_pr_ml1,'omitnan');
    sig_ss  = mean(rate_pr_ss,'omitnan');
    sig_b   = mean(rate_pr_b,'omitnan');
    sig_p   = mean(rate_pr_p,'omitnan');

    plot(sig_ml1, sig_ss, 'k', 'DisplayName', 'SS all');
    plot(sig_ml1, sig_b,  'r', 'DisplayName', 'SS burst');
    plot(sig_ml1, sig_p,  'b', 'DisplayName', 'SS pause');
    if numel(sig_ml1) >= 101
        scatter(sig_ml1(101), sig_ss(101), 8, 'k', 'filled');
        scatter(sig_ml1(101), sig_b(101),  8, 'r', 'filled');
        scatter(sig_ml1(101), sig_p(101),  8, 'b', 'filled');
    end

    axis equal
    xlim(proj_lim); ylim(proj_lim)
    xline(0); yline(0)
    xlabel('MLI1 (Hz)'); ylabel('PC (Hz)')
    title('PC vs MLI1')
%     legend('Location','best','Box','off')

ESN_Beautify_Plot(fig,[12,5])
print(fig, fullfile(img_save_path,'ss-rate-population-Projection-tag6.pdf'), '-dpdf', '-bestfit');


end
