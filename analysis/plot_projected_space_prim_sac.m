function plot_projected_space_prim_sac(data_path)

 img_save_path = fullfile(data_path,'population_figs\ephys');
 save_path     = fullfile(data_path,'population_data');
 
JDM_params_funcs
data_ss  = load(fullfile(save_path,'SS_population_clique_prim_sac.mat'), 'data');
data_ml1 = load(fullfile(save_path,'MLI_population_clique_prim_sac.mat'), 'data');
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
%% plot SS MLI1 and MLI2 plots
vm_color = [204, 174, 98]/256;
vm_lim   = [0,600];
time_ind = (-100:150)+250;
proj_lim = [-30,55];

cond_labels = {'High', 'Low'};

% SS rate and vm of high
data_ss_high  = mean(data_ss.rate_tot_sac(:,:,:, [1 2]), 4, 'omitnan');
data_vm_high = mean(data_ss.vm_tot_sac(:,:,:, [1 2]), 4, 'omitnan');
% SS rate and vm of low
data_ss_low =mean(data_ss.rate_tot_sac(:,:,:, [3 4]), 4, 'omitnan');
data_vm_low = mean(data_ss.vm_tot_sac(:,:,:, [3 4]), 4, 'omitnan');
% MLI1 rate and vm
data_ml1_high  = mean(data_ml1.rate_tot_sac(:,:,:, [1 2]), 4, 'omitnan');
data_ml1_low  = mean(data_ml1.rate_tot_sac(:,:,:, [3 4]), 4, 'omitnan');

% One figure with 2x2 tiles: row 1 = high, row 2 = low
fig = figure;
t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
title(t,'SS projections-primary Sac.(cs-on+90° and cs-on°)')

for i_group = 1:2
    if i_group == 1
        current_data = data_ss_high;
        current_vm   = data_vm_high;
        cond_label   = cond_labels{1};
    else
        current_data = data_ss_low;
        current_vm   = data_vm_low;
        cond_label   = cond_labels{2};
    end

    % =========================
    % (1) Projection onto sac dir + 90  -> p_y(t)
    % =========================
    ax = nexttile((i_group-1)*2 + 1);  
    hold(ax,'on')
    colororder({'k','k'})
    proj_w = -cosd((0:45:360-45)+90);

    % Bursters
    rate_pr_b = zeros(num_b,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_b = rate_pr_b + current_data(ind_b,time_ind,counter_dir).*proj_w(counter_dir);
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

    % Pausers
    rate_pr_p = zeros(num_p,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_p = rate_pr_p + current_data(ind_p,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_p = rate_pr_p * num_p .* ss_cs_rho(ind_p) / sum(ss_cs_rho(ind_p));

    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_p,'omitnan'), ...
        std(rate_pr_p,'omitnan')/sqrt(num_p), 'b','alpha');
    hl.DisplayName = ['SS pause (n = ' num2str(num_p) ')'];
    hp.HandleVisibility = 'off';

    % All SS
    rate_pr_ss = zeros(num_SS,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_ss = rate_pr_ss + current_data(:,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_ss = rate_pr_ss * num_SS .* ss_cs_rho / sum(ss_cs_rho);

    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_ss,'omitnan'), ...
        std(rate_pr_ss,'omitnan')/sqrt(num_SS), 'k','alpha');
    hl.DisplayName = ['SS all (n = ' num2str(num_SS) ')'];
    hp.HandleVisibility = 'off';

    xline(0,'--','HandleVisibility','off')
    yline(0,'--','HandleVisibility','off')
    xlim([-100,100]); ylim(proj_lim)
    xlabel('Time from deceleration')

    yyaxis right
    current_vel_sac = mean(mean(current_vm(:,time_ind),3,'omitnan'),1);
    area(time_ind-250,current_vel_sac,'FaceColor',vm_color, ...
        'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off')
    ylim(vm_lim)

    if i_group == 1
        legend('Box','off','Location','northwest')
    end
    title([cond_label ' — p_y(t)'])

    % =========================
    % (2) Projection onto sac dir       -> p_x(t)
    % =========================
    ax = nexttile((i_group-1)*2 + 2);  
    hold(ax,'on')
    colororder({'k','k'})
    proj_w = -cosd((0:45:360-45));

    % Bursters
    rate_pr_b = zeros(num_b,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_b = rate_pr_b + current_data(ind_b,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_b = rate_pr_b * num_b .* ss_cs_rho(ind_b) / sum(ss_cs_rho(ind_b));

    yyaxis left
    colororder({'k','k'})
    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_b,'omitnan'), ...
        std(rate_pr_b,'omitnan')/sqrt(num_b), 'r','alpha');
    hl.DisplayName = 'SS burst';
    hp.HandleVisibility = 'off';

    % Pausers
    rate_pr_p = zeros(num_p,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_p = rate_pr_p + current_data(ind_p,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_p = rate_pr_p * num_p .* ss_cs_rho(ind_p) / sum(ss_cs_rho(ind_p));

    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_p,'omitnan'), ...
        std(rate_pr_p,'omitnan')/sqrt(num_p), 'b','alpha');
    hl.DisplayName = 'SS pause';
    hp.HandleVisibility = 'off';

    % All SS
    rate_pr_ss = zeros(num_SS,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_ss = rate_pr_ss + current_data(:,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_ss = rate_pr_ss * num_SS .* ss_cs_rho / sum(ss_cs_rho);

    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_ss,'omitnan'), ...
        std(rate_pr_ss,'omitnan')/sqrt(num_SS), 'k','alpha');
    hl.DisplayName = 'SS all';
    hp.HandleVisibility = 'off';

    xline(0,'--','HandleVisibility','off')
    yline(0,'--','HandleVisibility','off')
    xlim([-100,100]); ylim(proj_lim)
    xlabel('Time from deceleration')

    yyaxis right
    current_vel_sac = mean(mean(current_vm(:,time_ind),3,'omitnan'),1);
    area(time_ind-250,current_vel_sac,'FaceColor',vm_color, ...
        'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off')
    ylim(vm_lim)
    ylabel('Velocity (deg/s)')
    title([cond_label ' — p_x(t)'])
end

ESN_Beautify_Plot(fig,[10,8])
print(fig, fullfile(img_save_path,'ss-rate-purkinje-Projection_prim-sac.pdf'), '-dpdf', '-bestfit');

% ========= One figure with 2x3 tiles for SS / MLI1 (and optional MLI2) =========
fig = figure;
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
title(t,'SS & MLI1 projection- primary sac.(cs-on+90° and cs-on°)')

for i_group = 1:2
    % ---------- pick group ----------
    if i_group == 1
        cond_label   = cond_labels{1};
        current_data = data_ss_high;
        current_vm   = data_vm_high;
        current_ml1  = data_ml1_high;
        if exist('data_ml2_high','var'), current_ml2 = data_ml2_high; else, current_ml2 = []; end
    else
        cond_label   = cond_labels{2};
        current_data = data_ss_low;
        current_vm   = data_vm_low;
        current_ml1  = data_ml1_low;
        if exist('data_ml2_low','var'), current_ml2 = data_ml2_low; else, current_ml2 = []; end
    end

    % =========================
    % (1) p_y(t): projection onto sac dir + 90°
    % =========================
    ax = nexttile((i_group-1)*3 + 1);
    hold(ax,'on')
    colororder({'k','k'})
    proj_w = -cosd((0:45:360-45)+90);

    % --- Bursters ---
    rate_pr_b = zeros(num_b,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_b = rate_pr_b + current_data(ind_b,time_ind,counter_dir).*proj_w(counter_dir);
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
        rate_pr_p = rate_pr_p + current_data(ind_p,time_ind,counter_dir).*proj_w(counter_dir);
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
        rate_pr_ss = rate_pr_ss + current_data(:,time_ind,counter_dir).*proj_w(counter_dir);
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
        rate_pr_ml1 = rate_pr_ml1 + current_ml1(:, time_ind, counter_dir) .* proj_w(counter_dir);
    end
    rate_pr_ml1 = rate_pr_ml1 .* ml1_cs_rho / sum(ml1_cs_rho) * num_MLI;
    [hl, hp] = boundedline(time_ind - 250, ...
        mean(rate_pr_ml1,'omitnan'), ...
        std(rate_pr_ml1,'omitnan') / sqrt(num_MLI), 'g', 'alpha');
    hl.DisplayName = ['MLI1 (n = ' num2str(num_MLI) ')'];
    hp.HandleVisibility = 'off';

    % --- Optional MLI2 ---
    have_mli2 = ~isempty(current_ml2) && exist('num_MLI2','var') && exist('ml2_cs_rho','var') && ~isempty(ml2_cs_rho);
    if have_mli2
        rate_pr_ml2 = zeros(num_MLI2, numel(time_ind));
        for counter_dir = 1:8
            rate_pr_ml2 = rate_pr_ml2 + current_ml2(:, time_ind, counter_dir) .* proj_w(counter_dir);
        end
        rate_pr_ml2 = rate_pr_ml2 .* ml2_cs_rho / sum(ml2_cs_rho) * num_MLI2;
        [hl, hp] = boundedline(time_ind - 250, ...
            mean(rate_pr_ml2,'omitnan'), ...
            std(rate_pr_ml2,'omitnan') / sqrt(num_MLI2), 'm', 'alpha');
        hl.DisplayName = ['MLI2 (n = ' num2str(num_MLI2) ')'];
        hp.HandleVisibility = 'off';
    end

    xline(0,'--','HandleVisibility','off')
    yline(0,'--','HandleVisibility','off')
    xlim([-100,100]); ylim(proj_lim)
    xlabel('Time from deceleration')

    yyaxis right
    current_vel_sac = mean(mean(current_vm(:,time_ind),3,'omitnan'),1);
    area(time_ind-250,current_vel_sac,'FaceColor',vm_color, ...
        'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off')
    ylim(vm_lim)
    if i_group == 1
        legend('Box','off','Location','northwest')
    end
    title([cond_label ' — p_y(t)'])

    % =========================
    % (2) p_x(t): projection onto sac dir (0°)
    % =========================
    ax = nexttile((i_group-1)*3 + 2);
    hold(ax,'on')

    proj_w = -cosd((0:45:360-45));

    % --- Bursters ---
    rate_pr_b = zeros(num_b,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_b = rate_pr_b + current_data(ind_b,time_ind,counter_dir).*proj_w(counter_dir);
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
        rate_pr_p = rate_pr_p + current_data(ind_p,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_p = rate_pr_p * num_p .* ss_cs_rho(ind_p) / sum(ss_cs_rho(ind_p));

    [hl,hp] = boundedline(time_ind-250, ...
        mean(rate_pr_p,'omitnan'), ...
        std(rate_pr_p,'omitnan')/sqrt(num_p), 'b','alpha');
    hl.HandleVisibility = 'off';

    % --- All SS ---
    rate_pr_ss = zeros(num_SS,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_ss = rate_pr_ss + current_data(:,time_ind,counter_dir).*proj_w(counter_dir);
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
        rate_pr_ml1 = rate_pr_ml1 + current_ml1(:, time_ind, counter_dir) .* proj_w(counter_dir);
    end
    rate_pr_ml1 = rate_pr_ml1 .* ml1_cs_rho / sum(ml1_cs_rho) * num_MLI;
    [hl, hp] = boundedline(time_ind - 250, ...
        mean(rate_pr_ml1,'omitnan'), ...
        std(rate_pr_ml1,'omitnan') / sqrt(num_MLI), 'g', 'alpha');
    hl.DisplayName = ['MLI1 (n = ' num2str(num_MLI) ')'];
    hp.HandleVisibility = 'off';

    % --- Optional MLI2 ---
    if have_mli2
        rate_pr_ml2 = zeros(num_MLI2, numel(time_ind));
        for counter_dir = 1:8
            rate_pr_ml2 = rate_pr_ml2 + current_ml2(:, time_ind, counter_dir) .* proj_w(counter_dir);
        end
        rate_pr_ml2 = rate_pr_ml2 .* ml2_cs_rho / sum(ml2_cs_rho) * num_MLI2;
        [hl, hp] = boundedline(time_ind - 250, ...
            mean(rate_pr_ml2,'omitnan'), ...
            std(rate_pr_ml2,'omitnan') / sqrt(num_MLI2), 'm', 'alpha');
        hl.DisplayName = ['MLI2 (n = ' num2str(num_MLI2) ')'];
        hp.HandleVisibility = 'off';
    end

    xline(0,'--','HandleVisibility','off')
    yline(0,'--','HandleVisibility','off')
    xlim([-100,100]); ylim(proj_lim)
    xlabel('Time from deceleration')

    yyaxis right
    current_vel_sac = mean(mean(current_vm(:,time_ind),3,'omitnan'),1);
    area(time_ind-250,current_vel_sac,'FaceColor',vm_color, ...
        'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off')
    ylim(vm_lim)
    ylabel('Velocity (deg/s)')
    title([cond_label ' — p_x(t)'])

    % =========================
    % (3) Phase plot: PC vs MLI1 (and optional MLI2) using p_x(t)
    % =========================
    ax = nexttile((i_group-1)*3 + 3);
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

    if have_mli2
        sig_ml2 = mean(rate_pr_ml2,'omitnan');
        plot(sig_ml2, sig_ss, 'm--', 'DisplayName', 'SS all vs MLI2');
    end

    axis equal
    xlim(proj_lim); ylim(proj_lim)
    xline(0); yline(0)
    xlabel('MLI1 (Hz)'); ylabel('PC (Hz)')
    title([cond_label ' — PC vs MLI1'])
%     legend('Location','best','Box','off')
end

ESN_Beautify_Plot(fig,[12,8])
print(fig, fullfile(img_save_path,'ss-rate-population-Projection-prim-sac.pdf'), '-dpdf', '-bestfit');


end
