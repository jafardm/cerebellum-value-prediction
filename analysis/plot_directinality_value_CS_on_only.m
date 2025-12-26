function plot_directinality_value_CS_on_only(data_path)
% ------------------------------------------------------------
% make_figC4_value_CS_on_only
%   Plots Figure C4: Value (High–Low) and Direction (CS-on–CS+180)
%   using cs_on_rate_spe_rpe_prim_sac.mat
%   Data shape: [cell × time × direction × reward × epoch]
%   epochs: 1 = visual, 2 = saccade offset
% ------------------------------------------------------------

      %% === Load & basic setup primary saccade ===
    JDM_params_funcs
    S = load(fullfile(data_path,'population_data','cs_on_rate_spe_rpe_prim_sac'));
    num_cs = numel(S.cs_on_data);

    % --- Extract preferred directions (ρ, angle) ---
    cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, S.cs_on_data, 'UniformOutput', false));
    cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, S.cs_on_data, 'UniformOutput', false));
    cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, S.cs_on_data, 'UniformOutput', false));
    cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, S.cs_on_data, 'UniformOutput', false));

    % --- Choose stronger vector for each cell ---
    ind_m = cs_on_rho_sac > cs_on_rho_vis;
    cs_on_ang = cs_on_ang_vis/180*pi;
    cs_on_ang(ind_m) = cs_on_ang_sac(ind_m)/180*pi;

    ang_bins = params.sac.ang_values/180*pi;
    [~, cs_on_bin] = min(abs(angdiff(repmat(ang_bins, num_cs, 1), cs_on_ang*ones(size(ang_bins)))), [], 2);

    order_ang = nan(num_cs, 8);
    for i = 1:num_cs
        order_ang(i,:) = circshift(1:8, 5 - cs_on_bin(i));
    end
    order_ang = [order_ang, order_ang(:,1)]; % wrap for continuity

    %% === Rotate data into CS-on frame primary saccades ===
    cs_rate_rot = nan(num_cs, 500, 8, 6, 2);
    cs_vm_rot   = nan(num_cs, 500, 8, 6, 2);
    for i = 1:num_cs
        cs_rate_rot(i,:,:,:,:) = S.cs_rate(i,:,order_ang(i,1:8),:,1:2);
        cs_vm_rot(i,:,:,:,:)   = S.cs_vm(i,:, order_ang(i,1:8),:,1:2);
    end


    %% === Config constants ===
    dir_on       = 5;
    t_ind_vis    = (-50:200) + 250;
    t_ind_sac    = (-100:100) + 250;
    high_conds   = [1 2];
    low_conds    = [3 4];
    vm_color     = [204,174,98]/256;
    vm_lim       = [0 600];
    print_size   = [10 4];

    %% === Average across reward conditions ===
    cs_rate_h = squeeze(mean(cs_rate_rot(:,:,:, high_conds, :), 4, 'omitnan'));
    cs_rate_l = squeeze(mean(cs_rate_rot(:,:,:, low_conds,  :), 4, 'omitnan'));
    cs_vm_h   = squeeze(mean(cs_vm_rot(:,:,:,   high_conds, :), 4, 'omitnan'));

    %% === Compute VALUE (High–Low) ===
    % Visual epoch (1)
    val_vis     = squeeze(cs_rate_h(:, t_ind_vis, dir_on, 1) - cs_rate_l(:, t_ind_vis, dir_on, 1));
    val_vis180  = squeeze(cs_rate_h(:, t_ind_vis, 1,       1) - cs_rate_l(:, t_ind_vis, 1,       1));
    value_vis_mean  = mean(val_vis, 1, 'omitnan');
    value_vis_sem   = std(val_vis, [], 1, 'omitnan') / sqrt(num_cs);
    value_vis180_mean = mean(val_vis180, 1, 'omitnan');
    value_vis180_sem  = std(val_vis180, [], 1, 'omitnan') / sqrt(num_cs);

    % Saccade epoch (2)
    val_sac     = squeeze(cs_rate_h(:, t_ind_sac, dir_on, 2) - cs_rate_l(:, t_ind_sac, dir_on, 2));
    val_sac180  = squeeze(cs_rate_h(:, t_ind_sac, 1,       2) - cs_rate_l(:, t_ind_sac, 1,       2));
    value_sac_mean  = mean(val_sac, 1, 'omitnan');
    value_sac_sem   = std(val_sac, [], 1, 'omitnan') / sqrt(num_cs);
    value_sac180_mean = mean(val_sac180, 1, 'omitnan');
    value_sac180_sem  = std(val_sac180, [], 1, 'omitnan') / sqrt(num_cs);

    %% === Compute DIRECTION (CS-on – CS+180, avg across reward) ===
    cs_comb = (cs_rate_h + cs_rate_l) / 2;
    dir_vis = squeeze(cs_comb(:, t_ind_vis, dir_on, 1) - cs_comb(:, t_ind_vis, 1, 1));
    dir_sac = squeeze(cs_comb(:, t_ind_sac, dir_on, 2) - cs_comb(:, t_ind_sac, 1, 2));
    dir_vis_mean = mean(dir_vis, 1, 'omitnan');
    dir_vis_sem  = std(dir_vis, [], 1, 'omitnan') / sqrt(num_cs);
    dir_sac_mean = mean(dir_sac, 1, 'omitnan');
    dir_sac_sem  = std(dir_sac, [], 1, 'omitnan') / sqrt(num_cs);

    %% === Onset detection ===
    time_vis = t_ind_vis - 250;
    time_sac = t_ind_sac - 250;

    iv = find(value_vis_mean > 3.3 * value_vis_sem, 1, 'first');
    iv2 = find(dir_vis_mean > 3.3 * dir_vis_sem, 1, 'first');
    value_onset_vis_ms = ternary(~isempty(iv), time_vis(iv), NaN);
    dir_onset_vis_ms   = ternary(~isempty(iv2), time_vis(iv2), NaN);

    is_ = find(value_sac_mean > 3.3 * value_sac_sem, 1, 'first');
    is2 = find(dir_sac_mean > 3.3 * dir_sac_sem, 1, 'first');
    value_onset_sac_ms = ternary(~isempty(is_), time_sac(is_), NaN);
    dir_onset_sac_ms   = ternary(~isempty(is2), time_sac(is2), NaN);

    %% === Plot ===
    fig = figure;
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % ----- Target-aligned -----
    nexttile; hold on;
    [hv,pv] = boundedline(time_vis, value_vis_mean, value_vis_sem, 'r','alpha');
    [hv180,pv180] = boundedline(time_vis, value_vis180_mean, value_vis180_sem, 'm','alpha');
    [hd,pd] = boundedline(time_vis, dir_vis_mean, dir_vis_sem, 'b','alpha');
    pv.HandleVisibility='off'; pv180.HandleVisibility='off'; pd.HandleVisibility='off';
    xline(0,'--k');
    if ~isnan(value_onset_vis_ms), xline(value_onset_vis_ms,'--r'); end
    if ~isnan(dir_onset_vis_ms),   xline(dir_onset_vis_ms,'--b'); end
    xlabel('Time from target onset (ms)');
    ylabel('CS rate difference (Hz)');
    title('Value & Direction — target-aligned');
    legend([hv hv180 hd], {'Value (H−L, θ)', 'Value (H−L, θ+180)', 'Direction (θ−θ+180)'}, ...
           'Box','off','Location','northeast');

    yyaxis right;
    area(time_vis, mean(mean(cs_vm_h(:, t_ind_vis, :, 1),3,'omitnan'),1,'omitnan'), ...
        'FaceColor',vm_color,'FaceAlpha',.25,'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim);

    % ----- Saccade-aligned -----
    nexttile; hold on;
    [hv2,pv2] = boundedline(time_sac, value_sac_mean, value_sac_sem, 'r','alpha');
    [hv2180,pv2180] = boundedline(time_sac, value_sac180_mean, value_sac180_sem, 'm','alpha');
    [hd2,pd2] = boundedline(time_sac, dir_sac_mean, dir_sac_sem, 'b','alpha');
    pv2.HandleVisibility='off'; pv2180.HandleVisibility='off'; pd2.HandleVisibility='off';
    xline(0,'--k');
    xlabel('Time from saccade onset (ms)');
    ylabel('CS rate difference (Hz)');
    title('Value & Direction — saccade-aligned');
    legend([hv2 hv2180 hd2], {'Value (H−L, θ)', 'Value (H−L, θ+180)', 'Direction (θ−θ+180)'}, ...
           'Box','off','Location','northeast');

    yyaxis right;
    area(time_sac, mean(mean(cs_vm_h(:, t_ind_sac, :, 2),3,'omitnan'),1,'omitnan'), ...
        'FaceColor',vm_color,'FaceAlpha',.25,'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim);

    %% === Beautify & save ===
    ESN_Beautify_Plot(fig, print_size);
    print(fig, fullfile(data_path,'population_figs','ephys','Fig1-E.pdf'), '-dpdf', '-bestfit');

    %% === Console reporting ===
    fprintf('----------------------------------------\n');
    fprintf('Value onset (target, θ):      %.1f ms\n', value_onset_vis_ms);
    fprintf('Direction onset (target):     %.1f ms\n', dir_onset_vis_ms);
    fprintf('Value onset (saccade, θ):     %.1f ms\n', value_onset_sac_ms);
    fprintf('Direction onset (saccade):    %.1f ms\n', dir_onset_sac_ms);
    fprintf('----------------------------------------\n');



    %% === Load corrective saccade dataset ===
    S2 = load(fullfile(data_path,'population_data','cs_on_rate_spe_rpe_corr_sac'));
    num_cs = numel(S2.cs_on_data);
    
    % --- Compute CS-on alignment ---
    cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg,S2.cs_on_data,'UniformOutput',false));
    cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg,S2.cs_on_data,'UniformOutput',false));
    cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg,S2.cs_on_data,'UniformOutput',false));
    cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg,S2.cs_on_data,'UniformOutput',false));
    
    ind_m = cs_on_rho_sac > cs_on_rho_vis;
    cs_on_ang = cs_on_ang_vis/180*pi;
    cs_on_ang(ind_m) = cs_on_ang_sac(ind_m)/180*pi;
    
    ang_bins = params.sac.ang_values/180*pi;
    [~, cs_on_bin] = min(abs(angdiff(repmat(ang_bins,num_cs,1), cs_on_ang*ones(size(ang_bins)))),[],2);
    
    order_ang = nan(num_cs,8);
    for i = 1:num_cs
        order_ang(i,:) = circshift(1:8, 5 - cs_on_bin(i));
    end
    order_ang = [order_ang, order_ang(:,1)];
    
    %% === Rotate data into CS-on frame ===
    cs_rate_rot_corr = nan(num_cs,500,8,4,2);
    cs_vm_rot_corr   = nan(num_cs,500,8,4,2);
    for i = 1:num_cs
        cs_rate_rot_corr(i,:,:,:,:) = S2.cs_rate(i,:,order_ang(i,1:8),:,:);
        cs_vm_rot_corr(i,:,:,:,:)   = S2.cs_vm(i,:, order_ang(i,1:8),:,:);
    end
    
    %% === Config ===
    dir_on     = 5;
    t_ind_vis  = (-50:200) + 250;
    time_vis   = t_ind_vis - 250;
    num_cs     = size(cs_rate_rot_corr,1);
    vm_color   = [204,174,98]/256;
    vm_lim     = [0 400];
    print_size = [12 4];
    
    % Reward indices: 1=HH, 2=HL, 3=LL, 4=LH
    
    %% === 1️⃣ SPE Response: Directionality averaged across all reward conditions ===
    dir_all_vis = squeeze(mean(cs_rate_rot_corr(:,t_ind_vis,dir_on,:,1) - cs_rate_rot_corr(:,t_ind_vis,1,:,1),4,'omitnan'));
    dir_all_mean = mean(dir_all_vis,1,'omitnan');
    dir_all_sem  = std(dir_all_vis,[],1,'omitnan')/sqrt(num_cs);
    
    %% === 2️⃣ Value Response: (High–Low) ignoring RPE ===
    % High = [HH, LH] (1,4); Low = [HL, LL] (2,3)
    high_conds = [1 4];
    low_conds  = [2 3];
    cs_rate_h  = squeeze(mean(cs_rate_rot_corr(:,:,:,high_conds,1),4,'omitnan'));
    cs_rate_l  = squeeze(mean(cs_rate_rot_corr(:,:,:,low_conds,1),4,'omitnan'));
    
    val_vis = squeeze(cs_rate_h(:,t_ind_vis,dir_on) - cs_rate_l(:,t_ind_vis,dir_on));
    val_mean = mean(val_vis,1,'omitnan');
    val_sem  = std(val_vis,[],1,'omitnan')/sqrt(num_cs);
    
    %% === 3️⃣ RPE Responses: (LH–HH) and (HL–LL) ===
    rpe_pos_vis = squeeze(cs_rate_rot_corr(:,t_ind_vis,dir_on,4,1) - cs_rate_rot_corr(:,t_ind_vis,dir_on,1,1)); % LH−HH
    rpe_neg_vis = squeeze(cs_rate_rot_corr(:,t_ind_vis,dir_on,2,1) - cs_rate_rot_corr(:,t_ind_vis,dir_on,3,1)); % HL−LL
    
    rpe_pos_mean = mean(rpe_pos_vis,1,'omitnan');
    rpe_pos_sem  = std(rpe_pos_vis,[],1,'omitnan')/sqrt(num_cs);
    rpe_neg_mean = mean(rpe_neg_vis,1,'omitnan');
    rpe_neg_sem  = std(rpe_neg_vis,[],1,'omitnan')/sqrt(num_cs);
    
    %% === 4️⃣ Onset (latency) detection ===
    lat_thresh = 3.3;  % detection threshold
    baseline_mask = time_vis < 0;
    
    % --- Baseline correction (remove mean baseline offset) ---
    dir_base = mean(dir_all_mean(baseline_mask), 'omitnan');
    val_base = mean(val_mean(baseline_mask), 'omitnan');
    dir_all_mean_bc = dir_all_mean - dir_base;
    val_mean_bc     = val_mean - val_base;
    
    % --- Only search after stimulus onset ---
    search_mask = time_vis >= 0;
    t_search = time_vis(search_mask);
    
    dir_mean_seg = dir_all_mean_bc(search_mask);
    dir_sem_seg  = dir_all_sem(search_mask);
    val_mean_seg = val_mean_bc(search_mask);
    val_sem_seg  = val_sem(search_mask);
    cons_win = 5;
    dir_logic = dir_mean_seg > lat_thresh * dir_sem_seg;
    val_logic = val_mean_seg > lat_thresh * val_sem_seg;
    i_dir = find(conv(double(dir_logic), ones(1,cons_win), 'same') >= cons_win, 1, 'first');
    i_val = find(conv(double(val_logic), ones(1,cons_win), 'same') >= cons_win, 1, 'first');
    
    dir_onset_vis_ms = ternary(~isempty(i_dir), t_search(i_dir), NaN);
    val_onset_vis_ms = ternary(~isempty(i_val), t_search(i_val), NaN);
    
    %% === Plot ===
    fig2 = figure;
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    
    % ----- (1) SPE: Directionality averaged across all reward conditions -----
    nexttile; hold on;
    [h1,p1] = boundedline(time_vis, dir_all_mean, dir_all_sem, 'b','alpha');
    xline(0,'--k');
    xline(dir_onset_vis_ms,'--b','LineWidth',1.2);
    xlabel('Time from target onset (ms)');
    ylabel('CS rate difference (Hz)');
    title('Direction Encoding');
    yyaxis right;
    area(time_vis, mean(mean(cs_vm_rot_corr(:,t_ind_vis,:,1,1),3,'omitnan'),1,'omitnan'), ...
        'FaceColor',vm_color,'FaceAlpha',.25,'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim);
    
    % ----- (2) VALUE: H−L ignoring RPE -----
    nexttile; hold on;
    [h2,p2] = boundedline(time_vis, val_mean, val_sem, 'r','alpha');
    xline(0,'--k');
    yline(0,'--k');
    xline(val_onset_vis_ms,'--r','LineWidth',1.2);
    xlabel('Time from target onset (ms)');
    title('Value Response — (LH HH) - (LL,HL) ');
    yyaxis right;
    area(time_vis, mean(mean(cs_vm_rot_corr(:,t_ind_vis,:,1,1),3,'omitnan'),1,'omitnan'), ...
        'FaceColor',vm_color,'FaceAlpha',.25,'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim);
    
    % ----- (3) RPE: LH−HH and HL−LL -----
    nexttile; hold on;
    [h3,p3] = boundedline(time_vis, rpe_pos_mean, rpe_pos_sem, 'm','alpha'); % RPE+
    [h4,p4] = boundedline(time_vis, rpe_neg_mean, rpe_neg_sem, 'c','alpha'); % RPE−
    xline(0,'--k');
    yline(0,'--k');
    xlabel('Time from target onset (ms)');
    title('LH−HH (RPE⁺) and HL−LL (RPE⁻)');
    legend([h3 h4],{'RPE⁺ (LH−HH)','RPE⁻ (HL−LL)'},'Box','off','Location','northeast');
    yyaxis right;
    area(time_vis, mean(mean(cs_vm_rot_corr(:,t_ind_vis,:,1,1),3,'omitnan'),1,'omitnan'), ...
        'FaceColor',vm_color,'FaceAlpha',.25,'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim);
    
    %% Save
    ESN_Beautify_Plot(fig2, print_size);
    print(fig2, fullfile(data_path,'population_figs','ephys','directionality-reward-corr-sacs.pdf'), '-dpdf', '-bestfit');
    fprintf('Saved: Fig_SPE_Value_RPE_corr.pdf — Direction (SPE), Value, and RPE (corrective saccades)\n');

    
end


%% --- Helper inline function ---
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
%% --- Helper function ---
function [m,s] = mean_sem(x, n)
    m = mean(x,1,'omitnan');
    s = std(x,[],1,'omitnan') / sqrt(n);
end
