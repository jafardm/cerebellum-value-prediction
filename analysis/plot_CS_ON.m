function plot_CS_ON(data_path)

   JDM_params_funcs
   load(fullfile(data_path,'population_data','cs_on_rate_spe_rpe_prim_sac'));
    % find the visual and motor rho for each cell
    cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg,cs_on_data,...
        'UniformOutput',false));
    cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg,cs_on_data,...
        'UniformOutput',false));
    
    cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg,cs_on_data,...
        'UniformOutput',false));
    cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg,cs_on_data,...
        'UniformOutput',false));
    
    num_cs = numel(cs_on_data);
    
    % select larger vector for each cell
    ind_m = cs_on_rho_sac > cs_on_rho_vis;
    
    cs_on_ang = cs_on_ang_vis/180*pi;
    cs_on_ang(ind_m) = cs_on_ang_sac(ind_m)/180*pi;
    
    cs_on_rho = cs_on_rho_vis;
    cs_on_rho(ind_m) = cs_on_rho_sac(ind_m);
    
    ang_bins = params.sac.ang_values/180*pi;
    
    [~,cs_on_bin] = min(abs(angdiff(repmat(ang_bins,num_cs,1),...
        cs_on_ang*ones(size(ang_bins)))),[],2);
    
    order_ang = nan(num_cs,8);
    
    for counter_cs = 1:num_cs
        order_ang(counter_cs,:) = circshift(1:8,5-cs_on_bin(counter_cs));
    end
    order_ang          = [order_ang,order_ang(:,1)];


    cs_on_rate_vis = cell2mat(cellfun(@(x) x.vis.fr_avg,...
        cs_on_data,'UniformOutput',false));
    cs_on_rate_sac = cell2mat(cellfun(@(x) x.sac.fr_avg,...
        cs_on_data,'UniformOutput',false));
    
    % organizing the tuning of the cs rates
    cs_on_rate_vis_rot = nan(num_cs,9);
    cs_on_rate_sac_rot = nan(num_cs,9);
    cs_rate_rot        = nan(num_cs,500,8,6,2);
    cs_vm_rot          = nan(num_cs,500,8,6,2);
 
    for counter_cs = 1:num_cs
        cs_on_rate_vis_rot(counter_cs,:) =...
            cs_on_rate_vis(counter_cs,order_ang(counter_cs,:));
        cs_on_rate_sac_rot(counter_cs,:) =...
            cs_on_rate_sac(counter_cs,order_ang(counter_cs,:));
        cs_rate_rot(counter_cs,:,:,:,:)    =...
            cs_rate(counter_cs,:,order_ang(counter_cs,1:8),:,1:2); % picked visual and sac dim
        cs_vm_rot(counter_cs,:,:,:,:)    =...
            cs_vm(counter_cs,:,order_ang(counter_cs,1:8),:,1:2);
    end

%% parameters
    leg_ = {'cs-on','cs±45','cs±90','cs±135','cs±180'};
    colors_h = cool(5); 
    colors_ = winter(5);
    colors_l = colors_(end:-1:1,:);  
    dirs = 1:5; 
    vm_color = [204, 174, 98]/256;
    vm_lim = [0 600];
    t_ind_vis = (-50:200) + 250;
    t_ind_sac = (-100:150) + 250;

    % Time window used for averaging and scatter plots
    t = -249:250;
    time_window_vis = [50 150];
    time_window_sac = [-50 50];
    dir = 5; % cs-on direction

%% CS rate for High and low conditions
    high_conds = [1 2 ];  % e.g., High-High, High-Low forced
    low_conds  = [3 4 ];  % e.g., Low-Low, Low-High forced
    
    cs_rate_rot_h = squeeze(mean(cs_rate_rot(:,:,:, high_conds, :), 4, 'omitnan'));
    cs_rate_rot_l = squeeze(mean(cs_rate_rot(:,:,:, low_conds, :), 4, 'omitnan'));
    
    cs_vm_rot_h = squeeze(mean(cs_vm_rot(:,:,:, high_conds, :), 4, 'omitnan'));
    cs_vm_rot_l = squeeze(mean(cs_vm_rot(:,:,:, low_conds, :), 4, 'omitnan'));

    %% Classify cells based on reward response
% Find indices for 50–150 ms
    % Find indices for the time windows
    t_inds_vis = find(t >= time_window_vis(1) & t <= time_window_vis(2));
    t_inds_sac = find(t >= time_window_sac(1) & t <= time_window_sac(2));

    % Average CS rates over time windows
    cs_high_vis = squeeze(cs_rate_rot_h(:,t_inds_vis,dir,1));
    cs_low_vis  = squeeze(cs_rate_rot_l(:,t_inds_vis,dir,1));
    cs_high_mean_vis = mean(cs_high_vis, 2, 'omitnan');
    cs_low_mean_vis  = mean(cs_low_vis, 2, 'omitnan');

    cs_high_sac = squeeze(cs_rate_rot_h(:,t_inds_sac,dir,2));
    cs_low_sac  = squeeze(cs_rate_rot_l(:,t_inds_sac,dir,2));
    cs_high_mean_sac = mean(cs_high_sac, 2, 'omitnan');
    cs_low_mean_sac  = mean(cs_low_sac, 2, 'omitnan');
    cs_high_tot = mean([cs_high_mean_vis cs_high_mean_sac],2);
    cs_low_tot = mean([cs_low_mean_vis cs_low_mean_sac],2);
    % Count how many cells are above the unity line
    above_unity = sum(cs_high_tot > cs_low_tot);
    total_cells = numel(cs_high_tot);
    
    fprintf('Cells above unity line: %d / %d\n', above_unity, total_cells);

    % Difference
    reward_diff = cs_high_tot - cs_low_tot;

    reward_positive_idx = reward_diff > 0;
    reward_negative_idx = reward_diff < 0;
    reward_positive_rho = mean(cs_on_rho_sac(reward_positive_idx));
    reward_negative_rho = mean(cs_on_rho_sac(reward_negative_idx));
%% Heatmap
    cs_rate_high = squeeze(cs_rate_rot_h(:,t_ind_vis,dir,1));  % [cells × time]
    cs_rate_low  = squeeze(cs_rate_rot_l(:,t_ind_vis,dir,1));  % same

    % Z-score each neuron's activity (across time) independently
    cs_rate_high_z = (cs_rate_high - mean(cs_rate_high, 2, 'omitnan')) ./ ...
                     std(cs_rate_high, 0, 2, 'omitnan');
    
    cs_rate_low_z  = (cs_rate_low - mean(cs_rate_low, 2, 'omitnan')) ./ ...
                     std(cs_rate_low, 0, 2, 'omitnan');

    % Sort neurons by time of peak response in high
    rng(42); % seed for reproducibility
    rand_sort_idx = randperm(size(cs_rate_high_z, 1));
    cs_rate_high_sorted = cs_rate_high_z(rand_sort_idx, :);
    cs_rate_low_sorted  = cs_rate_low_z(rand_sort_idx, :);  

%%  Difference: High − Low, averaged across directions
    % within cell differnce
    num_cells = size(cs_rate_rot_h,1);
    % Average over directions
    cs_high_avg_dir = squeeze(mean(cs_rate_rot_h(:,t_ind_vis,:,1), 3, 'omitnan'));  % [cells × time]
    cs_low_avg_dir  = squeeze(mean(cs_rate_rot_l(:,t_ind_vis,:,1), 3, 'omitnan'));  % same
    
    reward_diff = cs_high_avg_dir - cs_low_avg_dir;  % [cells × time]
    reward_diff_mean = mean(reward_diff, 1, 'omitnan');
    reward_diff_sem  = std(reward_diff, [], 1, 'omitnan') / sqrt(num_cells);

%% direction Difference: CS-on − CS+180, averaged across reward   

    % Combine high + low across reward
    cs_combined = (cs_rate_rot_h + cs_rate_rot_l) / 2;
    
    % Extract CS-on (dir=5) and CS+180 (dir=1)
    cs_on   = squeeze(cs_combined(:,t_ind_vis,5,1));   % [cells × time]
    cs_180  = squeeze(cs_combined(:,t_ind_vis,1,1));   % [cells × time]
    
    dir_diff = cs_on - cs_180;
    dir_diff_mean = mean(dir_diff, 1, 'omitnan');
    dir_diff_sem  = std(dir_diff, [], 1, 'omitnan') / sqrt(num_cells);
%% Stats of reward_diff and direction Difference

time_ms = t_ind_vis -250;  

% Define threshold for 99% confidence (mean > 3 * SEM)
reward_thresh = 3.3 * reward_diff_sem;
dir_thresh    = 3.3 * dir_diff_sem;

% Find first time point where mean exceeds threshold
reward_sig_idx = find(reward_diff_mean > reward_thresh, 1, 'first');
dir_sig_idx    = find(dir_diff_mean > dir_thresh, 1, 'first');

% Convert index to time
reward_onset_ms = time_ms(reward_sig_idx);
dir_onset_ms    = time_ms(dir_sig_idx);

% Display result
fprintf('Reward encoding begins at %d ms after target onset.\n', reward_onset_ms);
fprintf('Direction encoding begins at %d ms after target onset.\n', dir_onset_ms);

%%  reward-pos and reward neg cells
    
    % Preallocate
    n_pos = sum(reward_positive_idx);
    n_neg = sum(reward_negative_idx);
    cs_rate_rot_p = nan(n_pos, 500, 8, 2);  
    cs_vm_rot_p   = nan(n_pos, 500, 8, 2);
    cs_rate_rot_n = nan(n_neg, 500, 8, 2);
    cs_vm_rot_n   = nan(n_neg, 500, 8, 2);
    
    % Positive cells
    pos_ids = find(reward_positive_idx);
    for i = 1:n_pos
        cs_id = pos_ids(i);
        cs_rate_rot_p(i,:,:,:) = mean(cs_rate(cs_id,:,order_ang(cs_id,1:8),1:4,1:2), 4, 'omitnan');
        cs_vm_rot_p(i,:,:,:)   = mean(cs_vm(cs_id,:,order_ang(cs_id,1:8),1:4,1:2),   4, 'omitnan');
    end
    
    % Negative cells
    neg_ids = find(reward_negative_idx);
    for i = 1:n_neg
        cs_id = neg_ids(i);
        cs_rate_rot_n(i,:,:,:) = mean(cs_rate(cs_id,:,order_ang(cs_id,1:8),1:4,1:2), 4, 'omitnan');
        cs_vm_rot_n(i,:,:,:)   = mean(cs_vm(cs_id,:,order_ang(cs_id,1:8),1:4,1:2),   4, 'omitnan');
    end
%% Heatmap sorted on reward-poitive and rew-neagtive
% Get CS-on direction (dir=5), target-aligned
cs_rate_pos = squeeze(cs_rate_rot_p(:, t_ind_vis, dir, 1));  % [n_pos × time]
cs_rate_neg = squeeze(cs_rate_rot_n(:, t_ind_vis, dir, 1));  % [n_neg × time]

% Z-score across time for each neuron
cs_rate_pos_z = (cs_rate_pos - mean(cs_rate_pos, 2, 'omitnan')) ./ ...
                 std(cs_rate_pos, 0, 2, 'omitnan');
cs_rate_neg_z = (cs_rate_neg - mean(cs_rate_neg, 2, 'omitnan')) ./ ...
                 std(cs_rate_neg, 0, 2, 'omitnan');

% Concatenate: reward-positive cells first, reward-negative below
cs_rate_sorted = [cs_rate_pos_z; cs_rate_neg_z];
n_pos = size(cs_rate_pos_z, 1);
n_neg = size(cs_rate_neg_z, 1);

%% plot CS on across all cells

cs_data_tot_cells = squeeze(mean(cs_rate_rot,4,'omitnan'));
vm_data_tot_cells = squeeze(mean(cs_vm_rot,4,'omitnan'));

fig1 = figure;
colors_t = cool(5);
for counter_dir = numel(dirs):-1:1

    current_dir = dirs(counter_dir);

    if current_dir == 1 || current_dir == 5
        current_dir_ = current_dir;
    else
        current_dir_ = [current_dir, 10-current_dir];
    end
    current_cs_rate_sac = squeeze(mean(cs_data_tot_cells(:,t_ind_sac,current_dir_,2), 3,'omitnan'));
    current_cs_rate_vis = squeeze(mean(cs_data_tot_cells(:,t_ind_vis,current_dir_,1), 3,'omitnan'));

    if counter_dir == 1
        current_cs_vm_sac = squeeze(mean(vm_data_tot_cells(:,t_ind_sac,:,2), 3 ,'omitnan'));
        current_cs_vm_vis = squeeze(mean(vm_data_tot_cells(:,t_ind_vis,:,1),  3,'omitnan'));
    end

    subplot(1,2,1);
    title(sprintf('CS Cells (n = %d)', num_cs))
    colororder({'k','k'})
    yyaxis left
    [hl1,hl2] = boundedline(t_ind_vis - 250,...
        squeeze(mean(current_cs_rate_vis,1)),...
        squeeze(std(current_cs_rate_vis,1))/sqrt(num_cs),...
        'Color',colors_t(counter_dir,:),'alpha');
    hl1.DisplayName = leg_{6 - counter_dir};
    hl2.HandleVisibility = 'off';
    hold on;
    xlabel('time from target onset (ms)')
    legend('Box','off');
    xline(0,'--','HandleVisibility','off');
    
    if counter_dir == 1
        yyaxis right
        area(t_ind_vis-250,squeeze(mean(current_cs_vm_vis,1)),...
            'FaceColor',vm_color,...
            'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off')
        ylim(vm_lim);
    end

    subplot(1,2,2);
    title(sprintf('CS Cells (n = %d)', num_cs))
    colororder({'k','k'})
    yyaxis left
    [hl1,hl2] = boundedline(t_ind_sac-250,...
        squeeze(mean(current_cs_rate_sac,1)),...
        squeeze(std(current_cs_rate_sac,1))/sqrt(num_cs),...
        'Color',colors_t(counter_dir,:),'alpha');
    hl1.DisplayName = leg_{6-counter_dir};
    hl2.HandleVisibility = 'off';
    hold on
    xline(0,'--','HandleVisibility','off');
    xlabel('time from saccade onset (ms)')
    
    if counter_dir == 1
        yyaxis right
        area(t_ind_sac-250,squeeze(mean(current_cs_vm_sac,1)),...
            'FaceColor',vm_color,...
            'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off')
        ylim(vm_lim);
        ylabel('Velocity (deg/sec)');
    end
end

ESN_Beautify_Plot(fig1,[10,4])
    
print(fig1, fullfile(data_path,'population_figs\ephys',...
    'Fig1-C1.pdf'), '-dpdf', '-bestfit');

%% plot high low relitive to cs-on direction
    fig2 = figure;
    for counter_dir = numel(dirs):-1:1
        current_dir = dirs(counter_dir);
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            current_dir_ = [current_dir, 10-current_dir];  % e.g., [2 8] for ±45
        end
    
        % === High data ===
        cs_rate_vis_h = squeeze(mean(cs_rate_rot_h(:,t_ind_vis,current_dir_,1), 3));
        cs_rate_sac_h = squeeze(mean(cs_rate_rot_h(:,t_ind_sac,current_dir_,2), 3));
    
        if counter_dir == 1
            cs_vm_vis_h = squeeze(mean(cs_vm_rot_h(:,t_ind_vis,:,1), 3));
            cs_vm_sac_h = squeeze(mean(cs_vm_rot_h(:,t_ind_sac,:,2), 3));
            cs_vm_vis_l = squeeze(mean(cs_vm_rot_l(:,t_ind_vis,:,1), 3));
            cs_vm_sac_l = squeeze(mean(cs_vm_rot_l(:,t_ind_sac,:,2), 3));

        end
    
        % === Low data ===
        cs_rate_vis_l = squeeze(mean(cs_rate_rot_l(:,t_ind_vis,current_dir_,1), 3));
        cs_rate_sac_l = squeeze(mean(cs_rate_rot_l(:,t_ind_sac,current_dir_,2), 3));
      
        % VISUAL PANEL
        subplot(1,2,1); 
        title(sprintf('CS Cells (n = %d)', num_cs))

        % High
        colororder({'k','k'}); 
        yyaxis left;
        [hl1, hp1] = boundedline(t_ind_vis - 250,...
            mean(cs_rate_vis_h, 1, 'omitnan'),...
            std(cs_rate_vis_h, [], 1, 'omitnan') / sqrt(num_cs),'Color',colors_h(counter_dir,:),'alpha');
       
        hl1.DisplayName = ['High - ' leg_{6 - counter_dir}];
        hp1.HandleVisibility = 'off';
        % Low
        [hl2, hp2] = boundedline(t_ind_vis - 250,...
            mean(cs_rate_vis_l, 1, 'omitnan'),...
            std(cs_rate_vis_l, [], 1, 'omitnan') / sqrt(num_cs), 'Color', colors_h(counter_dir,:),'alpha');
        set(hl2, 'LineStyle', '--');
        hl2.DisplayName = ['Low - ' leg_{6 - counter_dir}];
        hp2.HandleVisibility = 'off';
    
        hold on; xlabel('time from stimulus onset (ms)');
        ylabel('CS rate(Hz)')
        ylabel('complex spike rate (Hz)');
        xline(0,'--','HandleVisibility','off');
    
       if counter_dir == 1
            yyaxis right;
            % --- High velocity ---
            area(t_ind_vis - 250, mean(cs_vm_vis_h, 1, 'omitnan'),...
                'FaceColor', vm_color, 'FaceAlpha', 0.3,...
                'EdgeColor', 'none', 'DisplayName','Vel High');
            hold on;
        
            % --- Low velocity ---
            area(t_ind_vis - 250, mean(cs_vm_vis_l, 1, 'omitnan'),...
                'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.3,...
                'EdgeColor', 'none', 'DisplayName','Vel Low');
        
            ylim(vm_lim);
       end

        % SACCADE PANEL
        subplot(1,2,2);  
        title(sprintf('CS Cells (n = %d)', num_cs))
        % High
        colororder({'k','k'});
        yyaxis left;
        [hl3, hp3] = boundedline(t_ind_sac - 250,...
            mean(cs_rate_sac_h, 1, 'omitnan'),...
            std(cs_rate_sac_h, [], 1, 'omitnan') / sqrt(num_cs),...
             'Color', colors_h(counter_dir,:),'alpha');
    
        hl3.DisplayName = ['High - ' leg_{6-counter_dir}];
        hp3.HandleVisibility = 'off';
    
        % Low
        [hl4, hp4] = boundedline(t_ind_sac - 250,...
            mean(cs_rate_sac_l, 1, 'omitnan'),...
            std(cs_rate_sac_l, [], 1, 'omitnan') / sqrt(num_cs),'Color', colors_h(counter_dir,:),'alpha');
        set(hl4, 'LineStyle', '--');
    
        hl4.DisplayName = ['Low - ' leg_{6-counter_dir}];
        hp4.HandleVisibility = 'off';
    
        hold on; xlabel('time from saccade onset (ms)');
        xline(0,'--','HandleVisibility','off');

       if counter_dir == 1
        yyaxis right;
        % --- High velocity ---
        area(t_ind_sac - 250, mean(cs_vm_sac_h, 1, 'omitnan'),...
            'FaceColor', vm_color, 'FaceAlpha', 0.3,...
            'EdgeColor', 'none', 'DisplayName','Vel High');
        hold on;
    
        % --- Low velocity ---
        area(t_ind_sac - 250, mean(cs_vm_sac_l, 1, 'omitnan'),...
            'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.3,...
            'EdgeColor', 'none', 'DisplayName','Vel Low');
    
        ylim(vm_lim);
        ylabel('Velocity (deg/sec)')
       end


        legend('Box','off', 'Location','northeast');
     end

    ESN_Beautify_Plot(fig2, [10, 6]);
    print(fig2, fullfile(data_path,'population_figs\ephys',...
    'Fig1-C2.pdf'), '-dpdf', '-bestfit');
%% plot High/Low heatmap
    fig3 = figure;
    % High value heatmap
    subplot(1,2,1);
    imagesc(t_ind_vis-250, 1:num_cs, cs_rate_high_sorted);
    xlabel('Time from target onset (ms)');
    ylabel('Neurons (sorted by peak time)');
    title('High Value — CS-on Direction');
    colormap(params.rwb_map)
    caxis([0 4]);
    s = colorbar;
    xline(0, '--w', 'LineWidth', 1.2);
    % Low value heatmap (same neuron order)
    subplot(1,2,2);
    imagesc(t_ind_vis-250, 1:num_cs, cs_rate_low_sorted);
    xlabel('Time from target onset (ms)');
    ylabel('Neurons (same order)');
    title('Low Value — CS-on Direction');
    caxis([0 4]);
    s = colorbar;
    xline(0, '--w', 'LineWidth', 1.2);
    ESN_Beautify_Plot(fig3, [10, 4]); 
    print(fig3, fullfile(data_path,'population_figs\ephys',...
        'Fig1-C3.pdf'), '-dpdf', '-bestfit');

%% plot CS-on − CS+180, averaged across reward
    fig4 = figure;
    hold on
    
    % Plot mean ± SEM for reward
    [hl1, hp1] = boundedline(t_ind_vis - 250, reward_diff_mean, reward_diff_sem, 'r','alpha');
    hl1.DisplayName = 'Reward (High - Low)';
    hp1.HandleVisibility = 'off';
    
    % Plot mean ± SEM for direction
    [hl2, hp2] = boundedline(t_ind_vis - 250, dir_diff_mean, dir_diff_sem, 'b','alpha');
    hl2.DisplayName = 'Direction (CS-on - CS+180)';
    hp2.HandleVisibility = 'off';
    
    % Vertical lines
    x0 = xline(0, '--k', 'DisplayName', 'Target Onset');
    x1 = xline(reward_onset_ms, '--r', 'DisplayName', 'Reward Onset (99.9% CI)');
    x2 = xline(dir_onset_ms, '--b', 'DisplayName', 'Direction Onset (99.9% CI)');
    
    % Axis labels and title
    xlabel('Time from target onset (ms)');
    ylabel('CS rate difference (mean ± SEM)');
    title('Time Course of Reward & Direction Encoding');
    legend('Box', 'off', 'Location', 'northeast');

    ESN_Beautify_Plot(fig4, [10 4]); 
    print(fig4, fullfile(data_path,'population_figs\ephys',...
      'Fig1-C4.pdf'), '-dpdf', '-bestfit');
    
%% Plot reward positive and negative cells    
    fig5 = figure;
    subplot(1,2,1)
    scatter(cs_low_tot, cs_high_tot, 25, 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.2);
    hold on;
    
    % Unity line
    max_val = max([cs_low_tot; cs_high_tot], [], 'omitnan');
    plot([0 max_val], [0 max_val], 'k--', 'LineWidth', 1.5);
    
    xlabel('Low reward CS rate');
    ylabel('High reward CS rate');
    title(sprintf('Reward Modulation of CS (N = %d, above unity = %d)', total_cells, above_unity));
    axis equal; xlim([0 max_val]); ylim([0 max_val]);

    subplot(1,2,2)
    % Plot reward-positive cells
    scatter(cs_low_tot(reward_positive_idx), cs_high_tot(reward_positive_idx), ...
        30, [0.2 0.6 1], 'filled', 'MarkerFaceAlpha', 0.6);  % blue
    hold on;
    
    % Plot reward-negative cells
    scatter(cs_low_tot(reward_negative_idx), cs_high_tot(reward_negative_idx), ...
        30, [1 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);  % red
    
    % Unity line
    max_val = max([cs_low_tot; cs_high_tot], [], 'omitnan');
    plot([0 max_val], [0 max_val], 'k--', 'LineWidth', 1.5);
    
    xlabel('Low reward CS rate');
    ylabel('High reward CS rate');
    legend({'Reward-Positive','Reward-Negative'}, 'Location', 'southeast', 'Box', 'off');
    title(sprintf('Reward Modulation (Pos = %d, Neg = %d)', ...
        sum(reward_positive_idx), sum(reward_negative_idx)));
    axis equal; xlim([0 max_val]); ylim([0 max_val]);
    ESN_Beautify_Plot(fig5, [10 4]);
    print(fig5, fullfile(data_path,'population_figs\ephys',...
          'Fig1-C5.pdf'), '-dpdf', '-bestfit');
%% Plot high/low direction
 fig6 = figure;
    for counter_dir = numel(dirs):-1:1
        current_dir = dirs(counter_dir);
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            current_dir_ = [current_dir, 10-current_dir];  % e.g., [2 8] for ±45
        end
    
        % === reward-pos data ===
        cs_rate_vis_p = squeeze(mean(cs_rate_rot_p(:,t_ind_vis,current_dir_,1), 3));
        cs_rate_sac_p = squeeze(mean(cs_rate_rot_p(:,t_ind_sac,current_dir_,2), 3));
    
        if counter_dir == 1
            cs_vm_vis_p = squeeze(mean(cs_vm_rot_p(:,t_ind_vis,:,1), 3));
            cs_vm_sac_p = squeeze(mean(cs_vm_rot_p(:,t_ind_sac,:,2), 3));
        end
    
        % === reward-neg ===
        cs_rate_vis_n = squeeze(mean(cs_rate_rot_n(:,t_ind_vis,current_dir_,1), 3));
        cs_rate_sac_n = squeeze(mean(cs_rate_rot_n(:,t_ind_sac,current_dir_,2), 3));
      
        % VISUAL PANEL
        subplot(1,2,1); 
        title(['rho rwd-pos = ' num2str(reward_positive_rho,...
            '%.2f') ',  rho rwd-neg = ' num2str(reward_negative_rho, '%.2f')], ...
          'FontSize', 8);

        % High
        colororder({'k','k'}); 
        yyaxis left;
        [hl1, hp1] = boundedline(t_ind_vis - 250,...
            mean(cs_rate_vis_p, 1, 'omitnan'),...
            std(cs_rate_vis_p, [], 1, 'omitnan') / sqrt(n_pos),'Color',colors_h(counter_dir,:),'alpha');
       
        hl1.DisplayName = ['pos - ' leg_{6 - counter_dir}];
        hp1.HandleVisibility = 'off';
        % Low
        [hl2, hp2] = boundedline(t_ind_vis - 250,...
            mean(cs_rate_vis_n, 1, 'omitnan'),...
            std(cs_rate_vis_n, [], 1, 'omitnan') / sqrt(n_neg),'alpha');
        set(hl2, 'LineStyle', '--', 'Color', colors_l(counter_dir,:));
        hl2.DisplayName = ['neg - ' leg_{6 - counter_dir}];
        hp2.HandleVisibility = 'off';
    
        hold on; xlabel('time from stimulus onset (ms)');
        ylabel('complex spike rate (Hz)');
        xline(0,'--','HandleVisibility','off');
    
        if counter_dir == 1
            yyaxis right;
            area(t_ind_vis - 250, mean(cs_vm_vis_p, 1, 'omitnan'),...
                'FaceColor', vm_color, 'FaceAlpha', 0.3,...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
            ylim(vm_lim);
        end
    
        % SACCADE PANEL
        subplot(1,2,2);  

        % positive
        colororder({'k','k'});
        yyaxis left;
        [hl3, hp3] = boundedline(t_ind_sac - 250,...
            mean(cs_rate_sac_p, 1, 'omitnan'),...
            std(cs_rate_sac_p, [], 1, 'omitnan') / sqrt(n_pos),...
             'Color', colors_h(counter_dir,:),'alpha');
    
        hl3.DisplayName = ['pos - ' leg_{6-counter_dir}];
        hp3.HandleVisibility = 'off';
    
        % Negative
        [hl4, hp4] = boundedline(t_ind_sac - 250,...
            mean(cs_rate_sac_n, 1, 'omitnan'),...
            std(cs_rate_sac_n, [], 1, 'omitnan') / sqrt(n_neg),'alpha');
        set(hl4, 'LineStyle', '--', 'Color', colors_l(counter_dir,:));
    
        hl4.DisplayName = ['Neg - ' leg_{6-counter_dir}];
        hp4.HandleVisibility = 'off';
    
        hold on; xlabel('time from saccade onset (ms)');
        xline(0,'--','HandleVisibility','off');
    
        if counter_dir == 1
            yyaxis right;
            area(t_ind_sac - 250, mean(cs_vm_sac_p, 1, 'omitnan'),...
                'FaceColor', vm_color, 'FaceAlpha', 0.3,...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
            ylim(vm_lim);
            ylabel('Velocity(deg/sec)')
        end
     
        legend('Box','off', 'Location','northeast');
     end

    ESN_Beautify_Plot(fig6, [10, 4]);
    print(fig6, fullfile(data_path,'population_figs\ephys',...
    'Fig1-C6.pdf'), '-dpdf', '-bestfit');
%% CS rate of reward-positive and negative cells
% Positive cells: extract from cs_rate_rot_h and _l using reward-positive indices
pos_cs_high_vis = squeeze(cs_rate_rot_h(reward_positive_idx, t_ind_vis, dir, 1));
pos_cs_low_vis  = squeeze(cs_rate_rot_l(reward_positive_idx, t_ind_vis, dir, 1));

pos_cs_high_sac = squeeze(cs_rate_rot_h(reward_positive_idx, t_ind_sac, dir, 2));
pos_cs_low_sac  = squeeze(cs_rate_rot_l(reward_positive_idx, t_ind_sac, dir, 2));

% Negative cells: extract from cs_rate_rot_h and _l using reward-negative indices
neg_cs_high_vis = squeeze(cs_rate_rot_h(reward_negative_idx, t_ind_vis, dir, 1));
neg_cs_low_vis  = squeeze(cs_rate_rot_l(reward_negative_idx, t_ind_vis, dir, 1));

neg_cs_high_sac = squeeze(cs_rate_rot_h(reward_negative_idx, t_ind_sac, dir, 2));
neg_cs_low_sac  = squeeze(cs_rate_rot_l(reward_negative_idx, t_ind_sac, dir, 2));

pos_vm_vis = squeeze(cs_vm_rot_h(reward_positive_idx, t_ind_vis, dir, 1));
neg_vm_vis = squeeze(cs_vm_rot_h(reward_negative_idx, t_ind_vis, dir, 1));  

pos_vm_sac = squeeze(cs_vm_rot_h(reward_positive_idx, t_ind_sac, dir, 2));
neg_vm_sac = squeeze(cs_vm_rot_h(reward_negative_idx, t_ind_sac, dir, 2));  
   % Time axis
  
    fig7 = figure;
    % Subplot for reward-positive cells
    subplot(1,2,1)
    hold on
    colororder({'k','k'})
    yyaxis left
    [hl1, hp1] = boundedline(t_ind_vis-250, ...
        mean(pos_cs_high_vis,1,'omitnan'), ...
        std(pos_cs_high_vis,[],1,'omitnan')/sqrt(size(pos_cs_high_vis,1)), ...
        'm','alpha');
    hp1.HandleVisibility = 'off';
    
    [hl2, hp2] = boundedline(t_ind_vis-250, ...
        mean(pos_cs_low_vis,1,'omitnan'), ...
        std(pos_cs_low_vis,[],1,'omitnan')/sqrt(size(pos_cs_low_vis,1)), ...
        'c--','alpha');
    hp2.HandleVisibility = 'off';
    ylabel('CS rate (Hz)')

    yyaxis right;
    area(t_ind_vis - 250, mean(pos_vm_vis, 1, 'omitnan'),...
        'FaceColor', vm_color, 'FaceAlpha', 0.3,...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylim(vm_lim);
     
    
    xline(0,'--k')
    xlabel('Time from target onset (ms)');
    title(sprintf('Reward-Positive Cells (n = %d)', size(pos_cs_high_vis,1)))
    legend([hl1, hl2], {'High', 'Low'}, 'Location','northeast','Box','off');

    
    % Subplot for reward-negative cells
    subplot(1,2,2)
    hold on
    colororder({'k','k'})
    yyaxis left
    [hl3, hp3] = boundedline(t_ind_vis-250, ...
        mean(neg_cs_high_vis,1,'omitnan'), ...
        std(neg_cs_high_vis,[],1,'omitnan')/sqrt(size(neg_cs_high_vis,1)), ...
        'm','alpha');
    hp3.HandleVisibility = 'off';
    
    [hl4, hp4] = boundedline(t_ind_vis-250, ...
        mean(neg_cs_low_vis,1,'omitnan'), ...
        std(neg_cs_low_vis,[],1,'omitnan')/sqrt(size(neg_cs_low_vis,1)), ...
        'c--','alpha');
    hp4.HandleVisibility = 'off';

    yyaxis right;
    area(t_ind_vis - 250, mean(neg_vm_vis, 1, 'omitnan'),...
        'FaceColor', vm_color, 'FaceAlpha', 0.3,...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylim(vm_lim);
    ylabel('velocity(deg/sec)')    
    xline(0,'--k')
    xlabel('Time from target onset (ms)');
    title(sprintf('Reward-Negative Cells (n = %d)', size(neg_cs_high_vis,1)))
    legend([hl3, hl4], {'High', 'Low'}, 'Location','northeast','Box','off');
    ESN_Beautify_Plot(fig7, [10 4])
    print(fig7, fullfile(data_path,'population_figs\ephys','Fig1-C7.pdf'),...
        '-dpdf', '-bestfit');

%% reward-positive and negative cells

    fig8 = figure;
    % Subplot for reward-positive cells
    subplot(1,2,1)
    hold on
    colororder({'k','k'})
    yyaxis left
    [hl1, hp1] = boundedline(t_ind_sac-250, ...
        mean(pos_cs_high_sac,1,'omitnan'), ...
        std(pos_cs_high_sac,[],1,'omitnan')/sqrt(size(pos_cs_high_sac,1)), ...
        'm','alpha');
    hp1.HandleVisibility = 'off';
    ylabel('CS rate (Hz)')
    [hl2, hp2] = boundedline(t_ind_sac-250, ...
        mean(pos_cs_low_sac,1,'omitnan'), ...
        std(pos_cs_low_sac,[],1,'omitnan')/sqrt(size(pos_cs_low_sac,1)), ...
        'c--','alpha');
    hp2.HandleVisibility = 'off';

    yyaxis right;
    area(t_ind_sac - 250, mean(pos_vm_sac, 1, 'omitnan'),...
        'FaceColor', vm_color, 'FaceAlpha', 0.3,...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylim(vm_lim);


    xline(0,'--k')
    xlabel('Time from saccade onset (ms)');
    title(sprintf('Reward-Positive Cells (n = %d)', size(pos_cs_high_sac,1)))
    legend([hl1, hl2], {'High', 'Low'}, 'Location','northeast','Box','off');

    
    % Subplot for reward-negative cells
    subplot(1,2,2)
    hold on
    colororder({'k','k'})
    yyaxis left
    [hl3, hp3] = boundedline(t_ind_sac-250, ...
        mean(neg_cs_high_sac,1,'omitnan'), ...
        std(neg_cs_high_sac,[],1,'omitnan')/sqrt(size(neg_cs_high_sac,1)), ...
        'm','alpha');
    hl3.DisplayName = 'High';
    hp3.HandleVisibility = 'off';
    
    [hl4, hp4] = boundedline(t_ind_sac-250, ...
        mean(neg_cs_low_sac,1,'omitnan'), ...
        std(neg_cs_low_sac,[],1,'omitnan')/sqrt(size(neg_cs_low_sac,1)), ...
        'c--','alpha');
    hl4.DisplayName = 'Low';
    hp4.HandleVisibility = 'off';

    yyaxis right;
    area(t_ind_sac - 250, mean(neg_vm_sac, 1, 'omitnan'),...
        'FaceColor', vm_color, 'FaceAlpha', 0.3,...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylim(vm_lim);
    ylabel('velocity(deg/sec)')      
    xline(0,'--k')
    xlabel('Time from saccade onset (ms)');
    title(sprintf('Reward-Negative Cells (n = %d)', size(neg_cs_high_sac,1)))
    legend([hl3, hl4], {'High', 'Low'}, 'Location','northeast','Box','off');
    ESN_Beautify_Plot(fig8, [10 4])
    print(fig8, fullfile(data_path,'population_figs\ephys','Fig1-C8.pdf'),...
        '-dpdf', '-bestfit');

%% Fig3: Heatmap — All cells sorted by (Low − High), shown side-by-side


% Get CS-on direction and time

cs_high = squeeze(cs_rate_rot_h(:, :, dir, 1));  % [cells × time]
cs_low  = squeeze(cs_rate_rot_l(:, :, dir, 1));  % same

% Sort based on mean difference in target window
t_inds_sort = find(t >= time_window_vis(1) & t <= time_window_vis(2));
mean_high = mean(cs_high(:, t_inds_sort), 2, 'omitnan');
mean_low  = mean(cs_low(:, t_inds_sort), 2, 'omitnan');
diff_sort = mean_low - mean_high;

[~, sort_idx] = sort(diff_sort, 'descend');  % Reward-negative at top

% Z-score each neuron's trace (for better contrast in heatmap)
cs_high_z = (cs_high - mean(cs_high, 2, 'omitnan')) ./ ...
             std(cs_high, 0, 2, 'omitnan');
cs_low_z  = (cs_low  - mean(cs_low, 2, 'omitnan')) ./ ...
             std(cs_low,  0, 2, 'omitnan');

% Apply sorting
cs_high_sorted = cs_high_z(sort_idx, :);
cs_low_sorted  = cs_low_z(sort_idx, :);

fig9 = figure;
% ------------- plot visual aligned ----------------
% Plot low reward
subplot(1,2,1);
imagesc(t_ind_vis - 250, 1:size(cs_low_sorted,1), cs_low_sorted);
xlabel('Time from target onset (ms)');
ylabel('Neurons (sorted by reward modulation)');
title('Low Reward — CS-on Direction');
xline(0, '--w', 'LineWidth', 1.2);
colormap(params.rwb_map);
caxis([0 4]);
colorbar;

% Plot high reward
subplot(1,2,2);
imagesc(t_ind_vis - 250, 1:size(cs_high_sorted,1), cs_high_sorted);
xlabel('Time from target onset (ms)');
ylabel('Neurons (same order)');
title('High Reward — CS-on Direction');
xline(0, '--w', 'LineWidth', 1.2);
colormap(params.rwb_map);
caxis([0 4]);
colorbar;

% Final formatting
ESN_Beautify_Plot(fig9, [10, 4]);
print(fig9, fullfile(data_path,'population_figs\ephys','Fig1-C9.pdf'), '-dpdf', '-bestfit');

% ------------- plot saccade aligend ----------------
% Saccade-aligned heatmaps (CS-on dir)
cs_high = squeeze(cs_rate_rot_h(:, t_ind_sac, dir, 2));  % [cells x time], sac-aligned only
cs_low  = squeeze(cs_rate_rot_l(:, t_ind_sac, dir, 2));  % same

% Time axis for plotting & sorting
time_axis_sac = t_ind_sac - 250;

% Sort based on mean difference in saccade window
t_inds_sort = find(time_axis_sac >= time_window_sac(1) & time_axis_sac <= time_window_sac(2));
mean_high = mean(cs_high(:, t_inds_sort), 2, 'omitnan');
mean_low  = mean(cs_low(:,  t_inds_sort), 2, 'omitnan');
[~, sort_idx] = sort(mean_low - mean_high, 'descend');  % Reward-negative at top

% Z-score each neuron's trace within the saccade-aligned window you plot
cs_high_z = (cs_high - mean(cs_high, 2, 'omitnan')) ./ std(cs_high, 0, 2, 'omitnan');
cs_low_z  = (cs_low  - mean(cs_low,  2, 'omitnan')) ./ std(cs_low,  0, 2, 'omitnan');

% Apply sorting (same order for both panels)
cs_high_sorted = cs_high_z(sort_idx, :);
cs_low_sorted  = cs_low_z( sort_idx, :);

fig10 = figure;

% Low reward
subplot(1,2,1);
imagesc(time_axis_sac, 1:size(cs_low_sorted,1), cs_low_sorted);
xlabel('Time from saccade onset (ms)');
ylabel('Neurons (sorted by reward modulation)');
title('Low Reward — CS-on Direction (saccade-aligned)');
xline(0, '--w', 'LineWidth', 1.2);
colormap(params.rwb_map); caxis([0 4]); colorbar;

% High reward (same neuron order)
subplot(1,2,2);
imagesc(time_axis_sac, 1:size(cs_high_sorted,1), cs_high_sorted);
xlabel('Time from saccade onset (ms)');
ylabel('Neurons (same order)');
title('High Reward — CS-on Direction (saccade-aligned)');
xline(0, '--w', 'LineWidth', 1.2);
colormap(params.rwb_map); caxis([0 4]); colorbar;

ESN_Beautify_Plot(fig10, [10, 4]);
print(fig10, fullfile(data_path,'population_figs\ephys','Fig1-C10.pdf'), '-dpdf', '-bestfit');


%%  Heatmap showing reward-positive vs reward-negative cells (CS-on dir, target-aligned)

% fig10 = figure;
% 
% % Sanity check for dimensions
% if size(cs_rate_rot_p, 4) < 1 || size(cs_rate_rot_n, 4) < 1
%     error('cs_rate_rot_p or cs_rate_rot_n missing proper reward condition dimension.');
% end
% 
% % Get CS-on direction (dir=5), visual epoch (target-aligned)
% 
% cs_rate_pos = squeeze(cs_rate_rot_p(:, t_ind_vis, dir, 1));  % [n_pos × time]
% cs_rate_neg = squeeze(cs_rate_rot_n(:, t_ind_vis, dir, 1));  % [n_neg × time]
% 
% % Z-score across time for each neuron
% cs_rate_pos_z = (cs_rate_pos - mean(cs_rate_pos, 2, 'omitnan')) ./ ...
%                  std(cs_rate_pos, 0, 2, 'omitnan');
% cs_rate_neg_z = (cs_rate_neg - mean(cs_rate_neg, 2, 'omitnan')) ./ ...
%                  std(cs_rate_neg, 0, 2, 'omitnan');
% 
% % Concatenate: reward-positive cells first, reward-negative below
% cs_rate_sorted = [cs_rate_pos_z; cs_rate_neg_z];
% n_pos = size(cs_rate_pos_z, 1);
% n_neg = size(cs_rate_neg_z, 1);
% 
% % Plot heatmap
% imagesc(t_ind_vis - 250, 1:(n_pos + n_neg), cs_rate_sorted);
% xlabel('Time from target onset (ms)');
% ylabel('Neurons (sorted by reward modulation)');
% title('Z-scored CS Rate — CS-on Direction (Visually Aligned)');
% 
% colormap(params.rwb_map);
% caxis([0 4]);
% xline(0, '--w', 'LineWidth', 1.2);
% line([t_ind_vis(1)-250, t_ind_vis(end)-250], [n_pos+0.5, n_pos+0.5], ...
%     'Color', 'w', 'LineStyle', '--', 'LineWidth', 1.2);  % white separator line
% 
% % Optional labels
% text(t_ind_vis(1)-240, n_pos/2, 'Reward-Positive', 'Color', 'w', 'FontSize', 8);
% text(t_ind_vis(1)-240, n_pos + n_neg/2, 'Reward-Negative', 'Color', 'w', 'FontSize', 8);
% 
% colorbar;
% ESN_Beautify_Plot(fig10, [10, 4]);
% print(fig10, fullfile(data_path,'population_figs\ephys','Fig1-C10.pdf'), '-dpdf', '-bestfit');


end 
