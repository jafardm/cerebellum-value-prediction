function plot_CS_Choice_RT(data_path)
% =============================================================
% Plot CS rate (Visual and Saccade epochs) for choice trials across RT bins
% → fig1: 4×2 subplots (CS + velocity) – VISUAL epoch, split by RT bin
% → fig2: 4×2 subplots (CS + velocity) – SACCADE epoch, split by RT bin
% → fig3: Avg across RT bins – CS-on & CS-on+180, Choice High vs Low
% =============================================================

JDM_params_funcs
S = load(fullfile(data_path,'population_data','cs_on_rate_choice_sac_RT.mat'));
cs_on_data = S.cs_on_data;
cs_rate    = S.cs_rate;   % [cell × time × dir × reward × RTbin × epoch]
cs_vm      = S.cs_vm;
rt_edges   = S.rt_edges;
n_RT_bins   = numel(rt_edges)-1;
num_cs     = numel(cs_on_data);
animal_list         = S.animal_list;
choice_perf_RT_sess = S.choice_perf_RT;
Choice_RTs = S.choice_RT_bins;

fprintf('Loaded %d CS cells\n', num_cs);

% === Find CS-on direction for each neuron ===
cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_data,'UniformOutput',false));
cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_data,'UniformOutput',false));
cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, cs_on_data,'UniformOutput',false));
cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, cs_on_data,'UniformOutput',false));

ind_m = cs_on_rho_sac > cs_on_rho_vis;
cs_on_ang = cs_on_ang_vis/180*pi;
cs_on_ang(ind_m) = cs_on_ang_sac(ind_m)/180*pi;

ang_bins = params.sac.ang_values/180*pi;
[~, cs_on_bin] = min(abs(angdiff(repmat(ang_bins,num_cs,1), cs_on_ang*ones(size(ang_bins)))),[],2);
order_ang = nan(num_cs,8);
for c = 1:num_cs
    order_ang(c,:) = circshift(1:8,5-cs_on_bin(c));
end

% === Rotate rate & velocity into CS-on frame ===
cs_rate_rot = nan(size(cs_rate));
cs_vm_rot   = nan(size(cs_vm));
for c = 1:num_cs
    cs_rate_rot(c,:,:,:,:,:) = cs_rate(c,:,order_ang(c,1:8),:,:,:);
    cs_vm_rot(c,:,:,:,:,:)   = cs_vm(c,:,order_ang(c,1:8),:,:,:);
end

% === Separate High vs Low ===
choice_high = 1;
choice_low  = 2;

% cs_rate_rot: [cell x time x dir x reward x RTbin x epoch]
cs_rate_rot_h = squeeze(mean(cs_rate_rot(:,:,:,choice_high,:,:),4,'omitnan')); % [cell x time x dir x RTbin x epoch]
cs_rate_rot_l = squeeze(mean(cs_rate_rot(:,:,:,choice_low,:,:),4,'omitnan'));  % [cell x time x dir x RTbin x epoch]

cs_vm_rot_h = squeeze(mean(cs_vm_rot(:,:,:,choice_high,:,:),4,'omitnan'));
cs_vm_rot_l = squeeze(mean(cs_vm_rot(:,:,:,choice_low,:,:),4,'omitnan'));

% === Parameters ===
leg_ = {'cs-on','cs+45','cs+90','cs+135','cs+180'};
colors_h = cool(5);
vm_color = [204,174,98]/256;
vm_lim_vis = [0 500];
vm_lim_sac = [0 600];
t_ind_vis = (-50:250)+250;
t_ind_sac = (-200:200)+250;
num_bins = size(cs_rate_rot,5);
dirs = 1:5;
dirs_to_plot   = [5 1];   % cs-on, cs±180 
leg_idx        = [1 5];   % indices into leg_

%% =============================================================
% === VISUAL EPOCH FIGURE (per RT bin)
% =============================================================

for a = 1:numel(animal_list)
    
    if a==1
        cell_idx1 = 1;
        cell_idx2 = 480;
    else
        cell_idx1 = 481;
        cell_idx2 = 906;
    end

    fig1 = figure('Color','w');
    colororder({'k','k'})
    t1 = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
    title(t1,[animal_list{a} '-CS Rate Choice Trials – Visual Epoch']);
    
    for rb = 1:num_bins
        nexttile
        hold on
        yyaxis left
    
        for k = 1:numel(dirs_to_plot)
            current_dir_ = dirs_to_plot(k);
            label_str    = leg_{leg_idx(k)};
    
            % Visual data (epoch index = 1)
            cs_rate_vis_h = squeeze(mean(cs_rate_rot_h(cell_idx1:cell_idx2,t_ind_vis,current_dir_,rb,1),3));
            cs_rate_vis_l = squeeze(mean(cs_rate_rot_l(cell_idx1:cell_idx2,t_ind_vis,current_dir_,rb,1),3));
    
            m_h = mean(cs_rate_vis_h,1,'omitnan');
            e_h = std(cs_rate_vis_h,[],1,'omitnan') / sqrt(num_cs);
            m_l = mean(cs_rate_vis_l,1,'omitnan');
            e_l = std(cs_rate_vis_l,[],1,'omitnan') / sqrt(num_cs);
    
            [hl1, hp1] = boundedline(t_ind_vis-250, m_h, e_h, ...
                'cmap', colors_h(current_dir_,:), 'alpha');
            set(hl1, 'DisplayName', sprintf('Choice-High - %s', label_str));
            set(hp1, 'HandleVisibility','off');
    
            [hl2, hp2] = boundedline(t_ind_vis-250, m_l, e_l, ...
                'cmap', colors_h(current_dir_,:), 'alpha');
            set(hl2, 'LineStyle','--', 'DisplayName', sprintf('Choice-Low - %s', label_str));
            set(hp2, 'HandleVisibility','off');
        end
    
        ylabel('CS rate (Hz)');
        xline(0,'--k','HandleVisibility','off');
        ylim([0, 6.5]);  
    
        yyaxis right
        hold on
        cs_vm_vis_h = squeeze(mean(cs_vm_rot_h(cell_idx1:cell_idx2,t_ind_vis,1:8,rb,1),[1 3],'omitnan'));
        cs_vm_vis_l = squeeze(mean(cs_vm_rot_l(cell_idx1:cell_idx2,t_ind_vis,1:8,rb,1),[1 3],'omitnan'));
        area(t_ind_vis-250, cs_vm_vis_h, 'FaceColor', vm_color, 'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', 'DisplayName','Vel Choice-High');
        area(t_ind_vis-250, cs_vm_vis_l, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', 'DisplayName','Vel Choice-Low');
        ylim(vm_lim_vis);
        ylabel('Velocity (deg/s)');
    
        title(sprintf('RT %.0f–%.0f ms', rt_edges(rb), rt_edges(rb+1)));
        xlabel('Time from stimulus onset (ms)');
        hold off
    end
    
    lgd1 = legend('Box','off','NumColumns',3,'Location','southoutside');
    lgd1.Layout.Tile = 'south';
    
    ESN_Beautify_Plot(fig1,[10,12]);
    save_path1 = fullfile(data_path,'population_figs','ephys',...
        [animal_list{a} '-CS_RT_Choice_visual.pdf']);
    %print(fig1,save_path1,'-dpdf','-bestfit')

end

%% =============================================================
% === SACCADE EPOCH FIGURE (per RT bin)
% =============================================================

for a = 1:numel(animal_list)
    
    if a==1
        cell_idx1 = 1;
        cell_idx2 = 480;
    else
        cell_idx1 = 481;
        cell_idx2 = 906;
    end
    fig2 = figure('Color','w');
    colororder({'k','k'})
    t2 = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
     title(t2,[animal_list{a} '-CS Rate Choice Trials – Saccade Epoch']);
    
    for rb = 1:num_bins
        nexttile
        hold on
        yyaxis left
    
        for k = 1:numel(dirs_to_plot)
            current_dir_ = dirs_to_plot(k);
            label_str    = leg_{leg_idx(k)};
    
            % Saccade data (epoch index = 2)
            cs_rate_sac_h = squeeze(mean(cs_rate_rot_h(cell_idx1:cell_idx2,t_ind_sac,current_dir_,rb,2),3));
            cs_rate_sac_l = squeeze(mean(cs_rate_rot_l(cell_idx1:cell_idx2,t_ind_sac,current_dir_,rb,2),3));
    
            m_h = mean(cs_rate_sac_h,1,'omitnan');
            e_h = std(cs_rate_sac_h,[],1,'omitnan') / sqrt(num_cs);
            m_l = mean(cs_rate_sac_l,1,'omitnan');
            e_l = std(cs_rate_sac_l,[],1,'omitnan') / sqrt(num_cs);
    
            [hl1, hp1] = boundedline(t_ind_sac-250, m_h, e_h, ...
                'cmap', colors_h(current_dir_,:), 'alpha');
            set(hl1, 'DisplayName', sprintf('Choice-High - %s', label_str));
            set(hp1, 'HandleVisibility','off');
    
            [hl2, hp2] = boundedline(t_ind_sac-250, m_l, e_l, ...
                'cmap', colors_h(current_dir_,:), 'alpha');
            set(hl2, 'LineStyle','--', 'DisplayName', sprintf('Choice-Low - %s', label_str));
            set(hp2, 'HandleVisibility','off');
        end
    
        ylabel('CS rate (Hz)');
        xline(0,'--k','HandleVisibility','off');
        ylim([0, 6.5]); 
    
        yyaxis right
        hold on
        cs_vm_sac_h = squeeze(mean(cs_vm_rot_h(cell_idx1:cell_idx2,t_ind_sac,1:8,rb,2),[1 3],'omitnan'));
        cs_vm_sac_l = squeeze(mean(cs_vm_rot_l(cell_idx1:cell_idx2,t_ind_sac,1:8,rb,2),[1 3],'omitnan'));
        area(t_ind_sac-250, cs_vm_sac_h, 'FaceColor', vm_color, 'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', 'DisplayName','Vel Choice-High');
        area(t_ind_sac-250, cs_vm_sac_l, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', 'DisplayName','Vel Choice-Low');
        ylim(vm_lim_sac);
        ylabel('Velocity (deg/s)');
    
        title(sprintf('RT %.0f–%.0f ms', rt_edges(rb), rt_edges(rb+1)));
        xlabel('Time from saccade onset (ms)');
        hold off
    end
    
    lgd2 = legend('Box','off','NumColumns',3,'Location','southoutside');
    lgd2.Layout.Tile = 'south';
    
    ESN_Beautify_Plot(fig2,[10,12]);
    save_path2 = fullfile(data_path,'population_figs','ephys',...
        [animal_list{a} '-CS_RT_Choice_saccade.pdf']);
    %print(fig2,save_path2,'-dpdf','-bestfit')
end

%% =============================================================
% ===  AVERAGE ACROSS RT BINS (cs-on & cs-on+180)
% =============================================================
% Collapse RT bins: [cell x time x dir x epoch]
cs_rate_rot_h_avgRT = squeeze(mean(cs_rate_rot_h,4,'omitnan'));
cs_rate_rot_l_avgRT = squeeze(mean(cs_rate_rot_l,4,'omitnan'));
cs_vm_rot_h_avgRT   = squeeze(mean(cs_vm_rot_h,4,'omitnan'));
cs_vm_rot_l_avgRT   = squeeze(mean(cs_vm_rot_l,4,'omitnan'));


%% =============================================================
%  PLOT: Chosen-direction tuning for HIGH vs LOW choices
% =============================================================

dirs_all = 1:8;                       
dir_labels = { ...
    'CS+180', ...  % index 1
    'CS-135', ...  % index 2
    'CS-90',  ...  % index 3
    'CS-45',  ...  % index 4
    'CS-on',  ...  % index 5
    'CS+45',  ...  % index 6
    'CS+90',  ...  % index 7
    'CS+135'};     % index 8


colors = [
    0.00 0.00 0.70;  % 1: CS+180  (dark blue)
    0.00 0.50 0.80;  % 2: CS-135  (teal)
    0.00 0.70 0.30;  % 3: CS-90   (green)
    0.40 0.80 0.00;  % 4: CS-45   (yellow-green)
    0.80 0.20 0.00;  % 5: CS-on   (red)
    0.85 0.55 0.00;  % 6: CS+45   (orange)
    0.60 0.20 0.80;  % 7: CS+90   (purple)
    0.30 0.00 0.40]; % 8: CS+135  (dark purple)



fig3 = figure;
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

%  SUBPLOT 1 — HIGH chosen

nexttile; hold on;
title('HIGH chosen: CS rate by target direction');


for d = dirs_all
    cs_high = squeeze(cs_rate_rot_h_avgRT(:, t_ind_vis, d, 1));
    m_h = mean(cs_high,1,'omitnan');
    e_h = std(cs_high,[],1,'omitnan')/sqrt(num_cs);

    [hl,hp] = boundedline(t_ind_vis-250, m_h, e_h, ...
        'cmap', colors(d,:), 'alpha');
    set(hl,'LineWidth',1.5,'DisplayName',dir_labels{d});
    set(hp,'HandleVisibility','off')
end

xlabel('Time from stimulus onset (ms)');
ylabel('CS rate (Hz)');
xline(0,'--k','HandleVisibility','off');
legend('Box','off','Location','northeast');
ylim([0 4]);

yyaxis right 
cs_vm_vis_h_avg = squeeze(mean(cs_vm_rot_h_avgRT(:,t_ind_vis,1:8,1),[1 3],'omitnan'));
area(t_ind_vis-250, cs_vm_vis_h_avg, 'FaceColor', vm_color, 'FaceAlpha', 0.3, ... 
     'EdgeColor','none','HandleVisibility','off'); ylim(vm_lim_vis);
ylabel('Velocity (deg/s)');
xlabel('Time from stimulus onset (ms)'); 
title('high chosen'); legend('Box','off');

%  SUBPLOT 2 — LOW chosen

nexttile; hold on;
title('LOW chosen: CS rate by target direction');

for d = dirs_all
    cs_low = squeeze(cs_rate_rot_l_avgRT(:, t_ind_vis, d, 1));
    m_l = mean(cs_low,1,'omitnan');
    e_l = std(cs_low,[],1,'omitnan')/sqrt(num_cs);

    [hl,hp] = boundedline(t_ind_vis-250, m_l, e_l, ...
        'cmap', colors(d,:), 'alpha');
    set(hl,'LineWidth',1.5,'DisplayName',dir_labels{d},'LineStyle','--');
    set(hp,'HandleVisibility','off');
end

xlabel('Time from stimulus onset (ms)');
ylabel('CS rate (Hz)');
xline(0,'--k','HandleVisibility','off');
legend('Box','off','Location','northeast');
ylim([0 4]);
yyaxis right 
cs_vm_vis_l_avg = squeeze(mean(cs_vm_rot_l_avgRT(:,t_ind_vis,1:8,1),[1 3],'omitnan'));
area(t_ind_vis-250, cs_vm_vis_l_avg, 'FaceColor', vm_color, 'FaceAlpha', 0.3, ... 
     'EdgeColor','none','HandleVisibility','off'); ylim(vm_lim_vis);
ylabel('Velocity (deg/s)');
xlabel('Time from stimulus onset (ms)'); 
title('low chosen'); legend('Box','off');


%% ----- Difference Plot: Visual epoch (tile 3) -----
nexttile
hold on

% === Extract CS-on (0°) and CS-on+180° (5th direction) ===
dir_on      = 5;
dir_on180   = 1;

% --- High and Low CS rates averaged over RT bins ---
H_on_vis    = squeeze(mean(cs_rate_rot_h_avgRT(:,t_ind_vis,dir_on,1),3));     % [cell × time]
L_on_vis    = squeeze(mean(cs_rate_rot_l_avgRT(:,t_ind_vis,dir_on,1),3));

H_180_vis   = squeeze(mean(cs_rate_rot_h_avgRT(:,t_ind_vis,dir_on180,1),3));
L_180_vis   = squeeze(mean(cs_rate_rot_l_avgRT(:,t_ind_vis,dir_on180,1),3));

% === Differences ===
% Low(cs-on) – Low(cs-on+180)
diff1 = (L_on_vis - L_180_vis);      % [cell × time]

% High(cs-on) – High(cs-on+180)
diff2 = (H_on_vis - H_180_vis);

% --- Mean and SEM ---
m1 = mean(diff1,1,'omitnan');
e1 = std(diff1,[],1,'omitnan')/sqrt(num_cs);

m2 = mean(diff2,1,'omitnan');
e2 = std(diff2,[],1,'omitnan')/sqrt(num_cs);

% --- Plot ---
[hl1,hp1] = boundedline(t_ind_vis-250, m1, e1, ...
    'cmap', [0.2 0.6 1], 'alpha'); 
set(hl1,'LineWidth',1.5,'DisplayName','Low \theta - Low\theta+180');
set(hp1,'HandleVisibility','off');

[hl2,hp2] = boundedline(t_ind_vis-250, m2, e2, ...
    'cmap', [1 0.3 0.3], 'alpha');
set(hl2,'LineWidth',1.5,'DisplayName','High\theta - High\theta+180');
set(hp2,'HandleVisibility','off');

xline(0,'--k','HandleVisibility','off');
ylabel('Δ CS rate (Hz)');
xlabel('Time from stimulus onset (ms)');
title('Visual epoch: directional difference');
ax = gca;                        
legend(ax, 'Box','off','Location','northeast'); 

ESN_Beautify_Plot(fig3,[10,4]);
save_path3 = fullfile(data_path,'population_figs','ephys','CS_Choice_trials.pdf');
% print(fig3,save_path3,'-dpdf','-bestfit');    
 %% ------ RT bin centers for plotting ----------------

    rt_centers  = (rt_edges(1:end-1) + rt_edges(2:end)) / 2; 
    num_animals = numel(choice_perf_RT_sess);
    
    fig4 = figure; clf;
    set(gcf,'Position',[200 200 1400 400]);
   
    for a = 1:num_animals
        subplot(1,num_animals,a); hold on;
        animal_name = animal_list{a};
    
        % ---- sessions for this animal ----
        perf_sessions = choice_perf_RT_sess{a};   % {num_sessions × 1}
        num_sessions  = numel(perf_sessions);
    
        % ---- Build padded session matrix ----
        perf_mat = nan(num_sessions, n_RT_bins);
    
        for s = 1:num_sessions
            vec = perf_sessions{s};
            perf_mat(s, 1:numel(vec)) = vec;
        end
    
        % ---- Plot each session trace ----
        for s = 1:num_sessions
            plot(rt_centers, perf_mat(s,:), '-', 'Color',[0.6 0.6 0.6 0.6], ... 
                 'LineWidth',1.2);  
        end
    
        % ---- Compute mean curve ----
        mean_perf = mean(perf_mat,1,'omitnan');
    
        % ---- Overlay mean curve ----
        plot(rt_centers, mean_perf, 'o-', 'LineWidth',3, 'Color','k');
    
        xlabel('Reaction Time (ms)');
        ylabel('P(High | RT bin)');
        title(sprintf('Subject %s', animal_name));
    
        ylim([0 1.2]);
        xlim([rt_centers(1) rt_centers(end)]);
    end
    
    ESN_Beautify_Plot(fig4,[10,4]);

%     save_path4 = fullfile(data_path,'population_figs','ephys','choice_RT_probability.pdf');
%     %print(fig4,save_path4,'-dpdf','-bestfit');


%% =============== Differnece high-low ===================
dir_cs_on = 5;
epoch_vis = 1;

t0 = 300;  % index of visual stimulus onset

rt_centers = (rt_edges(1:end-1) + rt_edges(2:end)) / 2;

for a = 1:numel(animal_list)

    % =============================
    % Determine cell indices per animal
    % =============================
    if a == 1
        cell_idx1 = 1;
        cell_idx2 = 480;
    else
        cell_idx1 = 481;
        cell_idx2 = size(cs_rate_rot_h,1);
    end

    n_cells_animal = cell_idx2 - cell_idx1 + 1;

    % =============================
    % Allocate per-animal diff array
    % =============================
    cs_diff = nan(n_cells_animal, num_bins);

    % =============================
    % Loop RT bins
    % =============================
    for b = 1:num_bins

        rt_max = rt_edges(b+1); % ms
        t_end  = t0 + round(rt_max);

        t_end = min(t_end, size(cs_rate_rot_h,2));

        % --- high ---
        h_tmp = squeeze(mean( ...
            cs_rate_rot_h(cell_idx1:cell_idx2, t0:t_end, dir_cs_on, b, epoch_vis), ...
            2, 'omitnan'));

        % --- low ---
        l_tmp = squeeze(mean( ...
            cs_rate_rot_l(cell_idx1:cell_idx2, t0:t_end, dir_cs_on, b, epoch_vis), ...
            2, 'omitnan'));

        % --- store High-Low per cell ---
        cs_diff(:, b) = h_tmp - l_tmp;
    end

    % =============================
    % Compute mean and SEM (per bin)
    % =============================
    mean_diff = mean(cs_diff, 1, 'omitnan');

    n_per_bin = sum(~isnan(cs_diff),1);
    sem_diff  = std(cs_diff,[],1,'omitnan') ./ sqrt(n_per_bin);

    % =============================
    % Plot per animal
    % =============================
    fig5 = figure; hold on;
    errorbar(rt_centers, mean_diff, sem_diff, 'k-o','LineWidth',2);

    xlabel('Reaction Time (ms)');
    ylabel('\Delta CS rate (High - Low)');
    title(sprintf('Delta CS (50ms to onset) – %s', animal_list{a}));

    ESN_Beautify_Plot(fig5,[8,4]);

    save_path5 = fullfile(data_path,'population_figs','ephys', ...
        sprintf('choice_RT_Delta_CS_rate_%s.pdf', animal_list{a}));

    % print(fig5, save_path5, '-dpdf','-bestfit');
end

%% distribution of the reaction times

RT_dist_all = cell(numel(animal_list), n_RT_bins);

for a = 1:numel(animal_list)
    sess_list = Choice_RTs{a};
    num_sess  = numel(sess_list);

    for b = 1:n_RT_bins
        RT_concat = [];

        for s = 1:num_sess
            RT_concat = [RT_concat; sess_list{s}{b}];  % concatenate RTs
        end

        RT_dist_all{a,b} = RT_concat;
    end
end


for a = 1:numel(animal_list)

    fig = figure('Color','w'); hold on;

   
    data_for_violin = RT_dist_all(a,:);  % 1×8 cell array

    % Find max length
    max_len = max(cellfun(@numel, data_for_violin));
    
    % Pre-allocate matrix
    mat = nan(max_len, n_RT_bins);
    
    % Fill each column
    for b = 1:n_RT_bins
        v = data_for_violin{b};
        mat(1:numel(v), b) = v;
    end

    violinplot(mat, ...
               arrayfun(@(x) sprintf('Bin %d',x), 1:n_RT_bins, 'UniformOutput',false));
    
    ylabel('Reaction Time (ms)');
    xlabel('RT bins');
    title(sprintf('RT Distribution per Bin -  %s', animal_list{a}));
    
    set(gca,'TickDir','out'); box off;
    ESN_Beautify_Plot(fig,[8,4]);
end
%% performace of 65F in choice trials H1/H2

high1_performance_32F = cell2mat(S.perf_H1_animal{1});
high2_performance_32F = cell2mat(S.perf_H2_animal{1});

% number of session
n_high1_32F = sum(~isnan(high1_performance_32F));
n_high2_32F = sum(~isnan(high2_performance_32F));

% Concatenate into a matrix for violinplot
data_matrix_32F = [high1_performance_32F, high2_performance_32F];


high1_performance_65F = cell2mat(S.perf_H1_animal{2});
high2_performance_65F = cell2mat(S.perf_H2_animal{2});
n_high1_65F = sum(~isnan(high1_performance_65F));
n_high2_65F = sum(~isnan(high2_performance_65F));

% Concatenate into a matrix for violinplot
data_matrix_65F = [high1_performance_65F, high2_performance_65F];


% Category labels
labels = {'High 1 ', 'High 2'};
fig = figure('Color','w'); 
% ============================
% 32F
% ============================
subplot(1,2,1)
violinplot(data_matrix_32F, labels);
ylabel('P(Choose High)');
xlabel('High-Stimuli');

title(sprintf('32F: H1 vs H2 (n = %d, %d)', ...
    n_high1_32F, n_high2_32F));

set(gca, 'FontSize', 12, 'TickDir', 'out');
ylim([0.5 1]);
box off;

% ============================
% 65F
% ============================
subplot(1,2,2)
violinplot(data_matrix_65F, labels);
ylabel('P(Choose High)');
xlabel('High-Stimuli');

title(sprintf('65F: H1 vs H2 (n = %d, %d)', ...
    n_high1_65F, n_high2_65F));

set(gca, 'FontSize', 12, 'TickDir', 'out');
ylim([0.5 1]);
box off;

ESN_Beautify_Plot(fig, [10, 4]);
save_path = fullfile(data_path,'population_figs','ephys', ...
'choice_performance_H1_H2.pdf');

 %print(fig, save_path, '-dpdf','-bestfit');


end
