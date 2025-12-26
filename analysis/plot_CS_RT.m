function plot_CS_RT(data_path)
% =============================================================
% Plot CS rate (Visual and Saccade epochs) across RT bins
% → Each: 4×2 subplots (CS + velocity), one shared legend
% =============================================================

JDM_params_funcs
S = load(fullfile(data_path,'population_data','cs_on_rate_spe_rpe_prim_sac_RT.mat'));
cs_on_data = S.cs_on_data;
cs_rate    = S.cs_rate;   % [cell × time × dir × reward × RTbin × epoch]
cs_vm      = S.cs_vm;
rt_edges   = S.rt_edges;
num_cs     = numel(cs_on_data);

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
high_conds = 1;
low_conds  = 2;

cs_rate_rot_h = squeeze(mean(cs_rate_rot(:,:,:,high_conds,:,:),4,'omitnan'));
cs_rate_rot_l = squeeze(mean(cs_rate_rot(:,:,:,low_conds,:,:),4,'omitnan'));

cs_vm_rot_h = squeeze(mean(cs_vm_rot(:,:,:,high_conds,:,:),4,'omitnan'));
cs_vm_rot_l = squeeze(mean(cs_vm_rot(:,:,:,low_conds,:,:),4,'omitnan'));

% === Parameters ===
leg_ = {'cs-on','cs±45','cs±90','cs±135','cs±180'};
colors_h = cool(5);
dirs = 1:5;
vm_color = [204,174,98]/256;
vm_lim_vis = [0 400];
vm_lim_sac = [0 600];
t_ind_vis = (-50:250)+250;
t_ind_sac = (-200:200)+250;
num_bins = size(cs_rate_rot,5);

fprintf('Plotting visual & saccade epochs for %d RT bins...\n', num_bins);

%% =============================================================
% === VISUAL EPOCH FIGURE
% =============================================================
fig1 = figure('Color','w');
colororder({'k','k'})
tiledlayout(4,2,'TileSpacing','compact','Padding','compact');

for rb = 1:num_bins
    nexttile
    hold on
    yyaxis left
% numel(dirs):-1:1
    for counter_dir = numel(dirs):-1:1
        current_dir = dirs(counter_dir);
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            continue;
        end

        % Visual data (epoch index = 1)
        cs_rate_vis_h = squeeze(mean(cs_rate_rot_h(:,t_ind_vis,current_dir_,rb,1),3));
        cs_rate_vis_l = squeeze(mean(cs_rate_rot_l(:,t_ind_vis,current_dir_,rb,1),3));

        m_h = mean(cs_rate_vis_h,1,'omitnan');
        e_h = std(cs_rate_vis_h,[],1,'omitnan') / sqrt(num_cs);
        m_l = mean(cs_rate_vis_l,1,'omitnan');
        e_l = std(cs_rate_vis_l,[],1,'omitnan') / sqrt(num_cs);

        [hl1, hp1] = boundedline(t_ind_vis-250, m_h, e_h, ...
            'cmap', colors_h(counter_dir,:), 'alpha');
        set(hl1, 'DisplayName', sprintf('High - %s', leg_{6-counter_dir}));
        set(hp1, 'HandleVisibility','off');

        [hl2, hp2] = boundedline(t_ind_vis-250, m_l, e_l, ...
            'cmap', colors_h(counter_dir,:), 'alpha');
        set(hl2, 'LineStyle','--', 'DisplayName', sprintf('Low - %s', leg_{6-counter_dir}));
        set(hp2, 'HandleVisibility','off');
     
    end

    ylabel('CS rate (Hz)');
    xline(0,'--k','HandleVisibility','off');
    ylim([0, 6]);  

    yyaxis right
    hold on
    cs_vm_vis_h = squeeze(mean(cs_vm_rot_h(:,t_ind_vis,1:8,rb,1),[1 3],'omitnan'));
    cs_vm_vis_l = squeeze(mean(cs_vm_rot_l(:,t_ind_vis,1:8,rb,1),[1 3],'omitnan'));
    area(t_ind_vis-250, cs_vm_vis_h, 'FaceColor', vm_color, 'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', 'DisplayName','Vel High');
    area(t_ind_vis-250, cs_vm_vis_l, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', 'DisplayName','Vel Low');
    ylim(vm_lim_vis);
    ylabel('Velocity (deg/s)');

    title(sprintf('RT %.0f–%.0f ms', rt_edges(rb), rt_edges(rb+1)));
    xlabel('Time from stimulus onset (ms)');
    hold off
end

lgd1 = legend('Box','off','NumColumns',3,'Location','southoutside');
lgd1.Layout.Tile = 'south';

ESN_Beautify_Plot(fig1,[10,12]);
save_path1 = fullfile(data_path,'population_figs','ephys','CS_RT_TAG1_visual.pdf');
print(fig1,save_path1,'-dpdf','-bestfit')
fprintf('Saved figure → %s\n', save_path1);

%% =============================================================
% === SACCADE EPOCH FIGURE
% =============================================================
fig2 = figure('Color','w');
colororder({'k','k'})
tiledlayout(4,2,'TileSpacing','compact','Padding','compact');

for rb = 1:num_bins
    nexttile
    hold on
    yyaxis left

    for counter_dir = numel(dirs):-1:1
        current_dir = dirs(counter_dir);
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            continue;
        end

        % Saccade data (epoch index = 2)
        cs_rate_sac_h = squeeze(mean(cs_rate_rot_h(:,t_ind_sac,current_dir_,rb,2),3));
        cs_rate_sac_l = squeeze(mean(cs_rate_rot_l(:,t_ind_sac,current_dir_,rb,2),3));

        m_h = mean(cs_rate_sac_h,1,'omitnan');
        e_h = std(cs_rate_sac_h,[],1,'omitnan') / sqrt(num_cs);
        m_l = mean(cs_rate_sac_l,1,'omitnan');
        e_l = std(cs_rate_sac_l,[],1,'omitnan') / sqrt(num_cs);

        [hl1, hp1] = boundedline(t_ind_sac-250, m_h, e_h, ...
            'cmap', colors_h(counter_dir,:), 'alpha');
        set(hl1, 'DisplayName', sprintf('High - %s', leg_{6-counter_dir}));
        set(hp1, 'HandleVisibility','off');

        [hl2, hp2] = boundedline(t_ind_sac-250, m_l, e_l, ...
            'cmap', colors_h(counter_dir,:), 'alpha');
        set(hl2, 'LineStyle','--', 'DisplayName', sprintf('Low - %s', leg_{6-counter_dir}));
        set(hp2, 'HandleVisibility','off');
    end

    ylabel('CS rate (Hz)');
    xline(0,'--k','HandleVisibility','off');
    ylim([0, 6]); 


    yyaxis right
    hold on
    cs_vm_sac_h = squeeze(mean(cs_vm_rot_h(:,t_ind_sac,1:8,rb,2),[1 3],'omitnan'));
    cs_vm_sac_l = squeeze(mean(cs_vm_rot_l(:,t_ind_sac,1:8,rb,2),[1 3],'omitnan'));
    area(t_ind_sac-250, cs_vm_sac_h, 'FaceColor', vm_color, 'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', 'DisplayName','Vel High');
    area(t_ind_sac-250, cs_vm_sac_l, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', 'DisplayName','Vel Low');
    ylim(vm_lim_sac);
    ylabel('Velocity (deg/s)');

    title(sprintf('RT %.0f–%.0f ms', rt_edges(rb), rt_edges(rb+1)));
    xlabel('Time from saccade onset (ms)');
    hold off
end

lgd2 = legend('Box','off','NumColumns',3,'Location','southoutside');
lgd2.Layout.Tile = 'south';

ESN_Beautify_Plot(fig2,[10,12]);
save_path2 = fullfile(data_path,'population_figs','ephys','CS_RT_TAG1-saccade.pdf');
print(fig2,save_path2,'-dpdf','-bestfit')
fprintf('Saved figure → %s\n', save_path2);

end
