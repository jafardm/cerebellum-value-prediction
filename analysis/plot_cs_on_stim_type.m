function plot_cs_on_stim_type(data_path)

    JDM_params_funcs
    S = load(fullfile(data_path,'population_data','cs_on_rate_prim_sac_stim_type.mat'));

    animal_list     = S.animal_list;
    cs_rate_all     = S.cs_rate;     % [cell x time x dir x 4 x 2]
    cs_vm_all       = S.cs_vm;       % [cell x time x dir x 4 x 2]
    cs_on_data_all  = S.cs_on_data;  % {cell}
    cs_rt_all       = S.cs_rt;       % [cell x dir x 4]

    % -------- 132F --------
    idx_132 = strcmp(animal_list,'132F');
    if any(idx_132)
        cs_rate_132    = cs_rate_all{idx_132};
        cs_vm_132      = cs_vm_all{idx_132};
        cs_rt_132      = cs_rt_all{idx_132};
        cs_on_data_132 = cs_on_data_all{idx_132};
        plot_one_animal('132F', cs_rate_132, cs_vm_132, cs_rt_132, ...
                        cs_on_data_132, data_path, params);
    end

    % -------- 65F --------
    idx_65 = strcmp(animal_list,'65F');
    if any(idx_65)
        cs_rate_65    = cs_rate_all{idx_65};
        cs_vm_65      = cs_vm_all{idx_65};
        cs_rt_65      = cs_rt_all{idx_65};
        cs_on_data_65 = cs_on_data_all{idx_65};
        plot_one_animal('65F', cs_rate_65, cs_vm_65, cs_rt_65, ...
                        cs_on_data_65, data_path, params);
    end
end


function plot_one_animal(sub_name, cs_rate, cs_vm, cs_rt, ...
                         cs_on_data, data_path, params)

    % cs_rate: [num_cs x 500 x 8 x 4(stimType) x 2(epoch: vis,sac)]
    % cs_vm  : same shape as cs_rate
    % cs_rt  : [num_cs x 8 x 4]  (RT in ms)
    % stimType: 1 = High1, 2 = High2, 3 = Low1, 4 = Low2

    % ===============================
    % 1) Compute CS-on alignment
    % ===============================
    cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg,cs_on_data,...
        'UniformOutput',false));
    cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg,cs_on_data,...
        'UniformOutput',false));

    cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg,cs_on_data,...
        'UniformOutput',false));
    cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg,cs_on_data,...
        'UniformOutput',false));

    num_cs = numel(cs_on_data);

    % choose larger vector (vis vs sac) per cell
    ind_m = cs_on_rho_sac > cs_on_rho_vis;

    cs_on_ang = cs_on_ang_vis/180*pi;
    cs_on_ang(ind_m) = cs_on_ang_sac(ind_m)/180*pi;

    ang_bins = params.sac.ang_values/180*pi;  % 8 bins

    [~,cs_on_bin] = min(abs(angdiff(repmat(ang_bins,num_cs,1),...
        cs_on_ang*ones(size(ang_bins)))),[],2);

    order_ang = nan(num_cs,8);
    for counter_cs = 1:num_cs
        order_ang(counter_cs,:) = circshift(1:8,5-cs_on_bin(counter_cs));
    end

    % ===============================
    % 2) Rotate cs_rate, cs_vm, cs_rt by CS-on
    % ===============================
    cs_rate_rot = nan(num_cs,500,8,4,2);
    cs_vm_rot   = nan(num_cs,500,8,4,2);
    cs_rt_rot   = nan(num_cs,8,4);

    for counter_cs = 1:num_cs
        ang_order = order_ang(counter_cs,1:8);

        % time x dir x cond x epoch
        cs_rate_rot(counter_cs,:,:,:,:) = ...
            cs_rate(counter_cs,:,ang_order,:,1:2); 

        cs_vm_rot(counter_cs,:,:,:,:) = ...
            cs_vm(counter_cs,:,ang_order,:,1:2);

        % dir x cond
        cs_rt_rot(counter_cs,:,:) = ...
            cs_rt(counter_cs,ang_order,:); 
    end

    % ===============================
    % 3) Global NaN cleanup (cells with all NaN rates)
    % ===============================
    all_nan_cell = all(all(all(all(isnan(cs_rate_rot),5),4),3),2);
    cs_rate_rot(all_nan_cell,:,:,:,:) = [];
    cs_vm_rot(all_nan_cell,:,:,:,:)   = [];
    cs_rt_rot(all_nan_cell,:,:)       = [];
    num_cs = size(cs_rate_rot,1);

    if num_cs == 0
        warning('No valid CS cells for %s', sub_name);
        return;
    end

    % ===============================
    % 4) Plot settings
    % ===============================
    % dir 5 = CS-on (after rotation)
    % dir 1 = CS-on+180 (after rotation)
    dir_on      = 5;
    dir_180     = 1;

    cond_names  = {'High1','High2','Low1','Low2'};
    n_conds     = 4;
    colors_cond = lines(n_conds);  % one color per condition

    t_ind_vis = (-50:200) + 250;   % 251 points
    t_ind_sac = (-100:150) + 250;  % 251 points
    t_vis     = t_ind_vis - 250;
    t_sac     = t_ind_sac - 250;

    % ===============================
    % 5) Figure: 3x2 layout per animal
    % Row 1: CS rate (vis, sac)
    % Row 2: eye velocity (vis, sac)
    % Row 3: RT bar plot (all conds; spans both columns)
    % ===============================
    fig = figure('Color','w');
    tl = tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('%s: CS rate, eye velocity, and RT (CS-on & CS-on+180)', sub_name));

    % ---------- VISUAL CS RATE ----------
    nexttile(tl,1);
    hold on;

    legend_entries_vis = {};
    for cond = 1:n_conds
        % CS-on, visual
        data_vis_on = squeeze(cs_rate_rot(:, t_ind_vis, dir_on, cond, 1));
        good_cells_on = ~all(isnan(data_vis_on),2);
        data_vis_on = data_vis_on(good_cells_on,:);

        % CS-on+180, visual
        data_vis_180 = squeeze(cs_rate_rot(:, t_ind_vis, dir_180, cond, 1));
        good_cells_180 = ~all(isnan(data_vis_180),2);
        data_vis_180 = data_vis_180(good_cells_180,:);

        % skip if nothing
        if isempty(data_vis_on) && isempty(data_vis_180)
            continue;
        end

        % CS-on
        if ~isempty(data_vis_on)
            m_on  = mean(data_vis_on,1,'omitnan');
            se_on = std(data_vis_on,0,1,'omitnan') ./ sqrt(size(data_vis_on,1));
            [hl1,hl2] = boundedline(t_vis, m_on, se_on, ...
                'Color', colors_cond(cond,:), 'alpha');
            hl1.DisplayName = 'cs-on';
            hl2.HandleVisibility = 'off';
            legend_entries_vis{end+1} = sprintf('%s CS-on', cond_names{cond}); %#ok<AGROW>
        end

        % CS-on+180 (dashed)
        if ~isempty(data_vis_180)
            m_180  = mean(data_vis_180,1,'omitnan');
            se_180 = std(data_vis_180,0,1,'omitnan') ./ sqrt(size(data_vis_180,1));
            [hl1,hl2] = boundedline(t_vis, m_180, se_180, ...
                'Color', colors_cond(cond,:), 'alpha');
            hl1.DisplayName = 'cs-on+180';
            hl2.HandleVisibility = 'off';
            set(hl1, 'LineStyle', '--');
            legend_entries_vis{end+1} = sprintf('%s CS-on+180', cond_names{cond}); %#ok<AGROW>
        end
    end

    xline(0,'k--');
    xlabel('Time from target onset (ms)');
    ylabel('CS rate (sp/s)');
    title('Visual epoch: CS rate');
    if ~isempty(legend_entries_vis)
        legend(legend_entries_vis,'Box','off','Location','best');
    end

    % ---------- SACCADE CS RATE ----------
    nexttile(tl,2);
    hold on;

    legend_entries_sac = {};
    for cond = 1:n_conds
        % CS-on, saccade
        data_sac_on = squeeze(cs_rate_rot(:, t_ind_sac, dir_on, cond, 2));
        good_cells_on = ~all(isnan(data_sac_on),2);
        data_sac_on = data_sac_on(good_cells_on,:);

        % CS-on+180, saccade
        data_sac_180 = squeeze(cs_rate_rot(:, t_ind_sac, dir_180, cond, 2));
        good_cells_180 = ~all(isnan(data_sac_180),2);
        data_sac_180 = data_sac_180(good_cells_180,:);

        if isempty(data_sac_on) && isempty(data_sac_180)
            continue;
        end

        % CS-on
        if ~isempty(data_sac_on)
            m_on  = mean(data_sac_on,1,'omitnan');
            se_on = std(data_sac_on,0,1,'omitnan') ./ sqrt(size(data_sac_on,1));
            [hl1,hl2] = boundedline(t_sac, m_on, se_on, ...
                'Color', colors_cond(cond,:), 'alpha');
            hl1.DisplayName = 'cs-on';
            hl2.HandleVisibility = 'off';
            legend_entries_sac{end+1} = sprintf('%s CS-on', cond_names{cond}); %#ok<AGROW>
        end

        % CS-on+180 (dashed)
        if ~isempty(data_sac_180)
            m_180  = mean(data_sac_180,1,'omitnan');
            se_180 = std(data_sac_180,0,1,'omitnan') ./ sqrt(size(data_sac_180,1));
            [hl1,hl2] = boundedline(t_sac, m_180, se_180, ...
                'Color', colors_cond(cond,:), 'alpha');
            hl1.DisplayName = 'cs-on+180';
            hl2.HandleVisibility = 'off';
            set(hl1, 'LineStyle', '--');
            legend_entries_sac{end+1} = sprintf('%s CS-on+180', cond_names{cond}); %#ok<AGROW>
        end
    end

    xline(0,'k--');
    xlabel('Time from saccade onset (ms)');
    ylabel('CS rate (sp/s)');
    title('Saccade epoch: CS rate');
    if ~isempty(legend_entries_sac)
        legend(legend_entries_sac,'Box','off','Location','best');
    end

    % ---------- VISUAL EYE VELOCITY ----------
    nexttile(tl,3);
    hold on;

    for cond = 1:n_conds
        data_vm_vis_on = squeeze(cs_vm_rot(:, t_ind_vis, dir_on, cond, 1));
        good_cells_vm_on = ~all(isnan(data_vm_vis_on),2);
        data_vm_vis_on = data_vm_vis_on(good_cells_vm_on,:);

        data_vm_vis_180 = squeeze(cs_vm_rot(:, t_ind_vis, dir_180, cond, 1));
        good_cells_vm_180 = ~all(isnan(data_vm_vis_180),2);
        data_vm_vis_180 = data_vm_vis_180(good_cells_vm_180,:);

        if isempty(data_vm_vis_on) && isempty(data_vm_vis_180)
            continue;
        end

        % CS-on velocity
        if ~isempty(data_vm_vis_on)
            m_vm_on  = mean(data_vm_vis_on,1,'omitnan');
            se_vm_on = std(data_vm_vis_on,0,1,'omitnan') ./ sqrt(size(data_vm_vis_on,1));
            [hl1,hl2] = boundedline(t_vis, m_vm_on, se_vm_on, ...
                'Color', colors_cond(cond,:), 'alpha');
            hl1.DisplayName = sprintf('%s cs-on',cond_names{cond});
            hl2.HandleVisibility = 'off';
        end

        % CS-on+180 velocity (dashed)
        if ~isempty(data_vm_vis_180)
            m_vm_180  = mean(data_vm_vis_180,1,'omitnan');
            se_vm_180 = std(data_vm_vis_180,0,1,'omitnan') ./ sqrt(size(data_vm_vis_180,1));
            [hl1,hl2] = boundedline(t_vis, m_vm_180, se_vm_180, ...
                'Color', colors_cond(cond,:), 'alpha');
            set(hl1,'LineStyle','--');
            hl1.DisplayName = sprintf('%s cs-on+180',cond_names{cond});
            hl2.HandleVisibility = 'off';
        end
    end

    xline(0,'k--');
    xlabel('Time from target onset (ms)');
    ylabel('Eye speed (deg/s)');
    title('Visual epoch: eye velocity');

    % ---------- SACCADE EYE VELOCITY ----------
    nexttile(tl,4);
    hold on;

    for cond = 1:n_conds
        data_vm_sac_on = squeeze(cs_vm_rot(:, t_ind_sac, dir_on, cond, 2));
        good_cells_vm_on = ~all(isnan(data_vm_sac_on),2);
        data_vm_sac_on = data_vm_sac_on(good_cells_vm_on,:);

        data_vm_sac_180 = squeeze(cs_vm_rot(:, t_ind_sac, dir_180, cond, 2));
        good_cells_vm_180 = ~all(isnan(data_vm_sac_180),2);
        data_vm_sac_180 = data_vm_sac_180(good_cells_vm_180,:);

        if isempty(data_vm_sac_on) && isempty(data_vm_sac_180)
            continue;
        end

        % CS-on velocity
        if ~isempty(data_vm_sac_on)
            m_vm_on  = mean(data_vm_sac_on,1,'omitnan');
            se_vm_on = std(data_vm_sac_on,0,1,'omitnan') ./ sqrt(size(data_vm_sac_on,1));
            [hl1,hl2] = boundedline(t_sac, m_vm_on, se_vm_on, ...
                'Color', colors_cond(cond,:), 'alpha');
            hl1.DisplayName = sprintf('%s cs-on',cond_names{cond});
            hl2.HandleVisibility = 'off';
        end

        % CS-on+180 velocity (dashed)
        if ~isempty(data_vm_sac_180)
            m_vm_180  = mean(data_vm_sac_180,1,'omitnan');
            se_vm_180 = std(data_vm_sac_180,0,1,'omitnan') ./ sqrt(size(data_vm_sac_180,1));
            [hl1,hl2] = boundedline(t_sac, m_vm_180, se_vm_180, ...
                'Color', colors_cond(cond,:), 'alpha');
            set(hl1,'LineStyle','--');
            hl1.DisplayName = sprintf('%s cs-on+180',cond_names{cond});
            hl2.HandleVisibility = 'off';
        end
    end

    xline(0,'k--');
    xlabel('Time from saccade onset (ms)');
    ylabel('Eye speed (deg/s)');
    title('Saccade epoch: eye velocity');

   % ---------- REACTION TIME VIOLIN PLOT ----------

    nexttile(tl,[1 2]);  % span two columns
    hold on;
    
    dir_for_rt = dir_on;   % use CS-on direction for RT
    
    % Collect RT values per condition, then pack into a matrix
    rt_vals_all = cell(1,n_conds);
    max_n = 0;
    for cond = 1:n_conds
        vals = squeeze(cs_rt_rot(:, dir_for_rt, cond));   % [cell]
        vals = vals(~isnan(vals));                        % drop NaNs
        rt_vals_all{cond} = vals;
        max_n = max(max_n, numel(vals));
    end
    
    if max_n == 0
        % nothing to plot
        text(0.5,0.5,'No RT data','HorizontalAlignment','center');
        axis off;
    else
        % Build [max_n x n_conds] matrix, padded with NaNs
        data_mat = nan(max_n, n_conds);
        for cond = 1:n_conds
            v = rt_vals_all{cond};
            if ~isempty(v)
                data_mat(1:numel(v),cond) = v;
            end
        end
    
        % Call Bastian's violinplot: each column = one violin
        violins = violinplot(data_mat, cond_names, ...
            'ViolinColor', colors_cond, ...
            'ViolinAlpha', 0.3, ...
            'ShowData', true, ...
            'ShowBox', true, ...
            'ShowMedian', true);
    
       
        xlim([0.5, n_conds+0.5]);
        xticks(1:n_conds);
        xticklabels(cond_names);
        ylabel('Reaction time (ms)');
        title(sprintf('%s: Reaction time by stim type', sub_name));
        set(gca,'FontSize',12);
    end

    ESN_Beautify_Plot(fig,[12,10]);

    % optional save
    out_dir = fullfile(data_path,'population_figs','ephys');
    if ~exist(out_dir,'dir'); mkdir(out_dir); end
    print(fig, fullfile(out_dir, ...
        sprintf('CS_stimType_%s.pdf', sub_name)), ...
        '-dpdf','-bestfit');

end
