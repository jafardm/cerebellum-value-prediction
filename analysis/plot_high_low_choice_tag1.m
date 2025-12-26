function plot_high_low_choice_tag1(data_path)

    % Load parameters and data
    MAF_params_funcs
    load(fullfile(data_path, 'population_data', 'cs_on_data'));

    % Extract visual and motor rho and angles
    cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_data, 'UniformOutput', false));
    cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_data, 'UniformOutput', false));
    cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, cs_on_data, 'UniformOutput', false));
    cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, cs_on_data, 'UniformOutput', false));
    
    num_cs = numel(cs_on_data);

    % Choose direction based on stronger tuning (rho)
    ind_m = cs_on_rho_sac > cs_on_rho_vis;
    cs_on_ang = deg2rad(cs_on_ang_vis);
    cs_on_ang(ind_m) = deg2rad(cs_on_ang_sac(ind_m));
    
    cs_on_rho = cs_on_rho_vis;
    cs_on_rho(ind_m) = cs_on_rho_sac(ind_m);

    % Bin angles to nearest preferred direction
    ang_bins = deg2rad(params.sac.ang_values);
    [~, cs_on_bin] = min(abs(angdiff(repmat(ang_bins, num_cs, 1), cs_on_ang .* ones(size(ang_bins)))), [], 2);

    % Align directions (rotating so CS-on is center)
    order_ang = nan(num_cs, 8);
    for counter_cs = 1:num_cs
        order_ang(counter_cs, :) = circshift(1:8, 5 - cs_on_bin(counter_cs));
    end

    % Extract and rotate average firing rates
    cs_on_rate_vis = cell2mat(cellfun(@(x) x.vis.fr_avg, cs_on_data, 'UniformOutput', false));
    cs_on_rate_sac = cell2mat(cellfun(@(x) x.sac.fr_avg, cs_on_data, 'UniformOutput', false));
    cs_on_rate_vis_rot = nan(num_cs, 9);
    cs_on_rate_sac_rot = nan(num_cs, 9);
    cs_rate_rot = nan(num_cs, 500, 8, 10, 2);  % time x dir x cond x epoch
    cs_vm_rot = nan(num_cs, 500, 8, 10, 2);

    % Repeat first direction to close the loop
    order_ang = [order_ang, order_ang(:,1)];

    for counter_cs = 1:num_cs
        cs_on_rate_vis_rot(counter_cs, :) = cs_on_rate_vis(counter_cs, order_ang(counter_cs,:));
        cs_on_rate_sac_rot(counter_cs, :) = cs_on_rate_sac(counter_cs, order_ang(counter_cs,:));
        cs_rate_rot(counter_cs, :, :, :, :) = cs_rate(counter_cs, :, order_ang(counter_cs,1:8), :, :);
        cs_vm_rot(counter_cs, :, :, :, :) = cs_vm(counter_cs, :, order_ang(counter_cs,1:8), :, :);
    end

    % Setup plotting
    choice_cond = params.reward_choice_conds;
    cond_labels = {'High','Low'};
    leg_ = {'cs-on','cs±45','cs±90','cs±135','cs±180'};
    colors_ = cool(5);
    dirs = 1:5;
    t_ind_vis = (-100:250) + 250;
    t_ind_sac = (-200:250) + 250;
    vm_color = [204, 174, 98]/256;
    vm_lim = [0,800];
    for counter_r = (1:size(choice_cond, 1)) + 8
        fig = figure;
        for counter_dir = numel(dirs):-1:1
            current_dir = dirs(counter_dir);
            if current_dir == 1 || current_dir == 5
                current_dir_ = current_dir;
            else
                current_dir_ = [current_dir, 10 - current_dir];
            end

            % Extract and average data across directions for each cell
            % Visual-aligned data
            current_cs_rate_vis_all = cs_rate_rot(:, t_ind_vis, current_dir_, counter_r, 1);
            current_cs_rate_vis_per_cell = squeeze(mean(current_cs_rate_vis_all, 3, 'omitnan'));
            mean_vis = mean(current_cs_rate_vis_per_cell, 1, 'omitnan');
            std_vis = std(current_cs_rate_vis_per_cell, 0, 1, 'omitnan') / sqrt(num_cs);

            % Saccade-aligned data
            current_cs_rate_sac_all = cs_rate_rot(:, t_ind_sac, current_dir_, counter_r, 2);
            current_cs_rate_sac_per_cell = squeeze(mean(current_cs_rate_sac_all, 3, 'omitnan'));
            mean_sac = mean(current_cs_rate_sac_per_cell, 1, 'omitnan');
            std_sac = std(current_cs_rate_sac_per_cell, 0, 1, 'omitnan') / sqrt(num_cs);

            if current_dir == 1
                current_cs_vm_sac = mean(mean(cs_vm_rot(:, t_ind_sac, :, :, 2), [3 4], 'omitnan'));
                current_cs_vm_vis = mean(mean(cs_vm_rot(:, t_ind_vis, :, :, 1), [3 4], 'omitnan'));
            end

            % === Visual-aligned Plot ===
            subplot(1,2,1);
            colororder({'k','k'})
            yyaxis left
            [hl1,hl2] = boundedline(t_ind_vis - 250, mean_vis, std_vis, 'Color', colors_(counter_dir,:));
            hl1.DisplayName = leg_{6 - counter_dir};
            hl2.HandleVisibility = 'off';
            hold on;
            xlabel('time from cue onset (ms)');
            ylabel(' CS rate (Hz)');
            legend('Box','off');
            xline(0, '--', 'HandleVisibility', 'off');

            if current_dir == 1
                yyaxis right
                area(t_ind_vis - 250, squeeze(mean(current_cs_vm_vis,1)), ...
                    'FaceColor', vm_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                ylim(vm_lim);
            end

            % === Saccade-aligned Plot ===
            subplot(1,2,2);
            colororder({'k','k'})
            yyaxis left
            [hl1,hl2] = boundedline(t_ind_sac - 250, mean_sac, std_sac, 'Color', colors_(counter_dir,:));
            hl1.DisplayName = leg_{6 - counter_dir};
            hl2.HandleVisibility = 'off';
            hold on;
            xlabel('time from sac onset (ms)');
            ylabel(' CS rate (Hz)');
            xline(0, '--', 'HandleVisibility', 'off');

            if current_dir == 1
                yyaxis right
                area(t_ind_sac - 250, squeeze(mean(current_cs_vm_sac,1)), ...
                    'FaceColor', vm_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                ylim(vm_lim);
                %plot([-70,30],[2.4,2.4], '-k', 'LineWidth', 2, 'HandleVisibility', 'off');
            end
        end        
 
        ESN_Beautify_Plot(fig,[10,4])  % Custom plot formatting
        sgtitle(['Choice Sac. : ' cond_labels{counter_r-8}], 'FontSize', 12, 'FontWeight', 'bold');
    end
end