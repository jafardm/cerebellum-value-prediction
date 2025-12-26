function plot_cs_spe_rpe_tag_1(data_path)

 load(fullfile(data_path,'population_data','cs_on_rate_spe_rpe_prim_sac'));
 img_save_path = fullfile(data_path,'population_figs\ephys');
 CS_poulation_vis = squeeze(cs_rate(:,:,:,:,1)); % 1 visual, 
 VM_population_vis = squeeze(cs_vm(:,:,:,:,1));

 CS_poulation_sac = squeeze(cs_rate(:,:,:,:,2)); %  2 saccade
 VM_population_sac = squeeze(cs_vm(:,:,:,:,2));

 CS_poulation_off = squeeze(cs_rate(:,:,:,:,3)); %  3 offset
 VM_population_off = squeeze(cs_vm(:,:,:,:,3));

 num_cs = numel(cs_on_data);
 t_ind_vis = (-50:200) + 250;
 t_ind_sac = (-100:200) + 250;
 t_ind_off = (-100:250) + 250;
 % Define colors
  colors = [
       190, 37, 43;    % DarkRed HH
       255,  0,  230;  % pink    HL
       0, 10, 255;     % blue    LL
       0, 204, 255;    % cyan    LH
    ]/ 255;

 vel_diff_colors = [
     34, 139,  34;    % for pair #1 (e.g., HL-HH or HLJ-HHJ)
    218, 165,  32;    % for pair #2 (e.g., LH-LL or LHJ-LLJ)
 ]/255;

 colors_c = [colors(1,:); colors(3,:)]; % Choice condition colors
 vm_color = [204, 174, 98]/256;
 cond_labels = {'HH','HL','LL','LH'};
 choice_labels = {'High','Low'};

% -------------- plot visual no SPE / with RPE ---------

 fig = figure;
 set(fig, 'PaperUnits', 'inches');
 set(fig, 'PaperPosition', [0 0 6 4]); % [left bottom width height]
 set(fig, 'PaperSize', [6 4]);         % match the figure size
 hold on
 colororder({'k','k'})  % Left/right axis text color
 vm_lim = [0,400];

   for ii = 1:4
        currect_cond = mean(squeeze(CS_poulation_vis(:, t_ind_vis, :, ii)), 3, 'omitnan');
        Avg_cond = mean(currect_cond, 1, 'omitnan');
        std_cond = std(currect_cond, 1, 'omitnan') ./ sqrt(num_cs);
        yyaxis left  
        [hl1, hl2] = boundedline(t_ind_vis-250, Avg_cond, std_cond, 'Color', colors(ii,:), 'alpha');
        hl1.DisplayName = cond_labels{ii};
        hl2.HandleVisibility = 'off';
        if ii == 2 || ii == 4
           hl1.LineStyle = '--';
        end
       
        if ii == 1
            current_cs_vm = mean(mean(VM_population_vis(:, t_ind_vis, :, :), 4, 'omitnan'), 3, 'omitnan'); 
            yyaxis right
            area(t_ind_vis-250, mean(current_cs_vm,1), ...
                'FaceColor', vm_color, ...
                'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off')
            ylim(vm_lim);
        end
    end
    legend('Box','off','Location','northwest')
    xline(0, '--', 'HandleVisibility', 'off');
    hold off
    xlabel('Time relative to stimulus onset (ms)')
    ylabel('CS rate (Hz)');
    title(['no SPE / with RPE visual  N = ' num2str(num_cs)])
    ESN_Beautify_Plot(fig,[6,4])
    print(fig,fullfile(img_save_path,'cs-f-tag-1-nspe-rpe-vis.pdf'),'-dpdf', '-bestfit')

    fig=figure;
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperPosition', [0 0 6 4]); 
    set(fig, 'PaperSize', [6 4]);         
    hold on
    colororder({'k','k'})
    yyaxis left
    for cond_idx = 5:6
        current_cond = mean(squeeze(CS_poulation_vis(:, t_ind_vis, :, cond_idx)), 3, 'omitnan');
        Avg_cond = mean(current_cond,'omitnan');
        std_cond = std(current_cond, 1, 'omitnan') ./ sqrt(num_cs);
        [hl1,hl2]= boundedline(t_ind_vis-250, Avg_cond, std_cond, 'Color', colors_c(cond_idx-4,:), 'alpha');
        hl1.DisplayName = choice_labels{cond_idx-4};
        hl2.HandleVisibility = 'off';
    end
    ylabel('CS rate (Hz)')
    xline(0, '--', 'HandleVisibility','off');
    xlabel('Time from stimulus onset (ms)') 
    legend('Box','off','Location','northwest')
    title(['no SPE / with RPE visual  N = ' num2str(num_cs)])
    yyaxis right
    current_vm = mean(mean(squeeze(VM_population_vis(:, t_ind_vis, :, [5 6])),[3 4],'omitnan')); % VM for cond 5 6
    area(t_ind_vis-250, current_vm, 'FaceColor', vm_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility','off')
    ylim(vm_lim)
    ylabel('VM (deg/s)')
    title(['choice/ visual  N = ' num2str(num_cs)])
    ESN_Beautify_Plot(fig,[6,4])
    print(fig,fullfile(img_save_path,'cs-c-tag-1-nspe-rpe-vis.pdf'),'-dpdf', '-bestfit')

% -------------- plot saccade no SPE / with RPE ---------
    fig = figure;
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperPosition', [0 0 6 4]); 
    set(fig, 'PaperSize', [6 4]);   
    hold on
    colororder({'k','k'})  % Left/right axis text color

   for ii = 1:4
        currect_cond = mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, ii)), 3, 'omitnan');
        Avg_cond = mean(currect_cond, 1, 'omitnan');
        std_cond = std(currect_cond, 1, 'omitnan') ./ sqrt(num_cs);
        yyaxis left  
        [hl1, hl2] = boundedline(t_ind_sac-250, Avg_cond, std_cond, 'Color', colors(ii,:), 'alpha');
        hl1.DisplayName = cond_labels{ii};
        hl2.HandleVisibility = 'off';
        if ii == 2 || ii == 4
           hl1.LineStyle = '--';
        end
       
        if ii == 1
            current_cs_vm = mean(mean(VM_population_sac(:, t_ind_sac, :, 1:4), 4, 'omitnan'), 3, 'omitnan'); 
            yyaxis right
            area(t_ind_sac-250, mean(current_cs_vm,1), ...
                'FaceColor', vm_color, ...
                'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off')
            ylim(vm_lim);
        end
    end
    legend('Box','off','Location','northwest')
    xline(0, '--', 'HandleVisibility', 'off');
    hold off
    xlabel('Time relative to saccade onset (ms)')
    ylabel('CS rate (Hz)');
    title(['no SPE / with RPE saccade  N = ' num2str(num_cs)])
    ESN_Beautify_Plot(fig,[6,4])
    print(fig,fullfile(img_save_path,'cs-f-tag-1-nspe-rpe-sac.pdf'),'-dpdf', '-bestfit')

    fig = figure;
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperPosition', [0 0 6 4]); 
    set(fig, 'PaperSize', [6 4]);   
    hold on
    colororder({'k','k'})
    yyaxis left
    for cond_idx = 5:6
        current_cond = mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, cond_idx)), 3, 'omitnan');
        Avg_cond = mean(current_cond,'omitnan');
        std_cond = std(current_cond, 1, 'omitnan') ./ sqrt(num_cs);
        [hl1,hl2]= boundedline(t_ind_sac-250, Avg_cond, std_cond, 'Color', colors_c(cond_idx-4,:), 'alpha');
        hl1.DisplayName = choice_labels{cond_idx-4};
        hl2.HandleVisibility = 'off';
    end
    ylabel('CS rate (Hz)')
    xline(0, '--', 'HandleVisibility','off');
    xlabel('Time from saccade onset (ms)') 
    legend('Box','off','Location','northwest')
    xlim([-100, 200])

    yyaxis right
    current_vm = mean(mean(squeeze(VM_population_sac(:, t_ind_sac, :, [5,6])),[3 4],'omitnan')); 
    area(t_ind_sac-250, current_vm, 'FaceColor', vm_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility','off')
    ylim(vm_lim)
    ylabel('VM (deg/s)')
    title(['choice saccade  N = ' num2str(num_cs)])
    ESN_Beautify_Plot(fig,[6,4])
    print(fig,fullfile(img_save_path,'cs-c-tag-1-nspe-rpe-sac.pdf'),'-dpdf', '-bestfit')

% -------------- plot sac offset no SPE / with RPE ---------
 fig = figure;
 set(fig, 'PaperUnits', 'inches');
 set(fig, 'PaperPosition', [0 0 6 4]); 
 set(fig, 'PaperSize', [6 4]);   
 hold on
 colororder({'k','k'})  % Left/right axis text color

   for ii = 1:4
        currect_cond = mean(squeeze(CS_poulation_off(:, t_ind_off, :, ii)), 3, 'omitnan');
        Avg_cond = mean(currect_cond, 1, 'omitnan');
        std_cond = std(currect_cond, 1, 'omitnan') ./ sqrt(num_cs);
        yyaxis left  
        [hl1, hl2] = boundedline(t_ind_off-250, Avg_cond, std_cond, 'Color', colors(ii,:), 'alpha');
        hl1.DisplayName = cond_labels{ii};
        hl2.HandleVisibility = 'off';
        if ii == 2 || ii == 4
           hl1.LineStyle = '--';
        end
       
        if ii == 1
            current_cs_vm = mean(mean(VM_population_off(:, t_ind_off, :, 1:4), 4, 'omitnan'), 3, 'omitnan'); 
            yyaxis right
            area(t_ind_off-250, squeeze(mean(current_cs_vm,1)), ...
                'FaceColor', vm_color, ...
                'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off')
            ylim(vm_lim);
        end
    end
    legend('Box','off','Location','northwest')
    xline(0, '--', 'HandleVisibility', 'off');
    hold off
    xlabel('Time relative to saccade offset (ms)')
    ylabel('CS rate (Hz)');
    title(['no SPE / with RPE saccade offset N = ' num2str(num_cs)])
    ESN_Beautify_Plot(fig,[6,4])
    print(fig,fullfile(img_save_path,'cs-f-tag-1-nspe-rpe-offset.pdf'),'-dpdf', '-bestfit')

%---------------- subtraction of No SPE RPE sac offset-------------------------------------
      neg_rpe =  mean(squeeze(CS_poulation_off(:, t_ind_off, :, 2)), 3, 'omitnan')...
    - mean(squeeze(CS_poulation_off(:, t_ind_off, :, 1)), 3, 'omitnan');
      avg_neg_rpe = mean(neg_rpe,1,'omitnan');
      std_neg_rpe = std(neg_rpe, 1, 'omitnan') ./ sqrt(num_cs);
      pos_rpe =  mean(squeeze(CS_poulation_off(:, t_ind_off, :, 4)), 3, 'omitnan')...
    - mean(squeeze(CS_poulation_off(:, t_ind_off, :, 3)), 3, 'omitnan');
      avg_pos_rpe = mean(pos_rpe,1,'omitnan');
      std_pos_rpe = std(pos_rpe, 1, 'omitnan') ./ sqrt(num_cs);

    % VM differences (+ SEM)
    neg_vm = mean(squeeze(VM_population_off(:, t_ind_off, :, 2)), 3, 'omitnan') ...
           - mean(squeeze(VM_population_off(:, t_ind_off, :, 1)), 3, 'omitnan'); % HL - HH
    avg_neg_vm = mean(neg_vm, 1, 'omitnan');
    sem_neg_vm = std(neg_vm, 1, 'omitnan') ./ sqrt(num_cs);
    
    pos_vm = mean(squeeze(VM_population_off(:, t_ind_off, :, 4)), 3, 'omitnan') ...
           - mean(squeeze(VM_population_off(:, t_ind_off, :, 3)), 3, 'omitnan'); % LH - LL
    avg_pos_vm = mean(pos_vm, 1, 'omitnan');
    sem_pos_vm = std(pos_vm, 1, 'omitnan') ./ sqrt(num_cs);

    fig = figure; 
    hold on
    colororder({'k','k'})
    
    yyaxis left
    [hl1, hpf1] = boundedline(t_ind_off-250, avg_neg_rpe, std_neg_rpe, 'Color', colors(2,:), 'alpha'); %#ok<ASGLU>
    hl1.DisplayName = 'RPE- (CS)';
    hpf1.HandleVisibility = 'off';
    [hl2, hpf2] = boundedline(t_ind_off-250, avg_pos_rpe, std_pos_rpe, 'Color', colors(4,:), 'alpha');
    hl2.DisplayName = 'RPE+ (CS)';
    hpf2.HandleVisibility = 'off';
    ylabel('CS rate (Hz)');
    legend('Box','off','Location','northwest')
    
    yyaxis right
    % VM diffs with boundedline (no area)
    [hvm1, hvpf1] = boundedline(t_ind_off-250, avg_neg_vm, sem_neg_vm, 'Color', vel_diff_colors(1,:), 'alpha'); %#ok<ASGLU>
    hvm1.DisplayName = 'VM diff (RPE-)';
    hvpf1.HandleVisibility = 'off';
    [hvm2, hvpf2] = boundedline(t_ind_off-250, avg_pos_vm, sem_pos_vm, 'Color', vel_diff_colors(2,:), 'alpha');
    hvm2.DisplayName = 'VM diff (RPE+)';
    hvpf2.HandleVisibility = 'off';
    ylabel('VM diff (deg/s)');
    
    xline(0, '--', 'HandleVisibility','off');
    yline(0, '--', 'HandleVisibility','off');
    xlabel('Time from saccade offset (ms)');
    title(['no SPE / with RPE diff  N= ' num2str(num_cs)])
    ESN_Beautify_Plot(fig,[6,4])
    print(fig,fullfile(img_save_path,'cs-f-tag-1-nspe-RPE_diff-offset.pdf'),'-dpdf', '-bestfit')

end