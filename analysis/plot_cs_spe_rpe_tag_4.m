function plot_cs_spe_rpe_tag_4(data_path)
 
 load(fullfile(data_path,'population_data','cs_on_rate_spe_rpe_corr_sac'));
 img_save_path = fullfile(data_path,'population_figs\ephys'); 

 CS_poulation_vis  = squeeze(cs_rate(:,:,:,:,1)); % 1 visual, 
 VM_population_vis = squeeze(cs_vm(:,:,:,:,1));

 CS_poulation_sac  = squeeze(cs_rate(:,:,:,:,2)); %  2 saccade
 VM_population_sac = squeeze(cs_vm(:,:,:,:,2));

 num_cs = numel(cs_on_data);
 t_ind_vis = (-100:200) + 250;
 t_ind_sac = (-100:200) + 250;

 % Define colors
  colors = [
       190, 37, 43;    % DarkRed
       255,  0,  230;  % pink
       0, 10, 255;     % blue
       0, 204, 255;    % cyan
    ]/ 255;

 vm_color = [204, 174, 98]/256;
 cond_labels = {'HHJ','HLJ','LLJ','LHJ'};


% -------------- plot visual ---------

 fig = figure;
 subplot(1,3,1)
 hold on
 colororder({'k','k'})  % Left/right axis text color
 vm_lim = [0,400];

   for ii = 1:4
        currect_cond = mean(squeeze(CS_poulation_vis(:, t_ind_vis, :, ii)), 3, 'omitnan');
        Avg_cond = mean(currect_cond, 1, 'omitnan');
        std_cond = std(currect_cond, 1, 'omitnan') ./ sqrt(num_cs);
        yyaxis left  
        [hl1, hl2] = boundedline(t_ind_vis-250, Avg_cond, std_cond, 'Color', colors(ii,:));
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
    legend('Box','off','Location','northeast')
    xline(0, '--', 'HandleVisibility', 'off');
    hold off
    xlabel('Time relative to visual onset (ms)')
    ylabel('CS rate (Hz)');
    title('visual onset')


    %---------------- plot RPE visual onset -------------------------------------
          neg_rpe =  mean(squeeze(CS_poulation_vis(:, t_ind_vis, :, 2)), 3, 'omitnan')...
    - mean(squeeze(CS_poulation_vis(:, t_ind_vis, :, 1)), 3, 'omitnan'); % HL-HH
      avg_neg_rpe = mean(neg_rpe,1,'omitnan');
      std_neg_rpe = std(neg_rpe, 1, 'omitnan') ./ sqrt(num_cs);

      pos_rpe =  mean(squeeze(CS_poulation_vis(:, t_ind_vis, :, 4)), 3, 'omitnan')...
    - mean(squeeze(CS_poulation_vis(:, t_ind_vis, :, 3)), 3, 'omitnan'); % LH-LL
      avg_pos_rpe = mean(pos_rpe,1,'omitnan');
      std_pos_rpe = std(pos_rpe, 1, 'omitnan') ./ sqrt(num_cs);
     

    subplot(1,3,2)
    hold on
    colororder({'k','k'})
    yyaxis left
     [hl1, hl2] = boundedline(t_ind_vis-250, avg_neg_rpe, std_neg_rpe, 'Color', colors(2,:));
        hl1.DisplayName = 'RPE-(HL-HH)';
        hl2.HandleVisibility = 'off';
     [hl1, hl2] = boundedline(t_ind_vis-250, avg_pos_rpe, std_pos_rpe, 'Color', colors(4,:));
        hl1.DisplayName = 'RPE+(LH-LL)';
        hl2.HandleVisibility = 'off';   
    xline(0, '--', 'HandleVisibility','off');
    xlabel('Time relative to visual onset (ms)') 
    legend('Box','off','Location','northeast')
    title('RPE-diff (visual onset)')
    yyaxis right
    current_vm = mean(mean(VM_population_vis(:, t_ind_vis, :, :), 4, 'omitnan'), 3, 'omitnan'); 
    area(t_ind_vis-250, squeeze(mean(current_vm,1)), 'FaceColor',...
        vm_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility','off')
    ylim(vm_lim)
      
  %---------------- plot LH-HH, HL-LL (visual) -------------------------------------
    LH_minus_HH =  mean(squeeze(CS_poulation_vis(:, t_ind_vis, :, 4)), 3, 'omitnan') ...
                 - mean(squeeze(CS_poulation_vis(:, t_ind_vis, :, 1)), 3, 'omitnan');  % 4-1
    avg_LH_minus_HH = mean(LH_minus_HH, 1, 'omitnan');
    std_LH_minus_HH = std(LH_minus_HH, 1, 'omitnan') ./ sqrt(num_cs);
    
    HL_minus_LL =  mean(squeeze(CS_poulation_vis(:, t_ind_vis, :, 2)), 3, 'omitnan') ...
                 - mean(squeeze(CS_poulation_vis(:, t_ind_vis, :, 3)), 3, 'omitnan');  % 2-3
    avg_HL_minus_LL = mean(HL_minus_LL, 1, 'omitnan');
    std_HL_minus_LL = std(HL_minus_LL, 1, 'omitnan') ./ sqrt(num_cs);
    
    subplot(1,3,3)
    hold on
    colororder({'k','k'})
    yyaxis left
    [hl1, hl2] = boundedline(t_ind_vis-250, avg_HL_minus_LL, std_HL_minus_LL, 'Color', colors(2,:)); % HL color
    hl1.DisplayName = 'HL-LL';
    hl2.HandleVisibility = 'off';
    
    [hl1, hl2] = boundedline(t_ind_vis-250, avg_LH_minus_HH, std_LH_minus_HH, 'Color', colors(4,:)); % LH color
    hl1.DisplayName = 'LH-HH';
    hl2.HandleVisibility = 'off';
    
    xline(0, '--', 'HandleVisibility','off');
    xlabel('Time relative to visual onset (ms)')
    legend('Box','off','Location','northeast')
    
    yyaxis right
    current_vm = mean(mean(VM_population_vis(:, t_ind_vis, :, :), 4, 'omitnan'), 3, 'omitnan');
    area(t_ind_vis-250, squeeze(mean(current_vm,1)), 'FaceColor', vm_color, ...
        'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility','off')
    ylim(vm_lim)
    ylabel('VM (deg/s)')
    title('HL-LL and LH-HH (visual onset)')

    ylabel('VM (deg/s)')
    title('HL-LL and LH-HH (visual onset)')
    sgtitle(sprintf('Corrective saccdes, CS Cells (n=%d)', num_cs));
    ESN_Beautify_Plot(fig,[18,6])
    print(fig,fullfile(img_save_path,'cs-f-tag-4-spe-rpe-vis.pdf'),'-dpdf', '-bestfit')

%% -------------- plot saccade  ---------
    fig2=figure;
    subplot(1,3,1)
    hold on
    colororder({'k','k'})  % Left/right axis text color
    vm_lim = [0,400];
   for ii = 1:4
        currect_cond = mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, ii)), 3, 'omitnan');
        Avg_cond = mean(currect_cond, 1, 'omitnan');
        std_cond = std(currect_cond, 1, 'omitnan') ./ sqrt(num_cs);
        yyaxis left  
        [hl1, hl2] = boundedline(t_ind_sac-250, Avg_cond, std_cond, 'Color', colors(ii,:));
        hl1.DisplayName = cond_labels{ii};
        hl2.HandleVisibility = 'off';
        if ii == 2 || ii == 4
           hl1.LineStyle = '--';
        end
       
        if ii == 1
            current_cs_vm = mean(mean(VM_population_sac(:, t_ind_sac, :, :), 4, 'omitnan'), 3, 'omitnan'); 
            yyaxis right
            area(t_ind_sac-250, mean(current_cs_vm,1), ...
                'FaceColor', vm_color, ...
                'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off')
            ylim(vm_lim);
        end
    end
    legend('Box','off','Location','northeast')
    xline(0, '--', 'HandleVisibility', 'off');
    hold off
    xlabel('Time relative to saccade onset (ms)')
    ylabel('CS rate (Hz)');
    title('saccade onset')
      
% -------------- plot RPE ---------

      neg_rpe =  mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, 2)), 3, 'omitnan')...
    - mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, 1)), 3, 'omitnan');
      avg_neg_rpe = mean(neg_rpe,1,'omitnan');
      std_neg_rpe = std(neg_rpe, 1, 'omitnan') ./ sqrt(num_cs);

      pos_rpe =  mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, 4)), 3, 'omitnan')...
    - mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, 3)), 3, 'omitnan');
      avg_pos_rpe = mean(pos_rpe,1,'omitnan');
      std_pos_rpe = std(pos_rpe, 1, 'omitnan') ./ sqrt(num_cs);
    

    subplot(1,3,2)
    hold on
    colororder({'k','k'})
    yyaxis left
     [hl1, hl2] = boundedline(t_ind_sac-250, avg_neg_rpe, std_neg_rpe, 'Color', colors(2,:));
        hl1.DisplayName = 'RPE-(HL-HH)';
        hl2.HandleVisibility = 'off';
     [hl1, hl2] = boundedline(t_ind_sac-250, avg_pos_rpe, std_pos_rpe, 'Color', colors(4,:));
        hl1.DisplayName = 'RPE+(LH-LL)';
        hl2.HandleVisibility = 'off';   
    xline(0, '--', 'HandleVisibility','off');
    xlabel('Time relative to saccade onset (ms)') 
    legend('Box','off','Location','northeast')

    yyaxis right
    current_vm = mean(mean(VM_population_sac(:, t_ind_sac, :, :), 4, 'omitnan'), 3, 'omitnan'); 
    area(t_ind_sac-250, squeeze(mean(current_vm,1)), 'FaceColor',...
        vm_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility','off')
    title('RPE-diff (saccade onset)')
   
    
      %---------------- plot LH-HH, HL-LL (saccade) -----------------------------------
    LH_minus_HH =  mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, 4)), 3, 'omitnan') ...
                 - mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, 1)), 3, 'omitnan');  % 4-1
    avg_LH_minus_HH = mean(LH_minus_HH, 1, 'omitnan');
    std_LH_minus_HH = std(LH_minus_HH, 1, 'omitnan') ./ sqrt(num_cs);
    
    HL_minus_LL =  mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, 2)), 3, 'omitnan') ...
                 - mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, 3)), 3, 'omitnan');  % 2-3
    avg_HL_minus_LL = mean(HL_minus_LL, 1, 'omitnan');
    std_HL_minus_LL = std(HL_minus_LL, 1, 'omitnan') ./ sqrt(num_cs);
    
    subplot(1,3,3)
    hold on
    colororder({'k','k'})
    yyaxis left
    [hl1, hl2] = boundedline(t_ind_sac-250, avg_HL_minus_LL, std_HL_minus_LL, 'Color', colors(2,:)); % HL color
    hl1.DisplayName = 'HL-LL';
    hl2.HandleVisibility = 'off';
    
    [hl1, hl2] = boundedline(t_ind_sac-250, avg_LH_minus_HH, std_LH_minus_HH, 'Color', colors(4,:)); % LH color
    hl1.DisplayName = 'LH-HH';
    hl2.HandleVisibility = 'off';
    
    xline(0, '--', 'HandleVisibility','off');
    xlabel('Time relative to saccade onset (ms)')
    legend('Box','off','Location','northeast')
    
    yyaxis right
    current_vm = mean(mean(VM_population_sac(:, t_ind_sac, :, :), 4, 'omitnan'), 3, 'omitnan');
    area(t_ind_sac-250, squeeze(mean(current_vm,1)), 'FaceColor', vm_color, ...
        'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility','off')
    ylim(vm_lim)
    ylabel('VM (deg/s)')
    title('HL-LL and LH-HH (saccade)')

    sgtitle(sprintf('corrective saccades CS Cells (n=%d)', num_cs));
    ESN_Beautify_Plot(fig2,[18,6])
    print(fig2,fullfile(img_save_path,'cs-f-tag-4-spe-rpe-sac.pdf'),'-dpdf', '-bestfit')
    

end   
