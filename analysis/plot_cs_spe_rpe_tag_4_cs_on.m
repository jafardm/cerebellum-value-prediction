function plot_cs_spe_rpe_tag_4_cs_on(data_path)
 JDM_params_funcs
 S = load(fullfile(data_path,'population_data','cs_on_rate_spe_rpe_corr_sac'));
 img_save_path = fullfile(data_path,'population_figs\ephys');
 cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg,S.cs_on_data,...
        'UniformOutput',false));
    cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg,S.cs_on_data,...
        'UniformOutput',false));
    
    cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg,S.cs_on_data,...
        'UniformOutput',false));
    cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg,S.cs_on_data,...
        'UniformOutput',false));
    
    num_cs = numel(S.cs_on_data);
    
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
        S.cs_on_data,'UniformOutput',false));
    cs_on_rate_sac = cell2mat(cellfun(@(x) x.sac.fr_avg,...
        S.cs_on_data,'UniformOutput',false));
    
    % organizing the tuning of the cs rates
    cs_on_rate_vis_rot = nan(num_cs,9);
    cs_on_rate_sac_rot = nan(num_cs,9);
    cs_rate_rot        = nan(num_cs,500,8,4,2);
    cs_vm_rot          = nan(num_cs,500,8,4,2);
 
    for counter_cs = 1:num_cs
        cs_on_rate_vis_rot(counter_cs,:) =...
            cs_on_rate_vis(counter_cs,order_ang(counter_cs,:));
        cs_on_rate_sac_rot(counter_cs,:) =...
            cs_on_rate_sac(counter_cs,order_ang(counter_cs,:));
        cs_rate_rot(counter_cs,:,:,:,:)    =...
            S.cs_rate(counter_cs,:,order_ang(counter_cs,1:8),:,:); % picked visual and sac dim
        cs_vm_rot(counter_cs,:,:,:,:)    =...
            S.cs_vm(counter_cs,:,order_ang(counter_cs,1:8),:,:);
    end

    CENTER = 250;
    t_vis  = (-100:200) + CENTER;   x_vis = t_vis - CENTER;
    t_sac  = (-100:200) + CENTER;   x_sac = t_sac - CENTER;
   
    
    E_VIS = 1; E_SAC = 2; 
    cs_on_dir = 5;          % compute_CS_ON rotates so CS-on is index 5
    cs_off_dir = 1;
    cs_off    = 1;
    forced = 1:4;           % HH,HL,LL,LH
    cond_labels = {'HHJ','HLJ','LLJ','LHJ'};
    vm_lim = [0 400];
    y_lim = [0 4];
    colors = [ ...
       190,  37,  43;  % HH
       255,   0, 230;  % HL
         0,  10, 255;  % LL
         0, 204, 255] / 255; % LH
    vm_color = [204, 174, 98]/256;


    fig = figure;
    % ---------- (1) VISUAL ----------
    subplot(1,2,1); hold on; colororder({'k','k'});
    
    for ii = forced
        [avg, sem] = mean_sem_cells(squeeze(cs_rate_rot(:, t_vis, cs_on_dir, ii, E_VIS)));
        yyaxis left
        [hl, hp] = boundedline(x_vis, avg, sem, 'Color', colors(ii,:), 'alpha'); %#ok<ASGLU>
        hl.DisplayName = cond_labels{ii}; hp.HandleVisibility = 'off';
        if ii==2 || ii==4, hl.LineStyle='--'; end
        ylim(y_lim)
    end

    yyaxis right
    vm_vis = squeeze(mean(cs_vm_rot(:, t_vis, cs_on_dir, forced, E_VIS), 4, 'omitnan')); % [cells x T]
    vm_vis = mean(vm_vis, 1, 'omitnan');
    area(x_vis, vm_vis, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
        'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim); 
    yyaxis left; ylabel('CS rate (Hz)');
    xline(0,'--','HandleVisibility','off'); title('Visual (CS-on)'); legend('Box','off','Location','northwest');

    % ---------- (2) SACCADE ----------
    subplot(1,2,2); hold on; colororder({'k','k'});
    for ii = forced
        [avg, sem] = mean_sem_cells(squeeze(cs_rate_rot(:, t_sac, cs_on_dir, ii, E_SAC)));
        yyaxis left
        [hl, hp] = boundedline(x_sac, avg, sem, 'Color', colors(ii,:), 'alpha'); %#ok<ASGLU>
        hl.DisplayName = cond_labels{ii}; hp.HandleVisibility = 'off';
        if ii==2 || ii==4, hl.LineStyle='--'; end
        ylim(y_lim)
    end
    yyaxis right
    vm_sac = squeeze(mean(cs_vm_rot(:, t_sac, cs_on_dir, forced, E_SAC), 4, 'omitnan')); % [cells x T]
    vm_sac = mean(vm_sac, 1, 'omitnan');
    area(x_sac, vm_sac, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
        'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim); ylabel('VM (deg/s)');
    yyaxis left;
    xline(0,'--','HandleVisibility','off'); title('Saccade (CS-on)');

    sgtitle(sprintf('CS-on direction — HHJ / HLJ / LLJ / LHJ   (N = %d cells)', num_cs));
    ESN_Beautify_Plot(fig,[10,4]);
    print(fig, fullfile(img_save_path,'Figure2B_theta.pdf'), '-dpdf', '-bestfit');

    % ---------- saccade ----------
    align_dir = 2;
    avg_all = cell(1,4); sem_all = cell(1,4);
    for ii = forced
        [avg, sem] = mean_sem_cells(squeeze(cs_rate_rot(:, t_sac, cs_on_dir, ii, align_dir)));
        avg_all{ii} = avg; sem_all{ii} = sem;
   
    end

        % ---------- (3) SACCADE — tag=4 at CS+180 ----------


    fig2 = figure;
    subplot(1,2,1); hold on; colororder({'k','k'});
    % --- Visual (CS+180) ---
    for ii = forced
        [avg, sem] = mean_sem_cells(squeeze(cs_rate_rot(:, t_vis, cs_off_dir, ii, E_VIS)));
        yyaxis left
        [hl, hp] = boundedline(x_vis, avg, sem, 'Color', colors(ii,:), 'alpha'); %#ok<ASGLU>
        hl.DisplayName = cond_labels{ii}; hp.HandleVisibility = 'off';
        if ii==2 || ii==4, hl.LineStyle='--'; end
        ylim([0 2])
    end
    yyaxis right
    vm_vis = squeeze(mean(cs_vm_rot(:, t_vis, cs_off_dir, forced, E_VIS), 4, 'omitnan'));
    vm_vis = mean(vm_vis, 1, 'omitnan');
    area(x_vis, vm_vis, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
        'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim);
    yyaxis left; ylabel('CS rate (Hz)');
    xline(0,'--','HandleVisibility','off');
    xlabel('Time from target onset (ms)')
    title('Visual (CS+180)');
    legend('Box','off','Location','northeast');

    % --- Saccade (CS+180) ---
    subplot(1,2,2); hold on; colororder({'k','k'});
    for ii = forced
        [avg, sem] = mean_sem_cells(squeeze(cs_rate_rot(:, t_sac, cs_off_dir, ii, E_SAC)));
        yyaxis left
        [hl, hp] = boundedline(x_sac, avg, sem, 'Color', colors(ii,:), 'alpha'); %#ok<ASGLU>
        hl.DisplayName = cond_labels{ii}; hp.HandleVisibility = 'off';
        if ii==2 || ii==4, hl.LineStyle='--'; end
        ylim([0 2])
    end
    yyaxis right
    vm_sac = squeeze(mean(cs_vm_rot(:, t_sac, cs_off_dir, forced, E_SAC), 4, 'omitnan'));
    vm_sac = mean(vm_sac, 1, 'omitnan');
    area(x_sac, vm_sac, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
        'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim); ylabel('VM (deg/s)');
    yyaxis left;
    xline(0,'--','HandleVisibility','off');
    xlabel('Time from saccade onset (ms)')
    title('Saccade (CS+180)');

    sgtitle(sprintf('CS+180 direction — HHJ / HLJ / LLJ / LHJ   (N = %d cells)', num_cs));
    ESN_Beautify_Plot(fig2,[10,4]);
    print(fig2, fullfile(img_save_path,'Figure2_theta_plus_pi.pdf'), '-dpdf', '-bestfit');

    
   % --- RPE differences at SACCADE CS-on direction---
        fig = figure; hold on; colororder({'k','k'});
        subplot(1,2,1)
        
        % Cell-by-cell matrices
        M_LH = squeeze(cs_rate_rot(:, t_sac, cs_on_dir, 4, E_SAC)); % LH
        M_HH = squeeze(cs_rate_rot(:, t_sac, cs_on_dir, 1, E_SAC)); % HH
        M_HL = squeeze(cs_rate_rot(:, t_sac, cs_on_dir, 2, E_SAC)); % HL
        M_LL = squeeze(cs_rate_rot(:, t_sac, cs_on_dir, 3, E_SAC)); % LL
        
        % Cell-level differences
        M_rpe_plus  = M_LH - M_HH; % LH - HH
        M_rpe_minus = M_HL - M_LL; % HL - LL
        
        % Mean ± SEM across cells
        [avg_plus, sem_plus]   = mean_sem_cells(M_rpe_plus);
        [avg_minus, sem_minus] = mean_sem_cells(M_rpe_minus);
        
        yyaxis left
        [hl1, hp1] = boundedline(x_sac, avg_plus, sem_plus, ...
            'cmap',[0 0.6 0], 'alpha'); %#ok<ASGLU>
        hl1.DisplayName = 'RPE+ (LH-HH)'; hp1.HandleVisibility = 'off';
        
        [hl2, hp2] = boundedline(x_sac, avg_minus, sem_minus, ...
            'cmap',[0.6 0.3 0], 'alpha'); %#ok<ASGLU>
        hl2.DisplayName = 'RPE- (HL-LL)'; hp2.HandleVisibility = 'off';
        
        yyaxis right
        VM_ = squeeze(mean(cs_vm_rot(:, t_sac, cs_on_dir, forced, E_SAC), 4, 'omitnan')); % [cells x T]
        VM_ = mean(VM_, 1, 'omitnan');
        area(x_sac, VM_, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
            'EdgeColor','none','HandleVisibility','off');
        ylim(vm_lim); ylabel('VM (deg/s)');
        
        yyaxis left; ylabel('CS rate (Hz)');
        xline(0,'--','HandleVisibility','off'); 
        yline(0,'--','HandleVisibility','off'); 
        title('RPE differences — Corrective Saccade (CS-on)'); xlabel('Time (ms)');
        legend('Box','off','Location','northeast');
        

        % --- RPE differences at SACCADE CS-off (cs+180) direction---
       subplot(1,2,2), hold on; colororder({'k','k'});
        
        % Cell-by-cell matrices
        M_LH = squeeze(cs_rate_rot(:, t_sac, cs_off, 4, E_SAC)); % LH
        M_HH = squeeze(cs_rate_rot(:, t_sac, cs_off, 1, E_SAC)); % HH
        M_HL = squeeze(cs_rate_rot(:, t_sac, cs_off, 2, E_SAC)); % HL
        M_LL = squeeze(cs_rate_rot(:, t_sac, cs_off, 3, E_SAC)); % LL
        
        % Cell-level differences
        M_rpe_plus  = M_LH - M_HH; % LH - HH
        M_rpe_minus = M_HL - M_LL; % HL - LL
        
        % Mean ± SEM across cells
        [avg_plus, sem_plus]   = mean_sem_cells(M_rpe_plus);
        [avg_minus, sem_minus] = mean_sem_cells(M_rpe_minus);
        
        yyaxis left
        [hl1, hp1] = boundedline(x_sac, avg_plus, sem_plus, ...
            'cmap',[0 0.6 0], 'alpha'); 
        hl1.DisplayName = 'RPE+ (LH-HH)'; hp1.HandleVisibility = 'off';
        
        [hl2, hp2] = boundedline(x_sac, avg_minus, sem_minus, ...
            'cmap',[0.6 0.3 0], 'alpha'); 
        hl2.DisplayName = 'RPE- (HL-LL)'; hp2.HandleVisibility = 'off';
        
        yyaxis right
        VM_ = squeeze(mean(cs_vm_rot(:, t_sac, cs_off, forced, E_SAC), 4, 'omitnan')); % [cells x T]
        VM_ = mean(VM_, 1, 'omitnan');
        area(x_sac, VM_, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
            'EdgeColor','none','HandleVisibility','off');
        ylim(vm_lim); ylabel('VM (deg/s)');
        
        yyaxis left; ylabel('CS rate (Hz)');
        xline(0,'--','HandleVisibility','off'); 
        yline(0,'--','HandleVisibility','off'); 
        title('RPE differences — Corrective Saccade (\theta+180)'); xlabel('Time (ms)');
        legend('Box','off','Location','northeast');
        
        ESN_Beautify_Plot(fig,[10,4]);
        print(fig, fullfile(img_save_path,'Figure2B-RPE-difference-cs-off.pdf'), '-dpdf','-bestfit');

        %% ---------- HH and LL differences (CS-on minus CS+180) ----------
   fig_new = figure; hold on; colororder({'k','k'});

    % Extract HH+LL combined (SPE) matrices at VISUAL epoch
    spe_on  = mean(squeeze(cs_rate_rot(:, t_vis, cs_on_dir, [1 3], E_VIS)), 3, 'omitnan');   % CS-on
    spe_off = mean(squeeze(cs_rate_rot(:, t_vis, cs_off_dir, [1 3], E_VIS)), 3, 'omitnan');  % CS+180
    
    % Cell-by-cell differences
    spe_diff = spe_on - spe_off;   % combined SPE diff: CS-on – CS+180
    
    % Mean ± SEM
    [avg_spe, sem_spe] = mean_sem_cells(spe_diff);
    
    % ---- Plot combined SPE difference ----
    [hl1, hp1] = boundedline(x_vis, avg_spe, sem_spe, ...
        'cmap', [190 37 43]/255, 'alpha');
    hl1.DisplayName = 'SPE: CS-on – CS+180';
    hp1.HandleVisibility = 'off';
    
    xline(0,'--','HandleVisibility','off');
    xlabel('Time from target onset (ms)')
    ylabel('CS rate (Hz)')
    
    % ---- Velocity trace ----
    yyaxis right
    VM_ = squeeze(mean(cs_vm_rot(:, t_vis, cs_on_dir, [1 3], E_VIS), 4, 'omitnan'));
    VM_ = mean(VM_, 1, 'omitnan');
    
    area(x_vis, VM_, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
        'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim);
    ylabel('VM (deg/s)');
    
    ESN_Beautify_Plot(fig_new,[8,4]);
    print(fig_new, fullfile(img_save_path,'CS-SPE-diff-CSon-vs-CS180.pdf'),'-dpdf','-bestfit');

end % end function


% -------- helper --------
function [avg, sem] = mean_sem_cells(M)
% M: [cells x time] 
    avg = mean(M, 1, 'omitnan');
    N   = sum(~isnan(M), 1);
    sem = std(M, 0, 1, 'omitnan') ./ max(sqrt(N), 1);
end