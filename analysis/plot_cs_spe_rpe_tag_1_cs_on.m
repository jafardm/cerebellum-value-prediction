function plot_cs_spe_rpe_tag_1_cs_on(data_path)
 JDM_params_funcs
 S = load(fullfile(data_path,'population_data','cs_on_rate_spe_rpe_prim_sac'));
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
    cs_rate_rot        = nan(num_cs,500,8,6,3);
    cs_vm_rot          = nan(num_cs,500,8,6,3);
 
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
    t_vis  = (-50:200)  + CENTER;   x_vis = t_vis - CENTER;
    t_sac  = (-100:200) + CENTER;   x_sac = t_sac - CENTER;
    t_off  = (-100:250) + CENTER;   x_off = t_off - CENTER;
    
    E_VIS = 1; E_SAC = 2; E_OFF = 3;
    cs_on_dir = 5;          % compute_CS_ON rotates so CS-on is index 5
    forced = 1:4;           % HH,HL,LL,LH
    cond_labels = {'HH','HL','LL','LH'};
    vm_lim = [0 600];
    
    colors = [ ...
       190,  37,  43;  % HH
       255,   0, 230;  % HL
         0,  10, 255;  % LL
         0, 204, 255] / 255; % LH
    vm_color = [204, 174, 98]/256;


    fig = figure;
    % ---------- (1) VISUAL ----------
    subplot(3,1,1); hold on; colororder({'k','k'});
    for ii = forced
        [avg, sem] = mean_sem_cells(squeeze(cs_rate_rot(:, t_vis, cs_on_dir, ii, E_VIS)));
        yyaxis left
        [hl, hp] = boundedline(x_vis, avg, sem, 'Color', colors(ii,:), 'alpha'); %#ok<ASGLU>
        hl.DisplayName = cond_labels{ii}; hp.HandleVisibility = 'off';
        if ii==2 || ii==4, hl.LineStyle='--'; end
    end

    yyaxis right
    vm_vis = squeeze(mean(cs_vm_rot(:, t_vis, cs_on_dir, forced, E_VIS), 4, 'omitnan')); % [cells x T]
    vm_vis = mean(vm_vis, 1, 'omitnan');
    area(x_vis, vm_vis, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
        'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim); ylabel('VM (deg/s)');
    yyaxis left; ylabel('CS rate (Hz)');
    xline(0,'--','HandleVisibility','off'); title('Visual (CS-on)'); legend('Box','off','Location','northwest');

    % ---------- (2) SACCADE ----------
    subplot(3,1,2); hold on; colororder({'k','k'});
    for ii = forced
        [avg, sem] = mean_sem_cells(squeeze(cs_rate_rot(:, t_sac, cs_on_dir, ii, E_SAC)));
        yyaxis left
        [hl, hp] = boundedline(x_sac, avg, sem, 'Color', colors(ii,:), 'alpha'); %#ok<ASGLU>
        hl.DisplayName = cond_labels{ii}; hp.HandleVisibility = 'off';
        if ii==2 || ii==4, hl.LineStyle='--'; end
    end
    yyaxis right
    vm_sac = squeeze(mean(cs_vm_rot(:, t_sac, cs_on_dir, forced, E_SAC), 4, 'omitnan')); % [cells x T]
    vm_sac = mean(vm_sac, 1, 'omitnan');
    area(x_sac, vm_sac, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
        'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim); ylabel('VM (deg/s)');
    yyaxis left; ylabel('CS rate (Hz)');
    xline(0,'--','HandleVisibility','off'); title('Saccade (CS-on)');

    % ---------- (3) OFFSET ----------
    subplot(3,1,3); hold on; colororder({'k','k'});
    for ii = forced
        [avg, sem] = mean_sem_cells(squeeze(cs_rate_rot(:, t_off, cs_on_dir, ii, E_OFF)));
        yyaxis left
        [hl, hp] = boundedline(x_off, avg, sem, 'Color', colors(ii,:), 'alpha'); %#ok<ASGLU>
        hl.DisplayName = cond_labels{ii}; hp.HandleVisibility = 'off';
        if ii==2 || ii==4, hl.LineStyle='--'; end
    end
    yyaxis right
    vm_off = squeeze(mean(cs_vm_rot(:, t_off, cs_on_dir, forced, E_OFF), 4, 'omitnan')); % [cells x T]
    vm_off = mean(vm_off, 1, 'omitnan');
    area(x_off, vm_off, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
        'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim); ylabel('VM (deg/s)');
    yyaxis left; ylabel('CS rate (Hz)');
    xline(0,'--','HandleVisibility','off'); title('Offset (CS-on)'); xlabel('Time (ms)');
 
    sgtitle(sprintf('CS-on direction — HH / HL / LL / LH   (N = %d cells)', num_cs));
    ESN_Beautify_Plot(fig,[6,9]);
    print(fig, fullfile(img_save_path,'Figure2A-theta.pdf'), '-dpdf', '-bestfit');


    % ---------- (3) OFFSET ----------
    
    avg_all = cell(1,4); sem_all = cell(1,4);
    for ii = forced
        [avg, sem] = mean_sem_cells(squeeze(cs_rate_rot(:, t_off, cs_on_dir, ii, E_OFF)));
        avg_all{ii} = avg; sem_all{ii} = sem;
   
    end
    
    % --- RPE differences at OFFSET ---
    fig = figure; hold on; colororder({'k','k'});
    
    % Cell-by-cell differences
    M_LH = squeeze(cs_rate_rot(:, t_off, cs_on_dir, 4, E_OFF)); % LH
    M_LL = squeeze(cs_rate_rot(:, t_off, cs_on_dir, 3, E_OFF)); % LL
    M_HL = squeeze(cs_rate_rot(:, t_off, cs_on_dir, 2, E_OFF)); % HL
    M_HH = squeeze(cs_rate_rot(:, t_off, cs_on_dir, 1, E_OFF)); % HH
    
    M_rpe_plus  = M_LH - M_LL; % LH - LL
    M_rpe_minus = M_HL - M_HH; % HL - HH
    
    % Compute mean ± SEM across cells
    [avg_plus, sem_plus]   = mean_sem_cells(M_rpe_plus);
    [avg_minus, sem_minus] = mean_sem_cells(M_rpe_minus);
    
    yyaxis left
    [hl1, hp1] = boundedline(x_off, avg_plus, sem_plus, ...
        'cmap',[0 0.6 0], 'alpha'); %#ok<ASGLU>
    hl1.DisplayName = 'RPE+ (LH-LL)'; hp1.HandleVisibility = 'off';
    
    [hl2, hp2] = boundedline(x_off, avg_minus, sem_minus, ...
        'cmap',[0.6 0.3 0], 'alpha'); %#ok<ASGLU>
    hl2.DisplayName = 'RPE- (HL-HH)'; hp2.HandleVisibility = 'off';
    
    yyaxis right
    vm_off = squeeze(mean(cs_vm_rot(:, t_off, cs_on_dir, forced, E_OFF), 4, 'omitnan')); % [cells x T]
    vm_off = mean(vm_off, 1, 'omitnan');
    area(x_off, vm_off, 'FaceColor', vm_color, 'FaceAlpha', .3, ...
        'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim); ylabel('VM (deg/s)');
    
    yyaxis left; ylabel('CS rate (Hz)');
    xline(0,'--','HandleVisibility','off'); 
    title('RPE differences — Primary Saccade Offset (CS-on)'); xlabel('Time (ms)');
    legend('Box','off','Location','northwest');
    
    ESN_Beautify_Plot(fig,[6,4]);
    print(fig, fullfile(img_save_path,'Figure2A-RPE-difference.pdf'), '-dpdf','-bestfit');

end


% -------- helper --------
function [avg, sem] = mean_sem_cells(M)
% M: [cells x time] 
    avg = mean(M, 1, 'omitnan');
    N   = sum(~isnan(M), 1);
    sem = std(M, 0, 1, 'omitnan') ./ max(sqrt(N), 1);
end