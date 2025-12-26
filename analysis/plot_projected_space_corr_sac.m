function plot_projected_space_corr_sac(data_path)

 img_save_path = fullfile(data_path,'population_figs','ephys');
 save_path     = fullfile(data_path,'population_data');
 if ~exist(img_save_path,'dir'); mkdir(img_save_path); end
 
 JDM_params_funcs
 data_ss  = load(fullfile(save_path,'SS_population_clique_corr_sac.mat'), 'data');
 data_ml1 = load(fullfile(save_path,'MLI_population_clique_corr_sac.mat'), 'data');
 load(fullfile(save_path,'purkinje_cell_ids.mat'),'ind_p','ind_b');

 data_ss  = data_ss.data;
 data_ml1 = data_ml1.data;

 sess_SS  = cell2mat(cellfun(@(x) str2double(x(1:6)), data_ss.cell_ids_tot,'UniformOutput',false));
 sess_MLI = cell2mat(cellfun(@(x) str2double(x(1:6)), data_ml1.cell_ids_tot,'UniformOutput',false));

 cliques_SS   = sess_SS.*10+data_ss.cliques;
 cliques_MLI  = sess_MLI.*10+data_ml1.cliques;

 num_SS   = numel(cliques_SS);
 num_MLI  = numel(cliques_MLI);
 num_b    = sum(ind_b);
 num_p    = sum(ind_p);

 ss_cs_rho   = data_ss.cs_on_rho_tot(:);
 ml1_cs_rho  = data_ml1.cs_on_rho_tot(:);

 %% plot settings
 vm_color = [204, 174, 98]/256;
 vm_lim   = [0,400];
 time_ind = (-100:150)+250;  % decel-aligned
 proj_lim = [-30,50];

 % Four explicit conditions in the 4th dim (exact order)
 cond_labels = {'HH','HL','LL','LH'};  % MUST match data 4th-dim order
 n_conds = numel(cond_labels);

 % Slice per-condition data (no averaging)
 data_ss_cond  = cell(1,n_conds);
 data_vm_cond  = cell(1,n_conds);
 data_ml1_cond = cell(1,n_conds);

 for ic = 1:n_conds
     data_ss_cond{ic}  = data_ss.rate_tot_sac(:,:,:, ic);
     data_vm_cond{ic}  = data_ss.vm_tot_sac(:,:,:,  ic);
     data_ml1_cond{ic} = data_ml1.rate_tot_sac(:,:,:,ic); 
 end

 %% ---------- Figure 1: 4x2 (each row = condition; cols = p_y(t), p_x(t)) ----------
 fig = figure('Color','w','WindowStyle','normal');
 t = tiledlayout(n_conds,2,'TileSpacing','compact','Padding','compact');
 title(t,'SS projections corrective Sac. (CS-on + 90° and cs-on°)')
 colororder({'k','k'})
 for ic = 1:n_conds
     current_data = data_ss_cond{ic};
     current_vm   = data_vm_cond{ic};
     cond_label   = cond_labels{ic};

     % =========================
     % (1) p_y(t): projection onto sac dir + 90°
     % =========================
     ax = nexttile((ic-1)*2 + 1);  hold(ax,'on')
     proj_w = -cosd((0:45:360-45)+90);

     % --- Bursters ---
     rate_pr_b = zeros(num_b,numel(time_ind));
     for d = 1:8
         rate_pr_b = rate_pr_b + current_data(ind_b,time_ind,d).*proj_w(d);
     end
     rate_pr_b = rate_pr_b .* (ss_cs_rho(ind_b) / sum(ss_cs_rho(ind_b))) * num_b;  % rho-weighted rows

     yyaxis left
     colororder({'k','k'})
     [hl,hp] = boundedline(time_ind-250, mean(rate_pr_b,'omitnan'), std(rate_pr_b,'omitnan')/sqrt(max(1,num_b)), 'r','alpha');
     hl.DisplayName = sprintf('SS burst (n = %d)', num_b);
     hp.HandleVisibility = 'off';
     ylabel('Firing rate (Hz)')

     % --- Pausers ---
     rate_pr_p = zeros(num_p,numel(time_ind));
     for d = 1:8
         rate_pr_p = rate_pr_p + current_data(ind_p,time_ind,d).*proj_w(d);
     end
     rate_pr_p = rate_pr_p .* (ss_cs_rho(ind_p) / sum(ss_cs_rho(ind_p))) * num_p;

     [hl,hp] = boundedline(time_ind-250, mean(rate_pr_p,'omitnan'), std(rate_pr_p,'omitnan')/sqrt(max(1,num_p)), 'b','alpha');
     hl.DisplayName = sprintf('SS pause (n = %d)', num_p);
     hp.HandleVisibility = 'off';

     % --- All SS ---
     rate_pr_ss = zeros(num_SS,numel(time_ind));
     for d = 1:8
         rate_pr_ss = rate_pr_ss + current_data(:,time_ind,d).*proj_w(d);
     end
     rate_pr_ss = rate_pr_ss .* (ss_cs_rho / sum(ss_cs_rho)) * num_SS;

     [hl,hp] = boundedline(time_ind-250, mean(rate_pr_ss,'omitnan'), std(rate_pr_ss,'omitnan')/sqrt(max(1,num_SS)), 'k','alpha');
     hl.DisplayName = sprintf('SS all (n = %d)', num_SS);
     hp.HandleVisibility = 'off';

     xline(0,'--','HandleVisibility','off'); yline(0,'--','HandleVisibility','off')
     xlim([-100,100]); ylim(proj_lim)
     xlabel('Time from deceleration (ms)')

     yyaxis right
     current_vel_sac = mean(mean(current_vm(:,time_ind),3,'omitnan'),1,'omitnan');
     area(time_ind-250,current_vel_sac,'FaceColor',vm_color,'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off')
     ylim(vm_lim); 

     if ic == 1, legend('Box','off','Location','northwest'); end
     title([cond_label ' — p_y(t)'])

     % =========================
     % (2) p_x(t): projection onto sac dir (0°)
     % =========================
     ax = nexttile((ic-1)*2 + 2);  hold(ax,'on')
     proj_w = -cosd((0:45:360-45));

     % --- Bursters ---
     rate_pr_b = zeros(num_b,numel(time_ind));
     for d = 1:8
         rate_pr_b = rate_pr_b + current_data(ind_b,time_ind,d).*proj_w(d);
     end
     rate_pr_b = rate_pr_b .* (ss_cs_rho(ind_b) / sum(ss_cs_rho(ind_b))) * num_b;

     yyaxis left
     colororder({'k','k'})
     [hl,hp] = boundedline(time_ind-250, mean(rate_pr_b,'omitnan'), std(rate_pr_b,'omitnan')/sqrt(max(1,num_b)), 'r','alpha');
     hl.DisplayName = 'SS burst'; hp.HandleVisibility = 'off';

     % --- Pausers ---
     rate_pr_p = zeros(num_p,numel(time_ind));
     for d = 1:8
         rate_pr_p = rate_pr_p + current_data(ind_p,time_ind,d).*proj_w(d);
     end
     rate_pr_p = rate_pr_p .* (ss_cs_rho(ind_p) / sum(ss_cs_rho(ind_p))) * num_p;

     [hl,hp] = boundedline(time_ind-250, mean(rate_pr_p,'omitnan'), std(rate_pr_p,'omitnan')/sqrt(max(1,num_p)), 'b','alpha');
     hl.DisplayName = 'SS pause'; hp.HandleVisibility = 'off';

     % --- All SS ---
     rate_pr_ss = zeros(num_SS,numel(time_ind));
     for d = 1:8
         rate_pr_ss = rate_pr_ss + current_data(:,time_ind,d).*proj_w(d);
     end
     rate_pr_ss = rate_pr_ss .* (ss_cs_rho / sum(ss_cs_rho)) * num_SS;

     [hl,hp] = boundedline(time_ind-250, mean(rate_pr_ss,'omitnan'), std(rate_pr_ss,'omitnan')/sqrt(max(1,num_SS)), 'k','alpha');
     hl.DisplayName = 'SS all'; hp.HandleVisibility = 'off';

     xline(0,'--','HandleVisibility','off'); yline(0,'--','HandleVisibility','off')
     xlim([-100,100]); ylim(proj_lim)
     xlabel('Time from deceleration (ms)')

     yyaxis right
     current_vel_sac = mean(mean(current_vm(:,time_ind),3,'omitnan'),1,'omitnan');
     area(time_ind-250,current_vel_sac,'FaceColor',vm_color,'FaceAlpha',.3,'EdgeColor','none','HandleVisibility','off')
     ylim(vm_lim); ylabel('Velocity (deg/s)')
     title([cond_label ' — p_x(t)'])
 end

 ESN_Beautify_Plot(fig,[8,12])
 print(fig, fullfile(img_save_path,'ss-rate-purkinje-Projection_corr-sac.pdf'), '-dpdf', '-bestfit');

 %% ---------- Figure 2: 4x3 SS & MLI1 (and optional MLI2) ----------
 fig = figure;
 t = tiledlayout(n_conds,3,'TileSpacing','compact','Padding','compact');
 title(t,'SS & MLI1 Projection per condition corrective Sacs.(CS-on+90° and cs-on°)')
 colororder({'k','k'})   
 for ic = 1:n_conds
    cond_label   = cond_labels{ic};
    current_data = data_ss_cond{ic};      % [Nss x T x 8]
    current_vm   = data_vm_cond{ic};      % [Nss x T x 8]
    current_ml1  = data_ml1_cond{ic};     % [Nml1 x T x 8]

    % --- shared helpers: velocity trace & safe denominators for rho-weights
    vel_trace = mean(mean(current_vm(:,time_ind,:),3,'omitnan'),1);      % [1 x T]
    vel_trace = squeeze(vel_trace);                                      % [T x 1]
    denom_b   = sum(ss_cs_rho(ind_b), 'omitnan'); if denom_b==0, denom_b = num_b; end
    denom_p   = sum(ss_cs_rho(ind_p), 'omitnan'); if denom_p==0, denom_p = num_p; end
    denom_ss  = sum(ss_cs_rho,        'omitnan'); if denom_ss==0, denom_ss = num_SS; end
    denom_ml1 = sum(ml1_cs_rho,       'omitnan'); if denom_ml1==0, denom_ml1 = num_MLI; end

    % =========================
    % (1) p_y(t): +90°
    % =========================
    ax = nexttile((ic-1)*3 + 1); hold(ax,'on')
    proj_w = -cosd((0:45:360-45)+90);

    % --- SS groups ---
    rate_pr_b  = zeros(num_b, numel(time_ind));
    rate_pr_p  = zeros(num_p, numel(time_ind));
    rate_pr_ss = zeros(num_SS,numel(time_ind));
    for d=1:8
        rate_pr_b  = rate_pr_b  + current_data(ind_b,time_ind,d).*proj_w(d);
        rate_pr_p  = rate_pr_p  + current_data(ind_p,time_ind,d).*proj_w(d);
        rate_pr_ss = rate_pr_ss + current_data(:,     time_ind,d).*proj_w(d);
    end
    rate_pr_b  = rate_pr_b  .* (ss_cs_rho(ind_b)/denom_b) * num_b;
    rate_pr_p  = rate_pr_p  .* (ss_cs_rho(ind_p)/denom_p) * num_p;
    rate_pr_ss = rate_pr_ss .* (ss_cs_rho       /denom_ss) * num_SS;

    yyaxis left
    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_b,'omitnan'),  std(rate_pr_b,'omitnan')/sqrt(max(1,num_b)), 'r','alpha'); hl.DisplayName = sprintf('SS burst (n=%d)',num_b); hp.HandleVisibility='off';
    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_p,'omitnan'),  std(rate_pr_p,'omitnan')/sqrt(max(1,num_p)), 'b','alpha'); hl.DisplayName = sprintf('SS pause (n=%d)',num_p); hp.HandleVisibility='off';
    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_ss,'omitnan'), std(rate_pr_ss,'omitnan')/sqrt(max(1,num_SS)),'k','alpha'); hl.DisplayName = sprintf('SS all (n=%d)',num_SS); hp.HandleVisibility='off';
    ylabel('Firing rate (Hz)')

    % --- MLI1 (+90°) ---
    rate_pr_ml1_y = zeros(num_MLI,numel(time_ind));
    for d=1:8
        rate_pr_ml1_y = rate_pr_ml1_y + current_ml1(:,time_ind,d).*proj_w(d);
    end
    rate_pr_ml1_y = rate_pr_ml1_y .* (ml1_cs_rho/denom_ml1) * num_MLI;
    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_ml1_y,'omitnan'), std(rate_pr_ml1_y,'omitnan')/sqrt(max(1,num_MLI)),'g','alpha'); hl.DisplayName = sprintf('MLI1 (n=%d)',num_MLI); hp.HandleVisibility='off';

    % --- Velocity overlay ---
    yyaxis right
    colororder({'k','k'})
    area(time_ind-250, vel_trace, 'FaceColor', vm_color, 'FaceAlpha', .3, 'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim); 
    
    yyaxis left
    xline(0,'--','HandleVisibility','off'); yline(0,'--','HandleVisibility','off')
    xlim([-100,100]); ylim(ax, proj_lim)  % left axis limits preserved
    xlabel('Time from deceleration (ms)')
    if ic == 1, legend('Box','off','Location','northwest'); end
    title([cond_label ' — p_y(t)'])

    % =========================
    % (2) p_x(t): 0°
    % =========================
    ax = nexttile((ic-1)*3 + 2); hold(ax,'on')
    proj_w = -cosd((0:45:360-45));

    % --- SS groups (0°) ---
    rate_pr_b  = zeros(num_b, numel(time_ind));
    rate_pr_p  = zeros(num_p, numel(time_ind));
    rate_pr_ss = zeros(num_SS,numel(time_ind));
    for d=1:8
        rate_pr_b  = rate_pr_b  + current_data(ind_b,time_ind,d).*proj_w(d);
        rate_pr_p  = rate_pr_p  + current_data(ind_p,time_ind,d).*proj_w(d);
        rate_pr_ss = rate_pr_ss + current_data(:,     time_ind,d).*proj_w(d);
    end
    rate_pr_b  = rate_pr_b  .* (ss_cs_rho(ind_b)/denom_b) * num_b;
    rate_pr_p  = rate_pr_p  .* (ss_cs_rho(ind_p)/denom_p) * num_p;
    rate_pr_ss = rate_pr_ss .* (ss_cs_rho       /denom_ss) * num_SS;

    yyaxis left
    colororder({'k','k'})
    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_b,'omitnan'),  std(rate_pr_b,'omitnan')/sqrt(max(1,num_b)), 'r','alpha'); hl.DisplayName = 'SS burst'; hp.HandleVisibility='off';
    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_p,'omitnan'),  std(rate_pr_p,'omitnan')/sqrt(max(1,num_p)), 'b','alpha'); hl.DisplayName = 'SS pause'; hp.HandleVisibility='off';
    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_ss,'omitnan'), std(rate_pr_ss,'omitnan')/sqrt(max(1,num_SS)),'k','alpha'); hl.DisplayName = 'SS all'; hp.HandleVisibility='off';
   

    % --- NEW: MLI1 (0°) ---
    rate_pr_ml1_x = zeros(num_MLI,numel(time_ind));
    for d=1:8
        rate_pr_ml1_x = rate_pr_ml1_x + current_ml1(:,time_ind,d).*proj_w(d);
    end
    rate_pr_ml1_x = rate_pr_ml1_x .* (ml1_cs_rho/denom_ml1) * num_MLI;
    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_ml1_x,'omitnan'), std(rate_pr_ml1_x,'omitnan')/sqrt(max(1,num_MLI)),'g','alpha'); hl.DisplayName = 'MLI1'; hp.HandleVisibility='off';

    % --- Velocity overlay ---
    yyaxis right
    area(time_ind-250, vel_trace, 'FaceColor', vm_color, 'FaceAlpha', .3, 'EdgeColor','none','HandleVisibility','off');
    ylim(vm_lim); ylabel('Velocity (deg/s)')
    yyaxis left
    xline(0,'--','HandleVisibility','off'); yline(0,'--','HandleVisibility','off')
    xlim([-100,100]); ylim(ax, proj_lim)
    xlabel('Time from deceleration (ms)')
    title([cond_label ' — p_x(t)'])

    % =========================
    % (3) Phase plot: use the 0° (p_x) signals
    % =========================
    ax = nexttile((ic-1)*3 + 3); hold(ax,'on')
    sig_ml1 = mean(rate_pr_ml1_x,'omitnan');  % 0° MLI
    sig_ss  = mean(rate_pr_ss,'omitnan');
    sig_b   = mean(rate_pr_b,'omitnan');
    sig_p   = mean(rate_pr_p,'omitnan');

    plot(sig_ml1, sig_ss, 'k', 'DisplayName', 'SS all');
    plot(sig_ml1, sig_b,  'r', 'DisplayName', 'SS burst');
    plot(sig_ml1, sig_p,  'b', 'DisplayName', 'SS pause');
    if numel(sig_ml1) >= 101
        scatter(sig_ml1(101), sig_ss(101), 8, 'k', 'filled');
        scatter(sig_ml1(101), sig_b(101),  8, 'r', 'filled');
        scatter(sig_ml1(101), sig_p(101),  8, 'b', 'filled');
    end

    axis equal
    xlim(proj_lim); ylim(proj_lim)
    xline(0); yline(0)
    xlabel('MLI1 (Hz)'); ylabel('PC (Hz)')
    title([cond_label ' — PC vs MLI1'])
end

 ESN_Beautify_Plot(fig,[10,12])
 print(fig, fullfile(img_save_path,'ss-rate-population-Projection-corr-sac.pdf'), '-dpdf', '-bestfit');

end
