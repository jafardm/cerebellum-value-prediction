function plot_projected_space_prim_sac4conds(data_path)

img_save_path = fullfile(data_path,'population_figs','ephys');
save_path     = fullfile(data_path,'population_data');

JDM_params_funcs
data_ss  = load(fullfile(save_path,'SS_population_clique_prim_sac.mat'), 'data');
load(fullfile(save_path,'purkinje_cell_ids.mat'),'ind_p','ind_b');

data_ss  = data_ss.data;

num_b = sum(ind_b);
num_p = sum(ind_p);
num_SS = size(data_ss.rate_tot_sac,1);
ss_cs_rho = data_ss.cs_on_rho_tot;

% ---- Params ----
vm_color = [204, 174, 98]/256;
vm_lim   = [0,600];
time_ind = (-100:150)+250;
proj_lim = [-30,55];

cond_labels = {'HH','HL','LL','LH'};

% === One figure with 4x2 tiles (each row = condition; 2 cols = py, px) ===
fig = figure;
t = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
title(t,'SS projections - Primary Saccade (cs-on+90° and cs-on°)')

for i_cond = 1:4
    current_data = data_ss.rate_tot_sac(:,:,:, i_cond);
    current_vm   = data_ss.vm_tot_sac(:,:,:, i_cond);
    cond_label   = cond_labels{i_cond};
    % =========================
    % (1) Projection onto sac dir + 90  -> p_y(t)
    % =========================
    ax = nexttile((i_cond-1)*2 + 1);  
    hold(ax,'on')
    proj_w = -cosd((0:45:360-45)+90);

   % Bursters
    rate_pr_b = zeros(num_b,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_b = rate_pr_b + current_data(ind_b,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_b = rate_pr_b * num_b .* ss_cs_rho(ind_b) / sum(ss_cs_rho(ind_b));

    yyaxis left
    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_b,'omitnan'), ...
        std(rate_pr_b,0,'omitnan')/sqrt(num_b),'r','alpha');
    hl.DisplayName = 'Bursters'; hp.HandleVisibility = 'off';

    % Pausers
    rate_pr_p = zeros(num_p,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_p = rate_pr_p + current_data(ind_p,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_p = rate_pr_p * num_p .* ss_cs_rho(ind_p) / sum(ss_cs_rho(ind_p));

    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_p,'omitnan'), ...
        std(rate_pr_p,0,'omitnan')/sqrt(num_p),'b','alpha');
    hl.DisplayName = 'Pausers'; hp.HandleVisibility = 'off';

    % --- All SS ---
   % All SS
    rate_pr_ss = zeros(num_SS,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_ss = rate_pr_ss + current_data(:,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_ss = rate_pr_ss * num_SS .* ss_cs_rho / sum(ss_cs_rho);


    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_ss,'omitnan'), ...
        std(rate_pr_ss,0,'omitnan')/sqrt(num_SS),'k','alpha');
    hl.DisplayName = 'All SS'; hp.HandleVisibility = 'off';

    xline(0,'--'); yline(0,'--')
    xlim([-100,100]); ylim(proj_lim)
    xlabel('Time (ms)'); ylabel('Rate (Hz)')
    yyaxis right
    current_vel = mean(mean(current_vm(:,time_ind),3,'omitnan'),1);
    current_vel(isnan(current_vel)) = 0;
    area(time_ind-250,current_vel,'FaceColor',vm_color,'FaceAlpha',.3,'EdgeColor','none')
    ylim(vm_lim)
    title([cond_label ' — p_y(t)'])
    if i_cond==1, legend('Location','northwest','Box','off'); end

    % =========================
    % (2) Projection onto sac dir -> p_x(t)
    % =========================
    ax = nexttile((i_cond-1)*2 + 2);  
    hold(ax,'on')
    proj_w = -cosd((0:45:360-45));

       % Bursters
    rate_pr_b = zeros(num_b,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_b = rate_pr_b + current_data(ind_b,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_b = rate_pr_b * num_b .* ss_cs_rho(ind_b) / sum(ss_cs_rho(ind_b));

    yyaxis left
    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_b,'omitnan'), ...
        std(rate_pr_b,0,'omitnan')/sqrt(num_b),'r','alpha');
    hl.DisplayName = 'Bursters'; hp.HandleVisibility = 'off';

    % Pausers
    rate_pr_p = zeros(num_p,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_p = rate_pr_p + current_data(ind_p,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_p = rate_pr_p * num_p .* ss_cs_rho(ind_p) / sum(ss_cs_rho(ind_p));

    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_p,'omitnan'), ...
        std(rate_pr_p,0,'omitnan')/sqrt(num_p),'b','alpha');
    hl.DisplayName = 'Pausers'; hp.HandleVisibility = 'off';

    % --- All SS ---
    rate_pr_ss = zeros(num_SS,numel(time_ind));
    for counter_dir = 1:8
        rate_pr_ss = rate_pr_ss + current_data(:,time_ind,counter_dir).*proj_w(counter_dir);
    end
    rate_pr_ss = rate_pr_ss * num_SS .* ss_cs_rho / sum(ss_cs_rho);


    [hl,hp] = boundedline(time_ind-250, mean(rate_pr_ss,'omitnan'), ...
        std(rate_pr_ss,0,'omitnan')/sqrt(num_SS),'k','alpha');
    hl.DisplayName = 'All SS'; hp.HandleVisibility = 'off';
    xline(0,'--'); yline(0,'--')
    xlim([-100,100]); ylim(proj_lim)
    xlabel('Time (ms)')
    yyaxis right
    current_vel = mean(mean(current_vm(:,time_ind),3,'omitnan'),1);
    current_vel(isnan(current_vel)) = 0;
    area(time_ind-250,current_vel,'FaceColor',vm_color,'FaceAlpha',.3,'EdgeColor','none')
    ylim(vm_lim); ylabel('Velocity (deg/s)')
    title([cond_label ' — p_x(t)'])
end
ESN_Beautify_Plot(fig,[8,10])

%print(fig, fullfile(img_save_path,'ss-rate-projection-prim-sac-4conds.pdf'),'-dpdf','-bestfit');

end
