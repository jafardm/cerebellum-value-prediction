function plot_SS_rate_amp_bins(data_path)

    % --- paths ---
    img_save_path = fullfile(data_path,'population_figs','ephys');
    save_path     = fullfile(data_path,'population_data');
    if ~exist(img_save_path,'dir'), mkdir(img_save_path); end

    % --- load ---
    JDM_params_funcs
    S = load(fullfile(save_path,'SS_population_clique_prim_sac.mat'), 'data','meta');
    data_ss = S.data;
    meta    = S.meta;

    % indices of bursters and pausers
    id_cells = load(fullfile(save_path,'purkinje_cell_ids.mat'),'ind_p','ind_b');
    ind_b = id_cells.ind_b;
    ind_p = id_cells.ind_p;
    num_b = sum(ind_b);
    num_p = sum(ind_p);

    % --- indices ---
    time_ind  = (-100:100) + 250;    % aligned to vmax
    time_axis = time_ind - 250;
    nBins     = numel(meta.amp_bin_centers);
    leg_      = {'cs±180','cs±135','cs±90','cs±45','cs-on'};
    colors_   = cool(5);
    vm_color  = [204, 174, 98]/256;
    vm_lim    = [0 600];
    high_inds = [1 2];
    low_inds  = [3 4];
    
    fig = figure;
    sgtitle('SS rate — Bursters vs Pausers × High vs Low reward', 'FontSize', 12)

    cell_groups = {ind_b, ind_p};
    group_names = {'Bursters','Pausers'};

    for g = 1:2 % bursters / pausers
        current_ids = cell_groups{g};
        
        % === High reward ===
        subplot(2,2,(g-1)*2 + 1); hold on;
        colororder({'k','k'});
        title([group_names{g} ' — High reward'])
        xlabel('Time (ms)'); ylabel('Rate (Hz)')
        yyaxis left
        for counter_dir = 5:-1:1
            if counter_dir==1 || counter_dir==5
                d_ = counter_dir;
            else
                d_ = [counter_dir, 10-counter_dir];
            end
            
            tmp = data_ss.rate_tot_sac(current_ids,time_ind,d_,high_inds);
            current_ss_rate = squeeze(mean(tmp,[3 4],'omitnan'));   % → [nCells × nTime]
    
            sig = mean(current_ss_rate,1,'omitnan');
            sig_se = squeeze(std(current_ss_rate,'omitnan'))./sqrt(sum(current_ids));
            
            [hl1,hl2] = boundedline(time_axis,sig,sig_se,'Color',colors_(counter_dir,:),'alpha');
            hl1.DisplayName = leg_{counter_dir};
            hl2.HandleVisibility = 'off';
        end
        xline(0,'--','HandleVisibility','off')
        
        % velocity trace (exclude from legend)
        yyaxis right
        current_vel_sac = squeeze(mean(data_ss.vm_tot_sac(:,time_ind,high_inds),3,'omitnan'));
        sig = mean(current_vel_sac,'omitnan');
        hVel = area(time_axis,sig,'FaceColor',vm_color,'FaceAlpha',.3,'EdgeColor','none');
        hVel.Annotation.LegendInformation.IconDisplayStyle = 'off';
        ylim(vm_lim); ylabel('Velocity (deg/s)')
      
        % === Low reward ===
        subplot(2,2,(g-1)*2 + 2); hold on;
        colororder({'k','k'});
        title([group_names{g} ' — Low reward'])
        xlabel('Time (ms)'); ylabel('Rate (Hz)')
        yyaxis left
        for counter_dir = 5:-1:1
            if counter_dir==1 || counter_dir==5
                d_ = counter_dir;
            else
                d_ = [counter_dir, 10-counter_dir];
            end
            
            tmp = data_ss.rate_tot_sac(current_ids,time_ind,d_,low_inds);
            current_ss_rate = squeeze(mean(tmp,[3 4],'omitnan'));   % → [nCells × nTime]
            sig = squeeze(mean(current_ss_rate,'omitnan'));
            sig_se = squeeze(std(current_ss_rate,'omitnan'))./sqrt(sum(current_ids));
            
            [hl1,hl2] = boundedline(time_axis,sig,sig_se,'Color',colors_(counter_dir,:),'alpha');
            hl1.DisplayName = leg_{counter_dir};
            hl2.HandleVisibility = 'off';
        end
        xline(0,'--','HandleVisibility','off')
        
        % velocity trace (exclude from legend)
        yyaxis right
        current_vel_sac = squeeze(mean(data_ss.vm_tot_sac(:,time_ind,low_inds),3,'omitnan'));
        sig = mean(current_vel_sac,'omitnan');
        hVel = area(time_axis,sig,'FaceColor',vm_color,'FaceAlpha',.3,'EdgeColor','none');
        hVel.Annotation.LegendInformation.IconDisplayStyle = 'off';
        ylim(vm_lim); ylabel('Velocity (deg/s)')
        if g == 1
            legend('Location','northeast','Box','off')
        end
    end
    
    ESN_Beautify_Plot(fig,[12,10])
    % print(fig, fullfile(img_save_path,'SS_rate_burst_pause.pdf'), '-dpdf','-bestfit');

    %% plot amplitude bins
   
    % --- average across directions ---
    % Collapse [cell x time x dir x cond x amp] → [cell x time x cond x amp]
  %% plot amplitude × direction bins

% --- dimensions ---
nDir  = 8; % number of directions
nAmp  = numel(meta.amp_bin_centers);

% Colors: High = red, Low = blue
colors = {'r','b'};

% Layout: amplitude bins × 2 groups (Bursters / Pausers), inside each panel plot directions
fig = figure;
t = tiledlayout(nAmp,2,'TileSpacing','compact','Padding','compact');
title(t,'SS rate by direction × amplitude × reward')

cell_groups = {ind_b, ind_p};
group_names = {'Bursters','Pausers'};

for b = 1:nAmp
    bin_lo = meta.amp_edges(b);
    bin_hi = meta.amp_edges(b+1);

    for g = 1:2
        current_ids = cell_groups{g};
        ax = nexttile((b-1)*2 + g); hold(ax,'on')

        % loop over directions
         for counter_dir = 5:-1:1
            if counter_dir==1 || counter_dir==5
                d = counter_dir;
            else
                d = [counter_dir, 10-counter_dir];
            end
            % --- High reward ---
            dat = squeeze(data_ss.rate_tot_sac_amp(current_ids,time_ind,d,[1 2],b));
            dat = mean(dat,3,'omitnan'); % collapse cond1+2
            m_rate = mean(dat,'omitnan');
            s_rate = std(dat,[],1,'omitnan') / sqrt(sum(current_ids));
            [hl1,hl2]= boundedline(time_axis, m_rate, s_rate, colors{1}, 'alpha');
            hl1.DisplayName = leg_{counter_dir};
            hl2.HandleVisibility = 'off';

            % --- Low reward ---
            dat = squeeze(data_ss.rate_tot_sac_amp(current_ids,time_ind,d,[3 4],b));
            dat = mean(dat,3,'omitnan'); % collapse cond3+4
            m_rate = mean(dat,'omitnan');
            s_rate = std(dat,[],1,'omitnan') / sqrt(sum(current_ids));
            boundedline(time_axis, m_rate, s_rate, colors{2}, 'alpha');
        end

        xline(0,'--k','HandleVisibility','off')
        xlabel('Time from decel (ms)')
        ylabel('Rate (Hz)')
        title(sprintf('%s — Amp %d–%d°',group_names{g}, bin_lo, bin_hi))

        if b == 1 && g == 1
            legend({'High','Low'},'Box','off','Location','best')
        end
    end
end

ESN_Beautify_Plot(fig,[14,10])

    % print(fig, fullfile(img_save_path,'SS_rate_amp_bins_burst_pause.pdf'), '-dpdf','-bestfit');

    

end
