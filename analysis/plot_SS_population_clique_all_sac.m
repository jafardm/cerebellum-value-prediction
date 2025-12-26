function plot_SS_population_clique_all_sac(data_path)
% SS, MLI rates for tags 1,4,6,7 saccdes
    img_save_path = fullfile(data_path,'population_figs\ephys');
    save_path = fullfile(data_path,'population_data');
     
    JDM_params_funcs
    data_ss  = load(fullfile(save_path,'SS_population_clique_all_sac.mat'), 'data');
    data_ml1 = load(fullfile(save_path,'MLI_population_clique_all_sac.mat'), 'data');
    data_ml2 = load(fullfile(save_path,'MLI2_population_clique_all_sac.mat'), 'data');
    %data_cs = load(fullfile(save_path,'CS_population_clique_prim_sac'), 'data');
    
    data_ss  = data_ss.data;
    data_ml1 = data_ml1.data;
    data_ml2 = data_ml2.data;
    
    sess_SS      = cell2mat(cellfun(@(x) str2double(x(1:6)),...
        data_ss.cell_ids_tot,'UniformOutput',false));
    sess_MLI     = cell2mat(cellfun(@(x) str2double(x(1:6)),...
        data_ml1.cell_ids_tot,'UniformOutput',false));
    sess_MLI2    = cell2mat(cellfun(@(x) str2double(x(1:6)),...
        data_ml2.cell_ids_tot,'UniformOutput',false));
    
    cliques_SS   = sess_SS.*10+data_ss.cliques;
    cliques_MLI  = sess_MLI.*10+data_ml1.cliques;
    cliques_MLI2 = sess_MLI2.*10+data_ml2.cliques;
    
    num_SS   = numel(cliques_SS);
    num_MLI  = numel(cliques_MLI);
    num_MLI2 = numel(cliques_MLI2);
    
    %burst-pause clustering using umap+kmeans: uncomment to rerun
    %(results might change slightly due to the randomness of umap and clustering)
%     ss_rate_tot_sac_cat_t = reshape(data_ss.rate_tot_sac,num_SS,500*8);
%     
%     red_dim_ss = run_umap(ss_rate_tot_sac_cat_t-mean(ss_rate_tot_sac_cat_t,2),...
%         'randomize','false','verbose','none');
%     clust = kmeans(red_dim_ss,2,'Replicates',1e3);
%     ind_b = clust == 2;
%     ind_p = ~ind_b;
% 
%     save(fullfile(data_path,'population_data', 'purkinje_cell_ids.mat'),'ind_p','ind_b');
    load(fullfile(save_path,'purkinje_cell_ids.mat'),'ind_p','ind_b')
    num_b = sum(ind_b);
    num_p = sum(ind_p);
    
    fprintf('\n================  number of the bursters %d \n',num_b)
    fprintf('\n================  number of the pausers %d \n',num_p)

    n_col = 4;
    n_row = 1;
    
    vm_color = [204, 174, 98]/256;
    vm_lim = [0,600];
    
    time_ind = (-100:100)+250;
    rate_lim = [35,90];
 
    leg_ = {'cs±180','cs±135','cs±90','cs±45','cs-on'};
    colors_ = cool(5);
    colors_ = colors_(end:-1:1,:);
  
    fig = figure;
    sgtitle('SS, MLI1, MLI2 population rates — Tags 1,4,6,7', 'FontSize', 12)
   
    for counter_dir = 5:-1:1
    
        current_dir = counter_dir;
    
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            current_dir_ = [current_dir, 10-current_dir];
        end
        num_b = sum(ind_b);
        current_ss_rate_sac = squeeze(mean(data_ss.rate_tot_sac(ind_b,time_ind,current_dir_),3));
        sig = squeeze(mean(current_ss_rate_sac));
        sig_se = squeeze(std(current_ss_rate_sac))./sqrt(num_b);

        subplot(n_row,n_col,1)
        colororder({'k','k'})
        yyaxis left
        [hl1,hl2] = boundedline(time_ind-250,sig,sig_se,'Color',colors_(counter_dir,:),'alpha');
        hl2.HandleVisibility = 'off';
        hl1.DisplayName = '';
        hold on;
        yyaxis left
        ylabel('Rate (Hz)')
        xline(0,'--','HandleVisibility','off');
    
        if counter_dir == 1
            yyaxis right
            colororder({'k','k'})
            hold on;
            current_vel_sac = squeeze(mean(data_ss.vm_tot_sac(:,time_ind,:),3));
            sig = squeeze(mean(current_vel_sac));
            area(time_ind-250,sig,'FaceColor',vm_color,...
                'FaceAlpha',.3,'EdgeColor','none');
            ylim(vm_lim);
        end
    end
    
    for counter_dir = 5:-1:1
    
        current_dir = counter_dir;
    
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            current_dir_ = [current_dir, 10-current_dir];
        end
        num_p = sum(ind_p);
        current_ss_rate_sac = squeeze(mean(data_ss.rate_tot_sac(ind_p,time_ind,current_dir_),3));
        sig = squeeze(mean(current_ss_rate_sac));
        sig_se = squeeze(std(current_ss_rate_sac))./sqrt(num_p);
    
        yyaxis left
        [hl1,hl2] = boundedline(time_ind-250,sig,sig_se,'Color',colors_(counter_dir,:),'alpha');
        hl2.HandleVisibility = 'off';
        hl1.DisplayName = leg_{counter_dir};
        hold on;
        xline(0,'--','HandleVisibility','off');
        ylim(rate_lim);
        xlabel('Time max vel.(ms)')
    end
    title('Pcells: burster/pauser')
    
    % P-Cells
    subplot(n_row,n_col,2)
    colororder({'k','k'})
    colors_ = cool(5);
    colors_ = colors_(end:-1:1,:);
    
    for counter_dir = 5:-1:1
        current_dir = counter_dir;
    
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            current_dir_ = [current_dir, 10-current_dir];
        end
        current_ss_rate_sac = squeeze(mean(data_ss.rate_tot_sac(:,time_ind,current_dir_),3));
        num_SS = size(current_ss_rate_sac,1);
        sig = squeeze(mean(current_ss_rate_sac));
        sig_se = std(current_ss_rate_sac)./sqrt(num_SS);
    
        yyaxis left
        [hl1,hl2] = boundedline(time_ind-250,sig,sig_se,'Color',colors_(counter_dir,:),'alpha');
        hl2.HandleVisibility = 'off';
        hl1.DisplayName = leg_{counter_dir};
        hold on;
        xline(0,'--','HandleVisibility','off');
        ylim(rate_lim)
        xlabel('Time max vel.(ms)')
        % plot velocities
        if counter_dir == 1
            yyaxis right
            hold on;
            current_vel_sac = squeeze(mean(data_ss.vm_tot_sac(:,time_ind,:),3));
            sig = squeeze(mean(current_vel_sac));
            area(time_ind-250,sig,'FaceColor',vm_color,...
                'FaceAlpha',.3,'EdgeColor','none')
            ylim(vm_lim);
        end
    end
    title('Pcells')
    
    % MLI1
    subplot(n_row,n_col,3)
    colororder({'k','k'})
    colors_ = cool(5);
    colors_ = colors_(end:-1:1,:);
    for counter_dir = 5:-1:1
    
        current_dir = counter_dir;
    
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            current_dir_ = [current_dir, 10-current_dir];
        end
        num_ml1 = size(data_ml1.rate_tot_sac,1);
        current_mli_rate_sac = squeeze(mean(data_ml1.rate_tot_sac(:,time_ind,current_dir_),3));
        sig = squeeze(mean(current_mli_rate_sac));
        sig_se = squeeze(std(current_mli_rate_sac))./sqrt(num_ml1);
    
        yyaxis left
        [hl1,hl2] = boundedline(time_ind-250,sig,sig_se,'Color',colors_(counter_dir,:),'alpha');
        hl2.HandleVisibility = 'off';
        hl1.DisplayName = leg_{counter_dir};
        hold on;
        xline(0,'--','HandleVisibility','off');
        ylim([20 90])
        xlabel('Time max vel.(ms)')
        if counter_dir == 1
            yyaxis right
            hold on;
            current_vel_sac = squeeze(mean(data_ml1.vm_tot_sac(:,time_ind,:),3));
            sig = squeeze(mean(current_vel_sac));
            area(time_ind-250,sig,'FaceColor',vm_color,...
                'FaceAlpha',.3,'EdgeColor','none')
            ylim(vm_lim);
        end
    end
    title('MLI1')
    
    % MLI2
    subplot(n_row,n_col,4)
    colororder({'k','k'})
    colors_ = cool(5);
    colors_ = colors_(end:-1:1,:);
    for counter_dir = 5:-1:1
    
        current_dir = counter_dir;
    
        if current_dir == 1 || current_dir == 5
            current_dir_ = current_dir;
        else
            current_dir_ = [current_dir, 10-current_dir];
        end
        num_ml2 = size(data_ml2.rate_tot_sac,1);
        current_mli2_rate_sac = squeeze(mean(data_ml2.rate_tot_sac(:,time_ind,current_dir_),3));
        sig = squeeze(mean(current_mli2_rate_sac));
        sig_se = squeeze(std(current_mli2_rate_sac))./sqrt(num_ml2);
    
        yyaxis left
        [hl1,hl2] = boundedline(time_ind-250,sig,sig_se,'Color',colors_(counter_dir,:),'alpha');
        hl2.HandleVisibility = 'off';
        hl1.DisplayName = leg_{counter_dir};
        hold on;
        ylim([0 50])
        xline(0,'--','HandleVisibility','off');
        xlabel('Time max vel.(ms)')
        if counter_dir == 1
            yyaxis right
            ylabel('Velocity (deg/sec)')
            hold on;
            current_vel_sac = squeeze(mean(data_ml2.vm_tot_sac(:,time_ind,:),3));
            sig = squeeze(mean(current_vel_sac));
            area(time_ind-250,sig,'FaceColor',vm_color,...
                'FaceAlpha',.3,'EdgeColor','none')
            ylim(vm_lim);
        end
    end
    title('MLI2')
    legend(leg_, 'Location','northwest', 'Box','off');
    ESN_Beautify_Plot(fig,[10 4])
    print(fig, fullfile(img_save_path,'ss-rate-all-sacs.pdf'), '-dpdf', '-bestfit');

end