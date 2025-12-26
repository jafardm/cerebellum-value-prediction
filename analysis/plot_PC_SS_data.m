function plot_PC_SS_data(data_path)

    img_save_path = fullfile(data_path,'population_figs','ephys');
    save_path     = fullfile(data_path,'population_data');
    if ~exist(img_save_path,'dir'), mkdir(img_save_path); end

    JDM_params_funcs
    S        = load(fullfile(save_path,'PC_SS_dataset_tag1.mat'),'data_out');
    data_ss  = S.data_out.SS;
    data_PC  = S.data_out.PC;
    
    PC_SS_data = cat(1,data_ss.rate_tot_sac,data_PC.rate_tot_sac);
    % concatenate SS and PC for clustering
    PC_SS_data_cat = [mean(data_ss.rate_tot_sac,4,'omitnan');...
        mean(data_PC.rate_tot_sac,4,'omitnan')];   % [cells × 500 × 8]
    num_cells  = size(PC_SS_data_cat,1);

    % --- Burst/Pause clustering ---
    X = reshape(PC_SS_data_cat,num_cells,500*8);          % flatten
    X = X - mean(X,2);                                % mean-center per cell
    red_dim = run_umap(X,'randomize','false','verbose','none');
    clust   = kmeans(red_dim,2,'Replicates',1e3);

    % assign bursters/pausers
    ind_b = clust==2;
    ind_p = clust==1;

    fprintf('Found %d bursters, %d pausers\n',sum(ind_b),sum(ind_p));

   % --- Plot ---
time_ind   = (-100:100)+250;
time_axis  = time_ind-250;   % adjust if your 500 ms window is different
dir_labels = arrayfun(@(x) sprintf('%d°',x),params.sac.ang_values,'UniformOutput',false);

fig = figure('Position',[100 100 1200 600]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');

for d = 1:8
    nexttile; hold on;

    % ===== High reward (conds 1,2) =====
    rate_b_high = squeeze(mean(PC_SS_data(ind_b,time_ind,d,[1 2]),[1 4],'omitnan'));
    sem_b_high  = squeeze(std(mean(PC_SS_data(ind_b,time_ind,d,[1 2]),4,'omitnan'),0,1,'omitnan')) ...
                  ./ sqrt(sum(ind_b));

    rate_p_high = squeeze(mean(PC_SS_data(ind_p,time_ind,d,[1 2]),[1 4],'omitnan'));
    sem_p_high  = squeeze(std(mean(PC_SS_data(ind_p,time_ind,d,[1 2]),4,'omitnan'),0,1,'omitnan')) ...
                  ./ sqrt(sum(ind_p));

    [hl1,hl2] = boundedline(time_axis,rate_b_high,sem_b_high,'r','alpha');
    hl2.HandleVisibility = 'off'; hl1.DisplayName = 'Bursters High';
    [hl1,hl2] = boundedline(time_axis,rate_p_high,sem_p_high,'b','alpha');
    hl2.HandleVisibility = 'off'; hl1.DisplayName = 'Pausers High';

    % ===== Low reward (conds 3,4) =====
    rate_b_low = squeeze(mean(PC_SS_data(ind_b,time_ind,d,[3 4]),[1 4],'omitnan'));
    sem_b_low  = squeeze(std(mean(PC_SS_data(ind_b,time_ind,d,[3 4]),4,'omitnan'),0,1,'omitnan')) ...
                 ./ sqrt(sum(ind_b));

    rate_p_low = squeeze(mean(PC_SS_data(ind_p,time_ind,d,[3 4]),[1 4],'omitnan'));
    sem_p_low  = squeeze(std(mean(PC_SS_data(ind_p,time_ind,d,[3 4]),4,'omitnan'),0,1,'omitnan')) ...
                 ./ sqrt(sum(ind_p));

    [hl1,hl2] = boundedline(time_axis,rate_b_low,sem_b_low,'--r','alpha');
    hl2.HandleVisibility = 'off'; hl1.DisplayName = 'Bursters Low';
    [hl1,hl2] = boundedline(time_axis,rate_p_low,sem_p_low,'--b','alpha');
    hl2.HandleVisibility = 'off'; hl1.DisplayName = 'Pausers Low';

    % formatting
    xline(0,'--k');
    title(['Dir ' dir_labels{d}]);
    xlabel('Time (ms)'); ylabel('Rate (Hz)');

    if d==1
        legend('Location','best','Box','off');
    end
end

ESN_Beautify_Plot(fig,[14 8])
sgtitle('SS + PC population: Bursters vs Pausers (High vs Low reward) across 8 directions');

    

end


