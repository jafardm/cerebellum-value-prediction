function plot_cs_spe_rpe_catch_up_sacs(data_path)
 JDM_params_funcs
 S = load(fullfile(data_path,'population_data','cs_on_rate_catch_up_sac_sac'));
 img_save_path = fullfile(data_path,'population_figs\ephys');  
 num_cs = numel(S.cs_on_data);  
% bulid pop[ulation data 
 CS_poulation_sac = S.cs_rate; %  2 saccade
 VM_population_sac = S.cs_vm;

 % --- extract amplitudes and counts ---
amp_all = S.N_amp_all;          % {animal}{session} 
% concatenate across animals/sessions
all_amps = cellfun(@(x) vertcat(x{:}), amp_all, 'UniformOutput', false);
all_amps = vertcat(all_amps{:});

N_lat_all = S.N_lat_all;
N_lat_cond_all = S.N_lat_cond_all;
% Concatenate across animals/sessions
all_lat     = cellfun(@(x) vertcat(x{:}), N_lat_all, 'UniformOutput', false);
all_lat     = vertcat(all_lat{:});
all_lat_cond = cellfun(@(x) vertcat(x{:}), N_lat_cond_all, 'UniformOutput', false);
all_lat_cond = vertcat(all_lat_cond{:});

t_ind_sac = (-100:100) + 250;
colors = [190, 37, 43; 255,  0,  230; 0, 10, 255;0, 204, 255]/ 255;
vm_color = [204, 174, 98]/256;
cond_labels = {'HH','HL','LL','LH'};
leg_ = {'cs-on','cs±45','cs±90','cs±135','cs±180'};
colors_ = cool(5);
dirs = 1:5;

% --- concatenate catch-up counts ---
nTrials_all = S.N_selected_sessions.N;         
N_concat = cell(size(nTrials_all));
for i_group = 1:numel(nTrials_all)
    trials = nTrials_all{i_group};
    n_trials = numel(trials);
    temp = nan(n_trials,4);
    for j = 1:n_trials
        mat = trials{j};     % [8x4]
        temp(j,:) = sum(mat,1,'omitnan');  % sum across directions
    end
    N_concat{i_group} = temp;
end
all_catchup_num = cat(1,N_concat{:});

% --- concatenate total trial counts ---
nTrials_total_all = S.N_selected_sessions.N_total;         
N_total_concat = cell(size(nTrials_total_all));
for i_group = 1:numel(nTrials_total_all)
    trials = nTrials_total_all{i_group};
    n_trials = numel(trials);
    temp = nan(n_trials,4);
    for j = 1:n_trials
        mat = trials{j};     % [8x4]
        temp(j,:) = sum(mat,1,'omitnan');  % sum across directions
    end
    N_total_concat{i_group} = temp;
end
all_total_num = cat(1,N_total_concat{:});

% --- compute per-session percentages ---
num_sessions = size(all_catchup_num,1);
perc_sessions = nan(num_sessions,4);

for cond = 1:4
    for s = 1:num_sessions
        if all_total_num(s,cond) > 0
            perc_sessions(s,cond) = 100 * all_catchup_num(s,cond) ./ all_total_num(s,cond);
        else
            perc_sessions(s,cond) = NaN;   % no trials → undefined %
        end
    end
end

% mean ± SEM, ignore NaNs
perc_mean = mean(perc_sessions,1,'omitnan');
perc_sem  = std(perc_sessions,0,1,'omitnan') ./ sqrt(sum(~isnan(perc_sessions),1));



fig = figure;
for cond = 1:4
    subplot(1,4,cond)
    
    % extract latencies for this condition
    lat_vals = all_lat(all_lat_cond == cond);
 
    % plot histogram of valid values only
    histogram(lat_vals(~isnan(lat_vals)), 'BinWidth', 50, ...
              'FaceColor',[0.3 0.6 0.9], 'FaceAlpha',0.7);
    xlim([0 1000]);   % adjust if needed
    xlabel('Latency (ms)');
    ylabel('Count');
    title(['Latency distribution - ' cond_labels{cond}]);
  
end
sgtitle('Catch-up latencies by condition')
ESN_Beautify_Plot(fig,[12,4])
fig_name1 = 'CS_catchup_population_latency.pdf';
print(fig,fullfile(img_save_path,fig_name1),'-dpdf', '-bestfit')

% -------------- plot no SPE / with RPE ---------

    fig = figure;
    colororder({'k','k'})
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
    
    % --- top left: CS rate vs cond (what you already had) ---
    nexttile
    hold on
    vm_lim = [0,200];
    
    for ii = 1:4
        currect_cond = mean(squeeze(CS_poulation_sac(:, t_ind_sac, :, ii)), 3, 'omitnan');
        Avg_cond = mean(currect_cond, 1, 'omitnan');
        std_cond = std(currect_cond, 1, 'omitnan') ./ sqrt(num_cs);
        yyaxis left
        ylim([0.3 2])
        [hl1, hl2] = boundedline(t_ind_sac-250, Avg_cond, std_cond, 'Color', colors(ii,:), 'alpha');
        hl1.DisplayName = cond_labels{ii};
        hl2.HandleVisibility = 'off';
        if ii == 2 || ii == 4
           hl1.LineStyle = '--';
        end
        if ii == 1
            current_cs_vm = mean(mean(VM_population_sac(:, t_ind_sac, :, :), 4, 'omitnan'), 3, 'omitnan'); 
            yyaxis right
            colororder({'k','k'})
            area(t_ind_sac-250, mean(current_cs_vm,1), ...
                'FaceColor', vm_color, ...
                'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off')
            ylim(vm_lim);
            ylabel('Velocity (deg/sec)')
        end
    end
    legend('Box','off','Location','northwest')
    xline(0, '--', 'HandleVisibility', 'off');
    xlabel('Time relative to sac onset (ms)')
    ylabel('CS rate (Hz)');
    title(['Catch-up sac. (N-{cells} = ' num2str(num_cs) ')'])

    nexttile
    histogram(all_amps, 'BinWidth',0.25, ...
        'FaceColor',[0.3 0.6 0.9], 'FaceAlpha',0.6)
    xlabel('Amplitude (deg)')
    ylabel('Count')
    title('Amplitude distribution (all conditions)')

    nexttile
    bar(1:4, perc_mean,'FaceColor','b'); hold on
    
    % add error bars
    errorbar(1:4, perc_mean, perc_sem, 'k', 'linestyle','none','LineWidth',1);
   
    set(gca,'XTick',1:4,'XTickLabel',cond_labels)
    ylabel('% of trials with catch-up')
    title('Catch-up % by condition')

    ESN_Beautify_Plot(fig,[12,4])
    fig_name1 = 'CS_catchup_population_overview.pdf';
    print(fig,fullfile(img_save_path,fig_name1),'-dpdf', '-bestfit')


end