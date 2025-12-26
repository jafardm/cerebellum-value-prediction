function plot_sac_kinematics_tag1(data_path,pltflag)
    
    if nargin <=1
        pltflag =0;
    end
    JDM_params_funcs;
    load(fullfile(data_path, 'population_data', 'sac_kinematic_prim_sac.mat'));
    img_save_path = fullfile('C:\Users\Jafar\Documents\reward\population_figs\behavior');
    subj_names = params.animal_list;
    num_subjects = numel(subj_names);
    
    labels = {'High','Low'};
    Line_Color =  [190, 37, 43;   % High - DarkRed
                    0, 10, 255] / 255;  % Low  - Blue
    bin_num = 50;

    % helper to append n to labels
    fmt_labels = @(Ns) arrayfun(@(i) sprintf('%s (n=%d)', labels{i}, Ns(i)), 1:numel(labels), 'UniformOutput', false);

    for jj = 1:num_subjects

%% ---- Vigor Force trials -----------
    data_vigor = tag1_vigor{jj,1}(1:4);
    combined_data = { [data_vigor{1}; data_vigor{2}], ...
                      [data_vigor{3}; data_vigor{4}] };
    N_primary = [numel(combined_data{1}), numel(combined_data{2})];

    all_data = cell2mat(data_vigor);
    bin_edges = linspace(min(all_data), max(all_data), bin_num);
    medians = cellfun(@median, combined_data);
    sems    = cellfun(@(x) std(x)/sqrt(numel(x)), combined_data);

    fig = figure;
    tiledlayout(1,4,'TileSpacing','compact');

    % Hist (primary)
    nexttile(1); hold on;
    h_dummy = gobjects(1,2);
    for i = 1:2
        curr = combined_data{i};
        histogram(curr, bin_edges,'Normalization','probability', ...
            'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
        h_dummy(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:),'LineWidth',1.5,'DisplayName',labels{i});
        xline(median(curr),'--','Color',Line_Color(i,:),'LineWidth',1.5);
    end
    legend(h_dummy, fmt_labels(N_primary), 'Location','best','FontSize',8); legend('boxoff');
    ylabel('Probability'); xlabel('Vigor'); title('primary sac.'); box off;

    % Bars (primary)
    nexttile(2); hold on;
    for i = 1:2
        bar(i, medians(i), 'FaceColor', Line_Color(i,:), 'EdgeColor', 'none');
    end
    errorbar(1:2, medians, sems, 'k','LineStyle','none','LineWidth',1);
    xticks(1:2); xticklabels(labels); ylabel('vigor'); title('primary sac.'); box off;
    ylim([min(medians-sems)*0.95, max(medians+sems)*1.08]);


    % Vigor choice
    data_vigor2 = tag1_vigor{jj,1}(5:6);
    N_choice = [numel(data_vigor2{1}), numel(data_vigor2{2})];
    all_data2 = cell2mat(data_vigor2);
    bin_edges2 = linspace(min(all_data2), max(all_data2), bin_num);
    medians2 = cellfun(@median, data_vigor2);
    sems2    = cellfun(@(x) std(x)/sqrt(numel(x)), data_vigor2);

    nexttile(3); hold on;
    h_dummy2 = gobjects(1,2);
    for i = 1:2
        histogram(data_vigor2{i}, bin_edges2,'Normalization','probability', ...
            'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
        h_dummy2(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:),'LineWidth',1.5,'DisplayName',labels{i});
        xline(medians2(i),'--','Color',Line_Color(i,:),'LineWidth',1.5);
    end
    legend(h_dummy2, fmt_labels(N_choice), 'Location','best','FontSize',8); legend('boxoff');
    ylabel('Probability'); xlabel('Vigor'); title('choice'); box off;

    nexttile(4); hold on;
    for i = 1:2
        bar(i, medians2(i), 'FaceColor', Line_Color(i,:), 'EdgeColor','none');
    end
    errorbar(1:2, medians2, sems2, 'k','LineStyle','none','LineWidth',1);
    xticks(1:2); xticklabels(labels); ylabel('vigor'); title('choice'); box off;
    ylim([min(medians2-sems2)*0.95, max(medians2+sems2)*1.08]);
   
    ESN_Beautify_Plot(fig,[10,4]);
    sgtitle([subj_names{jj} '-' 'vigor'],'FontSize',11);
    if pltflag == 1
        print(fig, fullfile(img_save_path,[subj_names{jj} '-vigor-tag1.pdf']), '-dpdf','-bestfit');
    end

% ----------------- Max velocity ---------------
    max_vel_data = tag1_max_vel{jj,1}(1:4);
    combined_data = { [max_vel_data{1}; max_vel_data{2}] , ...
                      [max_vel_data{3}; max_vel_data{4}] };
    N_primary = [numel(combined_data{1}), numel(combined_data{2})];

    all_data = cell2mat(max_vel_data);
    bin_edges = linspace(min(all_data), max(all_data), bin_num);
    medians = cellfun(@median, combined_data);
    sems    = cellfun(@(x) std(x)/sqrt(numel(x)), combined_data);

    fig = figure;
    tiledlayout(1,4,'TileSpacing','compact');

    nexttile(1); hold on;
    h_dummy = gobjects(1,2);
    for i = 1:2
        curr = combined_data{i};
        histogram(curr, bin_edges,'Normalization','probability', ...
            'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
        h_dummy(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:),'LineWidth',1.5);
        xline(median(curr),'--','Color',Line_Color(i,:),'LineWidth',1.5);
    end
    legend(h_dummy, fmt_labels(N_primary), 'Location','best','FontSize',8); legend('boxoff');
    ylabel('Probability'); xlabel('max velocity'); title('primary sac.'); box off;

    nexttile(2); hold on;
    for i = 1:2
        bar(i, medians(i), 'FaceColor', Line_Color(i,:), 'EdgeColor','none');
    end
    errorbar(1:2, medians, sems, 'k','LineStyle','none','LineWidth',1);
    xticks(1:2); xticklabels(labels); ylabel('max velocity'); title('primary sac.'); box off;
    ylim([min(medians-sems)*0.8, max(medians+sems)*1.08]);
   

    data_choice = tag1_max_vel{jj,1}(5:6);
    N_choice = [numel(data_choice{1}), numel(data_choice{2})];
    all_data2 = cell2mat(data_choice);
    bin_edges2 = linspace(min(all_data2), max(all_data2), bin_num);
    medians2 = cellfun(@median, data_choice);
    sems2    = cellfun(@(x) std(x)/sqrt(numel(x)), data_choice);

    nexttile(3); hold on;
    h_dummy = gobjects(1,2);
    for i = 1:2
        histogram(data_choice{i}, bin_edges2,'Normalization','probability', ...
            'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
        h_dummy(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:),'LineWidth',1.5);
        xline(medians2(i),'--','Color',Line_Color(i,:),'LineWidth',1.5);
    end
    legend(h_dummy, fmt_labels(N_choice), 'Location','best','FontSize',8); legend('boxoff');
    ylabel('Probability'); xlabel('max velocity'); title('choice'); box off;

    nexttile(4); hold on;
    for i = 1:2
        bar(i, medians2(i), 'FaceColor', Line_Color(i,:), 'EdgeColor','none');
    end
    errorbar(1:2, medians2, sems2, 'k','LineStyle','none','LineWidth',1);
    xticks(1:2); xticklabels(labels); ylabel('max velocity'); title('choice'); box off;
    ylim([min(medians2-sems2)*0.8, max(medians2+sems2)*1.08]);
  
    ESN_Beautify_Plot(fig,[10,4]);
    sgtitle([subj_names{jj} '-' 'max velocity'],'FontSize',11);
    if pltflag == 1
        print(fig, fullfile(img_save_path,[subj_names{jj} '-max-vel-tag1.pdf']), '-dpdf','-bestfit');
    end

%% ---------- End point error -------------
    data_error = tag1_end_point_error{jj,1}(1:4);
    data_error_mag = cellfun(@(x) sqrt(x(:,1).^2 + x(:,2).^2), data_error, 'UniformOutput', false);
    combined_data = { [data_error_mag{1}; data_error_mag{2}] , ...
                      [data_error_mag{3}; data_error_mag{4}] };
    N_primary = [numel(combined_data{1}), numel(combined_data{2})];

    all_data = cell2mat(data_error_mag);
    bin_edges = linspace(min(all_data)-1, max(all_data), bin_num);
    medians = cellfun(@median, combined_data);
    sems    = cellfun(@(x) std(x)/sqrt(numel(x)), combined_data);

    fig = figure;
    tiledlayout(1,4,'TileSpacing','compact');

    nexttile(1); hold on;
    h_dummy = gobjects(1,2);
    for i = 1:2
        curr = combined_data{i};
        histogram(curr, bin_edges,'Normalization','probability', ...
            'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
        h_dummy(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:),'LineWidth',1.5);
        xline(median(curr),'--','Color',Line_Color(i,:),'LineWidth',1.5);
    end
    legend(h_dummy, fmt_labels(N_primary), 'Location','best','FontSize',8); legend('boxoff');
    ylabel('Probability'); xlabel('End point error'); title('primary sac.'); box off;

    nexttile(2); hold on;
    for i = 1:2
        bar(i, medians(i), 'FaceColor', Line_Color(i,:), 'EdgeColor','none');
    end
    errorbar(1:2, medians, sems, 'k','LineStyle','none','LineWidth',1);
    xticks(1:2); xticklabels(labels); ylabel('End point error'); title('primary sac.'); box off;
    ylim([min(medians-sems)*0.7, max(medians+sems)*1.1]);


    data_choice = tag1_end_point_error{jj,1}(5:6);
    data_error_mag2 = cellfun(@(x) sqrt(x(:,1).^2 + x(:,2).^2), data_choice, 'UniformOutput', false);
    N_choice = [numel(data_error_mag2{1}), numel(data_error_mag2{2})];
    all_data2 = cell2mat(data_error_mag2);
    bin_edges2 = linspace(min(all_data2), max(all_data2), bin_num);
    medians2 = cellfun(@median, data_error_mag2);
    sems2    = cellfun(@(x) std(x)/sqrt(numel(x)), data_error_mag2);

    nexttile(3); hold on;
    h_dummy = gobjects(1,2);
    for i = 1:2
        histogram(data_error_mag2{i}, bin_edges2,'Normalization','probability', ...
            'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
        h_dummy(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:),'LineWidth',1.5);
        xline(medians2(i),'--','Color',Line_Color(i,:),'LineWidth',1.5);
    end
    legend(h_dummy, fmt_labels(N_choice), 'Location','best','FontSize',8); legend('boxoff');
    ylabel('Probability'); xlabel('End point error'); title('choice'); box off;

    nexttile(4); hold on;
    for i = 1:2
        bar(i, medians2(i), 'FaceColor', Line_Color(i,:), 'EdgeColor','none');
    end
    errorbar(1:2, medians2, sems2, 'k','LineStyle','none','LineWidth',1);
    xticks(1:2); xticklabels(labels); ylabel('End point error'); title('choice'); box off;
    ylim([min(medians2-sems2)*0.7, max(medians2+sems2)*1.1]);
  
    ESN_Beautify_Plot(fig,[10,4]);
    sgtitle([subj_names{jj} '-' 'end point error'],'FontSize',11);
    if pltflag == 1
        print(fig, fullfile(img_save_path,[subj_names{jj} '-end-point-error-tag1.pdf']), '-dpdf','-bestfit');
    end

%% ---------- Saccade reaction time -------------
    data_rt = tag1_sac_RT{jj,1}(1:4);
    combined_data = { [data_rt{1}; data_rt{2}] , ...
                      [data_rt{3}; data_rt{4}] };
    N_primary = [numel(combined_data{1}), numel(combined_data{2})];

    all_data = cell2mat(data_rt);
    bin_edges = linspace(min(all_data), max(all_data), bin_num);
    medians = cellfun(@median, combined_data);
    sems    = cellfun(@(x) std(x)/sqrt(numel(x)), combined_data);

    fig = figure;
    tiledlayout(1,4,'TileSpacing','compact');

    nexttile(1); hold on;
    h_dummy = gobjects(1,2);
    for i = 1:2
        curr = combined_data{i};
        histogram(curr, bin_edges,'Normalization','probability', ...
            'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
        h_dummy(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:),'LineWidth',1.5);
        xline(median(curr),'--','Color',Line_Color(i,:),'LineWidth',1.5);
    end
    legend(h_dummy, fmt_labels(N_primary), 'Location','best','FontSize',8); legend('boxoff');
    ylabel('Probability'); xlabel('Reaction time'); title('primary sac.'); box off;

    nexttile(2); hold on;
    for i = 1:2
        bar(i, medians(i), 'FaceColor', Line_Color(i,:), 'EdgeColor','none');
    end
    errorbar(1:2, medians, sems, 'k','LineStyle','none','LineWidth',1);
    xticks(1:2); xticklabels(labels); ylabel('Reaction time (ms)'); title('primary sac.'); box off;
    ylim([min(medians-sems)*0.6, max(medians+sems)*1.08]);
   
    data_choice = tag1_sac_RT{jj,1}(5:6);
    N_choice = [numel(data_choice{1}), numel(data_choice{2})];
    all_data2 = cell2mat(data_choice);
    bin_edges2 = linspace(min(all_data2), max(all_data2), bin_num);
    medians2 = cellfun(@median, data_choice);
    sems2    = cellfun(@(x) std(x)/sqrt(numel(x)), data_choice);

    nexttile(3); hold on;
    h_dummy = gobjects(1,2);
    for i = 1:2
        histogram(data_choice{i}, bin_edges2,'Normalization','probability', ...
            'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
        h_dummy(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:),'LineWidth',1.5);
        xline(medians2(i),'--','Color',Line_Color(i,:),'LineWidth',1.5);
    end
    legend(h_dummy, fmt_labels(N_choice), 'Location','best','FontSize',8); legend('boxoff');
    ylabel('Probability'); xlabel('Reaction time'); title('choice'); box off;

    nexttile(4); hold on;
    for i = 1:2
        bar(i, medians2(i), 'FaceColor', Line_Color(i,:), 'EdgeColor','none');
    end
    errorbar(1:2, medians2, sems2, 'k','LineStyle','none','LineWidth',1);
    xticks(1:2); xticklabels(labels); ylabel('Reaction time (ms)'); title('choice'); box off;
    ylim([min(medians2-sems2)*0.6, max(medians2+sems2)*1.08]);
   
    ESN_Beautify_Plot(fig,[10,4]);
    sgtitle([subj_names{jj} '-' 'Reaction time'],'FontSize',11);
    if pltflag == 1
        print(fig, fullfile(img_save_path,[subj_names{jj} '-Reaction-time-tag1.pdf']), '-dpdf','-bestfit');
    end

%% ------ Eye VM (saccade-aligned traces) ------
    data_vm_sac = tag1_eye_vm_sac{jj,1}(1:4);
    high_traces = cat(2, data_vm_sac{1:2});
    low_traces  = cat(2, data_vm_sac{3:4});
    N_primary = [size(high_traces,2), size(low_traces,2)];
    combined_tr = {high_traces, low_traces};

    t_ind = (-10:40) + 250;

    fig = figure;
    tiledlayout(1,2,'TileSpacing','compact');

    % primary
    nexttile(1); hold on
    L = gobjects(1,2);
    for i = 1:2
        vm = combined_tr{i}(t_ind,:);
        N = size(vm,2);
        [hl1,hl2] = boundedline(t_ind-250, mean(vm,2), std(vm,0,2)/sqrt(max(1,N)), 'Color', Line_Color(i,:));
        L(i) = hl1; L(i).DisplayName = sprintf('%s (n=%d)', labels{i}, N_primary(i));
        hl2.HandleVisibility = 'off';
    end
    legend(L, 'Location','best','FontSize',8); legend('boxoff');
    ylabel('Velocity (deg/s)'); xlabel('Time from sac onset (ms)'); box off; title('primary sac.');

    % choice
    data_vm_ch = tag1_eye_vm_sac{jj,1}(5:6);
    N_choice = [size(data_vm_ch{1},2), size(data_vm_ch{2},2)];

    nexttile(2); hold on
    L = gobjects(1,2);
    for i = 1:2
        vm = data_vm_ch{i}(t_ind,:);
        N = size(vm,2);
        [hl1,hl2] = boundedline(t_ind-250, mean(vm,2), std(vm,0,2)/sqrt(max(1,N)), 'Color', Line_Color(i,:));
        L(i) = hl1; L(i).DisplayName = sprintf('%s (n=%d)', labels{i}, N_choice(i));
        hl2.HandleVisibility = 'off';
    end
    legend(L, 'Location','best','FontSize',8); legend('boxoff');
    ylabel('Velocity (deg/s)'); xlabel('Time from sac onset (ms)'); box off; title('choice');

    sgtitle([subj_names{jj} '-' 'velocity'],'FontSize',11);
    ESN_Beautify_Plot(fig,[10,4]);
    if pltflag == 1
        print(fig, fullfile(img_save_path,[subj_names{jj} '-velocity-tag1.pdf']), '-dpdf','-bestfit');
    end

    end % subject loop
end
