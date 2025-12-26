function plot_sac_kinematics_tag4(data_path)
    JDM_params_funcs;
    load(fullfile(data_path, 'population_data', 'sac_kinematic_corr_sac.mat'));
    img_save_path = fullfile('C:\Users\Jafar\Documents\reward\population_figs\behavior');
    subj_names   = params.animal_list;
    num_subjects = numel(subj_names);
    
    labels = {'HH SPE+','HL SPE+ RPE-','LL SPE+','LH SPE+ RPE+'};
    Line_Color = lines(4);
    bin_num = 50;

    for jj = 1:num_subjects

        %% ---- Vigor (corrective) -----------
        data_vigor = tag4_vigor{jj,1};
        all_data   = cell2mat(data_vigor);
        bin_edges  = linspace(min(all_data), max(all_data), bin_num);
        medians    = cellfun(@median, data_vigor);
        sems       = cellfun(@(x) std(x)/sqrt(length(x)), data_vigor);
        Ns         = cellfun(@numel, data_vigor);     % counts for subplot(2)

        fig = figure; tiledlayout(1,2,'TileSpacing','compact');

        nexttile(1); hold on;
        h_dummy = gobjects(1,4);
        for i = 1:4
            histogram(data_vigor{i}, bin_edges, 'Normalization','probability', ...
                'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
            h_dummy(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:), ...
                'LineWidth',1.5,'DisplayName',labels{i});
            xline(medians(i),'--','Color',Line_Color(i,:),'LineWidth',1.5);
        end
        legend(h_dummy,'Location','best','FontSize',8); legend('boxoff');
        ylabel('Probability'); xlabel('Vigor'); title('corrective sac.'); box off;

        % (2) bars (legend WITH n)
        nexttile(2); hold on;
        h_dummy2 = gobjects(1,4);
        for i = 1:4
            bar(i, medians(i), 'FaceColor', Line_Color(i,:), 'EdgeColor','none');
            h_dummy2(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:), ...
                'LineWidth',1.5,'DisplayName', sprintf('%s (n=%d)', labels{i}, Ns(i)));
        end
        errorbar(1:4, medians, sems, 'k', 'LineStyle','none','LineWidth',1);
        xticks(1:4); xticklabels(labels); ylabel('vigor'); title('corrective sac.'); box off;
        ylim([mean(medians)*0.95, max(medians + sems)*1.01]);
        legend(h_dummy2,'Location','best','FontSize',8); legend('boxoff');

        ESN_Beautify_Plot(fig, [10, 4]);
        sgtitle([subj_names{jj} '-' 'vigor'],'FontSize',11);
        print(fig, fullfile(img_save_path,[subj_names{jj} '-vigor-tag4.pdf']), '-dpdf','-bestfit');

        %% ---- Max velocity (corrective) -----------
        max_vel_data = tag4_max_vel{jj,1};
        all_data     = cell2mat(max_vel_data);
        bin_edges    = linspace(min(all_data), max(all_data), bin_num);
        medians      = cellfun(@median, max_vel_data);
        sems         = cellfun(@(x) std(x)/sqrt(length(x)), max_vel_data);
        Ns           = cellfun(@numel, max_vel_data);

        fig = figure; tiledlayout(1,2,'TileSpacing','compact');

        % (1) histogram (legend WITHOUT n)
        nexttile(1); hold on;
        h_dummy = gobjects(1,4);
        for i = 1:4
            histogram(max_vel_data{i}, bin_edges, 'Normalization','probability', ...
                'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
            xline(medians(i),'--','Color',Line_Color(i,:),'LineWidth',1.5);
            h_dummy(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:), ...
                'LineWidth',1.5,'DisplayName',labels{i});
        end
        legend(h_dummy,'Location','best','FontSize',8); legend('boxoff');
        ylabel('Probability'); xlabel('Max-velocity'); title('corrective sac.'); box off;

        % (2) bars (legend WITH n)
        nexttile(2); hold on;
        h_dummy2 = gobjects(1,4);
        for i = 1:4
            bar(i, medians(i), 'FaceColor', Line_Color(i,:), 'EdgeColor','none');
            h_dummy2(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:), ...
                'LineWidth',1.5,'DisplayName', sprintf('%s (n=%d)', labels{i}, Ns(i)));
        end
        errorbar(1:4, medians, sems, 'k', 'LineStyle','none','LineWidth',1);
        xticks(1:4); xticklabels(labels); ylabel('velocity (deg/sec)'); title('corrective sac.'); box off;
        ylim([mean(medians)*0.5, max(medians + sems)*1.01]);
        legend(h_dummy2,'Location','best','FontSize',8); legend('boxoff');

        ESN_Beautify_Plot(fig, [10, 4]);
        sgtitle([subj_names{jj} '-' 'max-velocity'],'FontSize',11);
        print(fig, fullfile(img_save_path,[subj_names{jj} '-max-vel-tag4.pdf']), '-dpdf','-bestfit');

        %% ---- End point error (corrective) -----------
        data_error    = tag4_end_point_error{jj,1};
        data_error_mag = cellfun(@(x) sqrt(x(:,1).^2 + x(:,2).^2), data_error, 'UniformOutput', false);
        all_data      = cell2mat(data_error_mag);
        bin_edges     = linspace(min(all_data), max(all_data), bin_num);
        medians       = cellfun(@median, data_error_mag);
        sems          = cellfun(@(x) std(x)/sqrt(length(x)), data_error_mag);
        Ns            = cellfun(@numel, data_error_mag);

        fig = figure; tiledlayout(1,2,'TileSpacing','compact');

        % (1) histogram 
        nexttile(1); hold on;
        h_dummy = gobjects(1,4);
        for i = 1:4
            histogram(data_error_mag{i}, bin_edges, 'Normalization','probability', ...
                'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
            xline(medians(i),'--','Color',Line_Color(i,:),'LineWidth',1.5);
            h_dummy(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:), ...
                'LineWidth',1.5,'DisplayName',labels{i});
        end
        legend(h_dummy,'Location','best','FontSize',8); legend('boxoff');
        ylabel('Probability'); xlabel('End point error'); title('corrective sac.'); box off;

        % (2) bars 
        nexttile(2); hold on;
        h_dummy2 = gobjects(1,4);
        for i = 1:4
            bar(i, medians(i), 'FaceColor', Line_Color(i,:), 'EdgeColor','none');
            h_dummy2(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:), ...
                'LineWidth',1.5,'DisplayName', sprintf('%s (n=%d)', labels{i}, Ns(i)));
        end
        errorbar(1:4, medians, sems, 'k', 'LineStyle','none','LineWidth',1);
        xticks(1:4); xticklabels(labels); ylabel('end point error (deg)'); title('corrective sac.'); box off;
        ylim([mean(medians)*0.5, max(medians + sems)*1.01]);
        legend(h_dummy2,'Location','best','FontSize',8); legend('boxoff');

        ESN_Beautify_Plot(fig, [10, 4]);
        sgtitle([subj_names{jj} '-' 'end point er.'],'FontSize',11);
        print(fig, fullfile(img_save_path,[subj_names{jj} '-end-point-error-tag4.pdf']), '-dpdf','-bestfit');

        %% ---- Reaction time (corrective) -----------
        data_rt   = tag4_sac_RT{jj,1};
        all_data  = cell2mat(data_rt);
        bin_edges = linspace(min(all_data), max(all_data), bin_num);
        medians   = cellfun(@median, data_rt);
        sems      = cellfun(@(x) std(x)/sqrt(length(x)), data_rt);
        Ns        = cellfun(@numel, data_rt);

        fig = figure; tiledlayout(1,2,'TileSpacing','compact');

        % (1) histogram 
        nexttile(1); hold on;
        h_dummy = gobjects(1,4);
        for i = 1:4
            histogram(data_rt{i}, bin_edges, 'Normalization','probability', ...
                'FaceColor',Line_Color(i,:), 'EdgeColor','none','FaceAlpha',0.5);
            xline(medians(i),'--','Color',Line_Color(i,:),'LineWidth',1.5);
            h_dummy(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:), ...
                'LineWidth',1.5,'DisplayName',labels{i});
        end
        legend(h_dummy,'Location','best','FontSize',8); legend('boxoff');
        ylabel('Probability'); xlabel('Reaction time'); title('corrective sac.'); box off;

        % (2) bars 
        nexttile(2); hold on;
        h_dummy2 = gobjects(1,4);
        for i = 1:4
            bar(i, medians(i), 'FaceColor', Line_Color(i,:), 'EdgeColor','none');
            h_dummy2(i) = plot(NaN,NaN,'-','Color',Line_Color(i,:), ...
                'LineWidth',1.5,'DisplayName', sprintf('%s (n=%d)', labels{i}, Ns(i)));
        end
        errorbar(1:4, medians, sems, 'k', 'LineStyle','none','LineWidth',1);
        xticks(1:4); xticklabels(labels);
        ylabel('Reaction time (ms)'); title('corrective sac.'); box off;
        ylim([mean(medians)*0.5, max(medians + sems)*1.01]);
        legend(h_dummy2,'Location','best','FontSize',8); legend('boxoff');

        ESN_Beautify_Plot(fig, [10, 4]);
        sgtitle([subj_names{jj} '-' 'Reaction time'],'FontSize',11);
        print(fig, fullfile(img_save_path,[subj_names{jj} '-Reaction-time-tag4.pdf']), '-dpdf','-bestfit');

    end
end
