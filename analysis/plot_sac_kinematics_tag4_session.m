function plot_sac_kinematics_tag4_session(data_path)
    JDM_params_funcs;
    
    % NOTE: using the *session* version
    load(fullfile(data_path, 'population_data', 'sac_kinematic_corr_sac_session.mat'));
    
    img_save_path = fullfile('C:\Users\Jafar\Documents\reward\population_figs\behavior');
    subj_names   = params.animal_list;
    num_subjects = numel(subj_names);
    
    labels = {'HH:SPE+','HL:SPE+ RPE-','LL:SPE+','LH:SPE+ RPE+'};
    Line_Color = lines(4);
    bin_num = 50;

    for jj = 1:num_subjects
        fprintf('Subject %d / %d\n', jj, num_subjects);

        % =================================================================
        % 1) VIGOR
        % =================================================================
        data_vigor_sess = tag4_vigor{jj};   % {nSess x 1}, each {4 x 1}
        nSess = numel(data_vigor_sess);

        % ---- aggregate trials across sessions per condition ----
        cat_HH = []; cat_HL = []; cat_LL = []; cat_LH = [];

        % ---- session-level means for violin in subplot 2 & RPE in subplot 3 ----
        HH_sess_mean = NaN(nSess,1);
        HL_sess_mean = NaN(nSess,1);
        LL_sess_mean = NaN(nSess,1);
        LH_sess_mean = NaN(nSess,1);

        for ss = 1:nSess
            sess_data = data_vigor_sess{ss};
            if isempty(sess_data) || numel(sess_data) < 4, continue; end
            HH = sess_data{1};
            HL = sess_data{2};
            LL = sess_data{3};
            LH = sess_data{4};

            % all trials for hist
            cat_HH = [cat_HH; HH(:)];
            cat_HL = [cat_HL; HL(:)];
            cat_LL = [cat_LL; LL(:)];
            cat_LH = [cat_LH; LH(:)];

            % session means
            if ~isempty(HH), HH_sess_mean(ss) = mean(HH); end
            if ~isempty(HL), HL_sess_mean(ss) = mean(HL); end
            if ~isempty(LL), LL_sess_mean(ss) = mean(LL); end
            if ~isempty(LH), LH_sess_mean(ss) = mean(LH); end
        end

        all_data = [cat_HH; cat_HL; cat_LL; cat_LH];
        if isempty(all_data)
            warning('No vigor data for subject %d, skipping.', jj);
        else
            bin_edges = linspace(min(all_data), max(all_data), bin_num);
            medians   = [median(cat_HH), median(cat_HL), ...
                         median(cat_LL), median(cat_LH)];

            fig = figure; tiledlayout(1,3,'TileSpacing','compact');

            % ---- (1) histogram per condition (unchanged) ----
            nexttile(1); hold on;
            histogram(cat_HH, bin_edges,'FaceColor',Line_Color(1,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_HL, bin_edges,'FaceColor',Line_Color(2,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_LL, bin_edges,'FaceColor',Line_Color(3,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_LH, bin_edges,'FaceColor',Line_Color(4,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            xline(medians(1),'--','Color',Line_Color(1,:),'LineWidth',1.5);
            xline(medians(2),'--','Color',Line_Color(2,:),'LineWidth',1.5);
            xline(medians(3),'--','Color',Line_Color(3,:),'LineWidth',1.5);
            xline(medians(4),'--','Color',Line_Color(4,:),'LineWidth',1.5);
            legend(labels); legend('boxoff');
            ylabel('Probability'); xlabel('Vigor');
            title('Corrective sac.'); box off;

            % ---- (2) violin plot of SESSION MEANS per condition ----
            nexttile(2); hold on;
            vHH = HH_sess_mean(~isnan(HH_sess_mean));
            vHL = HL_sess_mean(~isnan(HL_sess_mean));
            vLL = LL_sess_mean(~isnan(LL_sess_mean));
            vLH = LH_sess_mean(~isnan(LH_sess_mean));

            vals   = [vHH; vHL; vLL; vLH];
            groups = [ones(numel(vHH),1); ...
                      2*ones(numel(vHL),1); ...
                      3*ones(numel(vLL),1); ...
                      4*ones(numel(vLH),1)];
            violinplot(vals, groups);
            xticks(1:4); xticklabels(labels);
            ylabel('Mean vigor per session'); title('Session means'); box off;

            % ---- (3) RPE differences per session: violin of Δmeans ----
            RPE_minus_sess = HL_sess_mean - HH_sess_mean;   % HL - HH
            RPE_plus_sess  = LH_sess_mean - LL_sess_mean;   % LH - LL

            vMinus = RPE_minus_sess(~isnan(RPE_minus_sess));
            vPlus  = RPE_plus_sess(~isnan(RPE_plus_sess));

            vals   = [vMinus; vPlus];
            groups = [ones(numel(vMinus),1); 2*ones(numel(vPlus),1)];

            nexttile(3); hold on;
            violinplot(vals, groups);
            xticks([1 2]); xticklabels({'RPE- (HL-HH)','RPE+ (LH-LL)'});
            ylabel('\Delta Vigor (per session)');
            title('RPE differences'); box off;

            ESN_Beautify_Plot(fig, [10, 4]);
            sgtitle([subj_names{jj} '-' 'vigor'],'FontSize',11);
            print(fig, fullfile(img_save_path,[subj_names{jj} '-vigor-tag4_session.pdf']), ...
                 '-dpdf','-bestfit');
        end


        % =================================================================
        % 2) MAX VELOCITY
        % =================================================================
        max_vel_sess = tag4_max_vel{jj};   % {nSess x 1}, each {4 x 1}

        cat_HH = []; cat_HL = []; cat_LL = []; cat_LH = [];
        HH_sess_mean = NaN(nSess,1);
        HL_sess_mean = NaN(nSess,1);
        LL_sess_mean = NaN(nSess,1);
        LH_sess_mean = NaN(nSess,1);

        for ss = 1:nSess
            sess_data = max_vel_sess{ss};
            if isempty(sess_data) || numel(sess_data) < 4, continue; end
            HH = sess_data{1};
            HL = sess_data{2};
            LL = sess_data{3};
            LH = sess_data{4};

            cat_HH = [cat_HH; HH(:)];
            cat_HL = [cat_HL; HL(:)];
            cat_LL = [cat_LL; LL(:)];
            cat_LH = [cat_LH; LH(:)];

            if ~isempty(HH), HH_sess_mean(ss) = mean(HH); end
            if ~isempty(HL), HL_sess_mean(ss) = mean(HL); end
            if ~isempty(LL), LL_sess_mean(ss) = mean(LL); end
            if ~isempty(LH), LH_sess_mean(ss) = mean(LH); end
        end

        all_data = [cat_HH; cat_HL; cat_LL; cat_LH];
        if isempty(all_data)
            warning('No max-vel data for subject %d, skipping.', jj);
        else
            bin_edges = linspace(min(all_data), max(all_data), bin_num);
            medians   = [median(cat_HH), median(cat_HL), ...
                         median(cat_LL), median(cat_LH)];

            fig = figure; tiledlayout(1,3,'TileSpacing','compact');

            % (1) histogram
            nexttile(1); hold on;
            histogram(cat_HH, bin_edges,'FaceColor',Line_Color(1,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_HL, bin_edges,'FaceColor',Line_Color(2,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_LL, bin_edges,'FaceColor',Line_Color(3,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_LH, bin_edges,'FaceColor',Line_Color(4,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            xline(medians(1),'--','Color',Line_Color(1,:),'LineWidth',1.5);
            xline(medians(2),'--','Color',Line_Color(2,:),'LineWidth',1.5);
            xline(medians(3),'--','Color',Line_Color(3,:),'LineWidth',1.5);
            xline(medians(4),'--','Color',Line_Color(4,:),'LineWidth',1.5);
            legend(labels); legend('boxoff');
            ylabel('Probability'); xlabel('Max-velocity');
            title('Corrective sac.'); box off;

            % (2) violin of SESSION MEANS
            nexttile(2); hold on;
            vHH = HH_sess_mean(~isnan(HH_sess_mean));
            vHL = HL_sess_mean(~isnan(HL_sess_mean));
            vLL = LL_sess_mean(~isnan(LL_sess_mean));
            vLH = LH_sess_mean(~isnan(LH_sess_mean));

            vals   = [vHH; vHL; vLL; vLH];
            groups = [ones(numel(vHH),1); ...
                      2*ones(numel(vHL),1); ...
                      3*ones(numel(vLL),1); ...
                      4*ones(numel(vLH),1)];
            violinplot(vals, groups);
            xticks(1:4); xticklabels(labels);
            ylabel('Mean velocity (deg/s) per session'); title('Session means'); box off;

            % (3) RPE differences per session
            RPE_minus_sess = HL_sess_mean - HH_sess_mean;
            RPE_plus_sess  = LH_sess_mean - LL_sess_mean;

            vMinus = RPE_minus_sess(~isnan(RPE_minus_sess));
            vPlus  = RPE_plus_sess(~isnan(RPE_plus_sess));

            vals   = [vMinus; vPlus];
            groups = [ones(numel(vMinus),1); 2*ones(numel(vPlus),1)];

            nexttile(3); hold on;
            violinplot(vals, groups);
            xticks([1 2]); xticklabels({'RPE- (HL-HH)','RPE+ (LH-LL)'});
            ylabel('\Delta velocity (deg/s, per session)');
            title('RPE differences'); box off;

            ESN_Beautify_Plot(fig, [10, 4]);
            sgtitle([subj_names{jj} '-' 'max-velocity'],'FontSize',11);
            print(fig, fullfile(img_save_path,[subj_names{jj} '-max-vel-tag4_session.pdf']), ...
                  '-dpdf','-bestfit');
        end


        % =================================================================
        % 3) ENDPOINT ERROR (magnitude)
        % =================================================================
        data_err_sess = tag4_end_point_error{jj};   % {nSess x 1}, each {4 x 1 of [N x 2]}
        
        cat_HH = []; cat_HL = []; cat_LL = []; cat_LH = [];
        HH_sess_mean = NaN(nSess,1);
        HL_sess_mean = NaN(nSess,1);
        LL_sess_mean = NaN(nSess,1);
        LH_sess_mean = NaN(nSess,1);

        for ss = 1:nSess
            sess_data = data_err_sess{ss};
            if isempty(sess_data) || numel(sess_data) < 4, continue; end
            % raw endpoint vectors
            eHH = sess_data{1};
            eHL = sess_data{2};
            eLL = sess_data{3};
            eLH = sess_data{4};

            mag_HH = sqrt(eHH(:,1).^2 + eHH(:,2).^2);
            mag_HL = sqrt(eHL(:,1).^2 + eHL(:,2).^2);
            mag_LL = sqrt(eLL(:,1).^2 + eLL(:,2).^2);
            mag_LH = sqrt(eLH(:,1).^2 + eLH(:,2).^2);

            cat_HH = [cat_HH; mag_HH(:)];
            cat_HL = [cat_HL; mag_HL(:)];
            cat_LL = [cat_LL; mag_LL(:)];
            cat_LH = [cat_LH; mag_LH(:)];

            if ~isempty(mag_HH), HH_sess_mean(ss) = mean(mag_HH); end
            if ~isempty(mag_HL), HL_sess_mean(ss) = mean(mag_HL); end
            if ~isempty(mag_LL), LL_sess_mean(ss) = mean(mag_LL); end
            if ~isempty(mag_LH), LH_sess_mean(ss) = mean(mag_LH); end
        end

        all_data = [cat_HH; cat_HL; cat_LL; cat_LH];
        if isempty(all_data)
            warning('No endpoint error data for subject %d, skipping.', jj);
        else
            bin_edges = linspace(min(all_data), max(all_data), bin_num);
            medians   = [median(cat_HH), median(cat_HL), ...
                         median(cat_LL), median(cat_LH)];

            fig = figure; tiledlayout(1,3,'TileSpacing','compact');

            % (1) histogram
            nexttile(1); hold on;
            histogram(cat_HH, bin_edges,'FaceColor',Line_Color(1,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_HL, bin_edges,'FaceColor',Line_Color(2,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_LL, bin_edges,'FaceColor',Line_Color(3,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_LH, bin_edges,'FaceColor',Line_Color(4,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            xline(medians(1),'--','Color',Line_Color(1,:),'LineWidth',1.5);
            xline(medians(2),'--','Color',Line_Color(2,:),'LineWidth',1.5);
            xline(medians(3),'--','Color',Line_Color(3,:),'LineWidth',1.5);
            xline(medians(4),'--','Color',Line_Color(4,:),'LineWidth',1.5);
            legend(labels); legend('boxoff');
            ylabel('Probability'); xlabel('End point error (deg)');
            title('Corrective sac.'); box off;

            % (2) violin of SESSION MEANS
            nexttile(2); hold on;
            vHH = HH_sess_mean(~isnan(HH_sess_mean));
            vHL = HL_sess_mean(~isnan(HL_sess_mean));
            vLL = LL_sess_mean(~isnan(LL_sess_mean));
            vLH = LH_sess_mean(~isnan(LH_sess_mean));

            vals   = [vHH; vHL; vLL; vLH];
            groups = [ones(numel(vHH),1); ...
                      2*ones(numel(vHL),1); ...
                      3*ones(numel(vLL),1); ...
                      4*ones(numel(vLH),1)];
            violinplot(vals, groups);
            xticks(1:4); xticklabels(labels);
            ylabel('Mean endpoint error (deg) per session');
            title('Session means'); box off;

            % (3) RPE differences per session
            RPE_minus_sess = HL_sess_mean - HH_sess_mean;
            RPE_plus_sess  = LH_sess_mean - LL_sess_mean;

            vMinus = RPE_minus_sess(~isnan(RPE_minus_sess));
            vPlus  = RPE_plus_sess(~isnan(RPE_plus_sess));

            vals   = [vMinus; vPlus];
            groups = [ones(numel(vMinus),1); 2*ones(numel(vPlus),1)];

            nexttile(3); hold on;
            violinplot(vals, groups);
            xticks([1 2]); xticklabels({'RPE- (HL-HH)','RPE+ (LH-LL)'});
            ylabel('\Delta endpoint error (deg, per session)');
            title('RPE differences'); box off;

            ESN_Beautify_Plot(fig, [10, 4]);
            sgtitle([subj_names{jj} '-' 'end point er.'],'FontSize',11);
            print(fig, fullfile(img_save_path,[subj_names{jj} '-end-point-error-tag4_session.pdf']), ...
                  '-dpdf','-bestfit');
        end


        % =================================================================
        % 4) REACTION TIME
        % =================================================================
        data_rt_sess = tag4_sac_RT{jj};   % {nSess x 1}

        cat_HH = []; cat_HL = []; cat_LL = []; cat_LH = [];
        HH_sess_mean = NaN(nSess,1);
        HL_sess_mean = NaN(nSess,1);
        LL_sess_mean = NaN(nSess,1);
        LH_sess_mean = NaN(nSess,1);

        for ss = 1:nSess
            sess_data = data_rt_sess{ss};
            if isempty(sess_data) || numel(sess_data) < 4, continue; end

            HH = sess_data{1};
            HL = sess_data{2};
            LL = sess_data{3};
            LH = sess_data{4};

            cat_HH = [cat_HH; HH(:)];
            cat_HL = [cat_HL; HL(:)];
            cat_LL = [cat_LL; LL(:)];
            cat_LH = [cat_LH; LH(:)];

            if ~isempty(HH), HH_sess_mean(ss) = mean(HH); end
            if ~isempty(HL), HL_sess_mean(ss) = mean(HL); end
            if ~isempty(LL), LL_sess_mean(ss) = mean(LL); end
            if ~isempty(LH), LH_sess_mean(ss) = mean(LH); end
        end

        all_data = [cat_HH; cat_HL; cat_LL; cat_LH];
        if isempty(all_data)
            warning('No RT data for subject %d, skipping.', jj);
        else
            bin_edges = linspace(min(all_data), max(all_data), bin_num);
            medians   = [median(cat_HH), median(cat_HL), ...
                         median(cat_LL), median(cat_LH)];

            fig = figure; tiledlayout(1,3,'TileSpacing','compact');

            % (1) histogram
            nexttile(1); hold on;
            histogram(cat_HH, bin_edges,'FaceColor',Line_Color(1,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_HL, bin_edges,'FaceColor',Line_Color(2,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_LL, bin_edges,'FaceColor',Line_Color(3,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            histogram(cat_LH, bin_edges,'FaceColor',Line_Color(4,:), ...
                      'FaceAlpha',0.5,'EdgeColor','none');
            xline(medians(1),'--','Color',Line_Color(1,:),'LineWidth',1.5);
            xline(medians(2),'--','Color',Line_Color(2,:),'LineWidth',1.5);
            xline(medians(3),'--','Color',Line_Color(3,:),'LineWidth',1.5);
            xline(medians(4),'--','Color',Line_Color(4,:),'LineWidth',1.5);
            legend(labels); legend('boxoff');
            ylabel('Probability'); xlabel('Reaction time (ms)');
            title('Corrective sac.'); box off;

            % (2) violin per condition (SESSION MEANS)
            nexttile(2); hold on;
            vHH = HH_sess_mean(~isnan(HH_sess_mean));
            vHL = HL_sess_mean(~isnan(HL_sess_mean));
            vLL = LL_sess_mean(~isnan(LL_sess_mean));
            vLH = LH_sess_mean(~isnan(LH_sess_mean));

            vals   = [vHH; vHL; vLL; vLH];
            groups = [ones(numel(vHH),1); ...
                      2*ones(numel(vHL),1); ...
                      3*ones(numel(vLL),1); ...
                      4*ones(numel(vLH),1)];
            violinplot(vals, groups);
            xticks(1:4); xticklabels(labels);
            ylabel('Mean RT (ms) per session');
            title('Session means'); box off;

            % (3) RPE differences per session
            RPE_minus_sess = HL_sess_mean - HH_sess_mean;
            RPE_plus_sess  = LH_sess_mean - LL_sess_mean;

            vMinus = RPE_minus_sess(~isnan(RPE_minus_sess));
            vPlus  = RPE_plus_sess(~isnan(RPE_plus_sess));

            vals   = [vMinus; vPlus];
            groups = [ones(numel(vMinus),1); 2*ones(numel(vPlus),1)];

            nexttile(3); hold on;
            violinplot(vals, groups);
            xticks([1 2]); xticklabels({'RPE- (HL-HH)','RPE+ (LH-LL)'});
            ylabel('\Delta RT (ms, per session)');
            title('RPE differences'); box off;

            ESN_Beautify_Plot(fig, [10, 4]);
            sgtitle([subj_names{jj} '-' 'Reaction time'],'FontSize',11);
            print(fig, fullfile(img_save_path,[subj_names{jj} '-Reaction-time-tag4_session.pdf']), ...
                  '-dpdf','-bestfit');
        end

    end % jj
end
