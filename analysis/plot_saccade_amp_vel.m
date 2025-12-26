function plot_saccade_amp_vel(data_path)
JDM_params_funcs;
load(fullfile(data_path, 'population_data', 'behave_data.mat'));
alpha = .9;
sz    = .2; 
amp_edges = params.pop.sac.amp_edges;
vel_edges = params.pop.sac.vel_edges;
acc_edges = params.pop.sac.acc_edges;
subject_names = params.animal_list;
fig = figure;
tag_cond = {10,[1,4,6,7]};
clr      = {'b','r','g','k'};

n_col = 4;
n_row = 2;

lim_acc = [-2,3.5];

for counter_animal = 1:2
    current_sac_amp = data_behav.sac_amp{counter_animal};
    current_sac_vel = data_behav.sac_vel{counter_animal};
    current_sac_tag = data_behav.sac_tag{counter_animal};
    current_sac_acc = data_behav.sac_acc{counter_animal};
    current_sac_dec = data_behav.sac_dec{counter_animal};

    b_ms = gmregress(log10(current_sac_amp),log10(current_sac_vel));

    b_dec = gmregress(log10(current_sac_vel),log10(current_sac_dec));
    b_acc = gmregress(log10(current_sac_vel(current_sac_acc>0)),log10(current_sac_acc(current_sac_acc>0)));
    b_acc_dec = gmregress(log10(current_sac_dec(current_sac_acc>0)),log10(current_sac_acc(current_sac_acc>0)));

    for counter_cond = 1:numel(tag_cond)
        current_cond = tag_cond{counter_cond};
        ind_ = ismember(current_sac_tag,current_cond);
        current_clr = clr{counter_cond};

        subplot(n_row,n_col,n_col*(counter_animal-1)+1)
        title(subject_names{counter_animal})
        scatter(log10(current_sac_amp(ind_)),log10(current_sac_vel(ind_)),sz,'filled',current_clr,'MarkerFaceAlpha',alpha);
        hold on;
        x_ = [-1,1.5];
        y_ = 1:3;
        if counter_cond == numel(tag_cond)
            plot(x_,x_*b_ms(2)+b_ms(1),'k');
        end
        xline(log10(amp_edges(2:end-1)),'--k')
        yline(log10(vel_edges(2:end-1)),'--k')
        xticks(x_(1):x_(2));
        xticklabels(10.^(x_(1):x_(2)));
        yticks(y_);
        yticklabels(10.^y_);
        xlabel('amp (deg)')
        ylabel('vel (deg/s)')
        xlim([-1,1.5])
        ylim(y_([1,end]))

        subplot(n_row,n_col,n_col*(counter_animal-1)+2)
        scatter(log10(current_sac_vel(ind_)),log10(current_sac_dec(ind_)),sz,'filled',current_clr,'MarkerFaceAlpha',alpha);
        hold on;
        x_ = 1:3;
        y_ = lim_acc(1):lim_acc(end);
        xline(log10(vel_edges(2:end-1)),'--k')
        yline(log10(acc_edges(2:end-1)),'--k')
        xticks(x_);
        xticklabels(10.^x_);
        if counter_cond == numel(tag_cond)
            plot(x_,x_*b_dec(2)+b_dec(1),'k');
        end
        yticks(y_);
        yticklabels(10.^y_);
        xlabel('vel (deg/s)')
        ylabel('dec (deg/s^2)')
        xlim([1,3])
        ylim(lim_acc)

        subplot(n_row,n_col,n_col*(counter_animal-1)+3)
        scatter(log10(current_sac_vel(ind_)),log10(current_sac_acc(ind_)),sz,'filled',current_clr,'MarkerFaceAlpha',alpha);
        hold on;
        x_ = 1:3;
        y_ = lim_acc(1):lim_acc(end);
        xline(log10(vel_edges(2:end-1)),'--k')
        yline(log10(acc_edges(2:end-1)),'--k')
        xticks(x_);
        xticklabels(10.^x_);
        if counter_cond == numel(tag_cond)
            plot(x_,x_*b_acc(2)+b_acc(1),'k');
        end
        yticks(y_);
        yticklabels(10.^y_);
        xlabel('vel (deg/s)')
        ylabel('acc (deg/s^2)')
        xlim([x_(1) x_(end)])
        ylim(lim_acc)

        subplot(n_row,n_col,n_col*(counter_animal-1)+4)
        scatter(log10(current_sac_dec(ind_)),log10(current_sac_acc(ind_)),sz,'filled',current_clr,'MarkerFaceAlpha',alpha);
        hold on;
        x_ = -1:3;
        y_ = lim_acc(1):lim_acc(end);
        xline(log10(acc_edges(2:end-1)),'--k')
        yline(log10(acc_edges(2:end-1)),'--k')
        xticks(x_);
        xticklabels(10.^x_);
        if counter_cond == numel(tag_cond)
            plot(x_,x_*b_acc_dec(2)+b_acc_dec(1),'k');
        end
        yticks(y_);
        yticklabels(10.^y_);
        xticks(y_);
        xticklabels(10.^y_);
        xlabel('dec (deg/s^2)')
        ylabel('acc (deg/s^2)')
        axis equal
        xlim(lim_acc)
        ylim(lim_acc)

    end
end

ESN_Beautify_Plot(fig,[20,8]);
end