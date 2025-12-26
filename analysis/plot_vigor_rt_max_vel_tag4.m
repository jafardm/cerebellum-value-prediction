function plot_vigor_rt_max_vel_tag4(data_path)
JDM_params_funcs;
load(fullfile(data_path, 'population_data', 'sac_kinematic_corr_sac.mat'));
img_save_path = fullfile('C:\Users\Jafar\Documents\reward\population_figs\behavior');
subj_names = params.animal_list;
num_subjects = size(subj_names,2);
labels = {'HH SPE+', 'HL SPE+ RPE-', 'LL SPE+', 'LH SPE+ RPE+'};
Line_Color = lines(4);

for jj = 1:num_subjects
    %-----------------   vigor ---------------
    data_vigor = tag4_vigor{jj,1};
    all_data = [];  % Flatten data and assign group labels
    group_labels = {};
    
    for i = 1:length(data_vigor)
        this_data = data_vigor{i};
        all_data = [all_data; this_data];  % concatenate all vigor values
        group_labels = [group_labels; repmat(labels(i), length(this_data), 1)];
    end
    % Convert group_labels to a categorical array for plotting
    group_labels = categorical(group_labels, labels, 'Ordinal', true);
    
    
    % Plot
    fig = figure;
    v = violinplot(all_data, group_labels);
    hold on
    cats = categories(group_labels);
    for i = 1:length(cats)
        % Get all data in this group
        this_group_data = all_data(group_labels == cats{i});
        med = median(this_group_data);
        % Draw median line across the width of the violin
        plot([i - 0.2, i + 0.2], [med, med], 'k-', 'LineWidth', 1.5);
    end
    hold off
    
    % Labels and formatting
    ylabel('Vigor');
    xtickangle(45);
    title('Vigor corrective sac');
    box off
    ESN_Beautify_Plot(fig, [10, 4]);
    %print(fig, fullfile(img_save_path, 'vigor-tag4.pdf'), '-dpdf', '-bestfit');
    
    %-----------------   max velocity   ---------------
    
    max_vel_data = tag4_max_vel{jj,1};
    all_data = [];
    group_labels = {};
    
    for i = 1:length(max_vel_data)
        this_data = max_vel_data{i};
        all_data = [all_data; this_data];  % concatenate all vigor values
        group_labels = [group_labels; repmat(labels(i), length(this_data), 1)];
    end
    
    % Convert group_labels to a categorical array for plotting
    group_labels = categorical(group_labels, labels, 'Ordinal', true);
    
    fig = figure;
    v = violinplot(all_data, group_labels);
    
    hold on
    cats = categories(group_labels);
    for i = 1:length(cats)
        this_group_data = all_data(group_labels == cats{i});
        med = median(this_group_data);
        plot([i - 0.2, i + 0.2], [med, med], 'k-', 'LineWidth', 1.5);
    end
    hold off
    
    for i = 1:length(v)
        v(i).ViolinColor = {Line_Color(i, :)};
    end
    xtickangle(45);
    ylabel('Velocity (deg/sec)');
    title('Max velocity corrective sac');
    box off
    ESN_Beautify_Plot(fig,[10,4]);
    %print(fig,fullfile(img_save_path,'Max_vel-tag4.pdf'),'-dpdf','-bestfit')
    
    %-----------------   Reaction times ---------------
    
    data_RT = tag4_sac_RT{jj,1};
    all_data = [];  % Flatten data and assign group labels
    group_labels = {};
    
    for i = 1:length(data_RT)
        this_data = data_RT{i};
        all_data = [all_data; this_data];  % concatenate all vigor values
        group_labels = [group_labels; repmat(labels(i), length(this_data), 1)];
    end
    
    % Convert group_labels to a categorical array for plotting
    group_labels = categorical(group_labels, labels, 'Ordinal', true);
    fig = figure;
    v = violinplot(all_data, group_labels);
    
    hold on
    cats = categories(group_labels);
    for i = 1:length(cats)
        this_group_data = all_data(group_labels == cats{i});
        med = median(this_group_data);
        plot([i - 0.2, i + 0.2], [med, med], 'k-', 'LineWidth', 1.5);
    end
    hold off
    
    for i = 1:length(v)
        v(i).ViolinColor = {Line_Color(i, :)};
    end
    
    xtickangle(45);
    ylabel('Reaction times(ms)');
    title('RTs corrective sac');
    box off
    ESN_Beautify_Plot(fig, [10, 4]);
    %print(fig, fullfile(img_save_path, 'reactionTime-tag4.pdf'), '-dpdf', '-bestfit');
    
    %-----------------   End point error ---------------
    
    data_enderror = tag4_end_point_error{jj,1};
    all_data = [];
    group_labels = {};
    
    for i = 1:length(data_enderror)
        this_data = sqrt(data_enderror{i}(:,1).^2 + data_enderror{i}(:,2).^2);  
        all_data = [all_data; this_data];
        group_labels = [group_labels; repmat(labels(i), length(this_data), 1)];
    end
    
    group_labels = categorical(group_labels, labels, 'Ordinal', true);
    
    fig = figure;
    v = violinplot(all_data, group_labels);
    hold on
    cats = categories(group_labels);
    for i = 1:length(cats)
        this_group_data = all_data(group_labels == cats{i});
        med = median(this_group_data);
        plot([i - 0.2, i + 0.2], [med, med], 'k-', 'LineWidth', 1.5);
    end
    hold off
    
    for i = 1:length(v)
        v(i).ViolinColor = {Line_Color(i, :)};
    end
    xtickangle(45);
    ylabel('End point error (deg)');
    title('End point error corrective sac');
    box off
    ESN_Beautify_Plot(fig, [10, 4]);
    %print(fig, fullfile(img_save_path, 'enderror-tag4.pdf'), '-dpdf', '-bestfit');
end
end