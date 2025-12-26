function plot_prob_view_adj_matrix_exp(data_path)
JDM_params_funcs;

current_path = params.path_data_monkey_sorted{2};
current_sess = '2024-12-30';
subj_name    = '65F';
out_dir      = fullfile(data_path,'population_figs','ephys');
%% probe view and CS-on clique
fig = figure;
sgtitle([subj_name '-' 'recording session: ' current_sess])
hold on
subplot(1,5,1)
res_view = MAF_plot_sess_elec_view(current_path,current_sess,...
    {'CS','PC','SS','MLI','MLI2','axonal','GLI'},1);
axis off

cell_clique = res_view.units_info.cell_clique;
cs_on_tot   = res_view.cs_on_tot;
cell_type   = res_view.units_info.cell_type;
ind_cs      = find(ismember(cell_type,{'PC'}));
cs_clique = cell_clique(ind_cs);
cs_on_tot   = cs_on_tot(ind_cs);
cs_rate = cell2mat(cellfun(@(x) x.vis.fr_avg, cs_on_tot,'UniformOutput',false));
cs_on_ang = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_tot,'UniformOutput',false));
cs_on_rho = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_tot,'UniformOutput',false));
vec = cs_on_rho.*exp(1i*cs_on_ang*pi/180);
vec_clique = grpstats(vec,cs_clique);
[cs_rate_clique_m,cs_rate_clique_se] = grpstats(cs_rate,cs_clique,{'mean','sem'});
clique_ang = angle(vec_clique);
clique_rho = abs(vec_clique);

cliques_cs = unique(cs_clique);
num_clique = numel(cliques_cs);

colors_ = lines(num_clique+1);
colors_ = colors_(2:end,:);

for counter_clique = 1:num_clique
    subplot(5,5,5*(num_clique-counter_clique)+3)
    current_clique_rho  = clique_rho(counter_clique);
    current_clique_ang  = clique_ang(counter_clique);
    current_cs_rate_avg = cs_rate_clique_m(counter_clique,:);
    current_cs_rate_se  = cs_rate_clique_se(counter_clique,:);
    % Von Mises
    if current_clique_rho < 0.53
        vonM_k = 2*current_clique_rho + current_clique_rho + 5*current_clique_rho/6;
    elseif current_clique_rho>=0.53 && current_clique_rho<0.85
        vonM_k = -.4 + 1.39*current_clique_rho + 0.43/(1-current_clique_rho);
    else
        vonM_k = 1/(current_clique_rho^3 - 4*current_clique_rho^2 + 3*current_clique_rho);
    end
    % evaluate pdf
    vonM_var = 1 - (besseli(1,vonM_k) / besseli(0,vonM_k));
    vonM_std = wrapTo360(rad2deg(sqrt(vonM_var)));

    fr_amplitude = current_cs_rate_avg;
    vonMises_std = vonM_std;
    ang_avg = current_clique_ang*180/pi;
    rho_avg = current_clique_rho;
    std_curv_ang = (ang_avg-vonMises_std) : 2 : (ang_avg+vonMises_std);
    std_curv_amp = repmat(rho_avg, length(std_curv_ang), 1);

    ang_values   = params.sac.ang_values;

    plot_data_amp_mean = [fr_amplitude,fr_amplitude(1), nan]';
    plot_data_amp_se = [current_cs_rate_se,current_cs_rate_se(1), nan]';
    plot_data_deg_mean = [ang_values, ang_values(1), nan]';
    polarplot(deg2rad(plot_data_deg_mean),plot_data_amp_mean, 'color', colors_(counter_clique,:));
    hold on;
    polarplot(deg2rad(plot_data_deg_mean),plot_data_amp_mean+plot_data_amp_se, 'color', colors_(counter_clique,:));
    polarplot(deg2rad(std_curv_ang),std_curv_amp, 'color', colors_(counter_clique,:));
    polarplot([0 deg2rad(ang_avg)],[0 rho_avg], 'color', colors_(counter_clique,:));
    set(findall(gcf,'Type','polaraxes'),'RLim',[0 1.5]);
    thetaticks([]);

end
ESN_Beautify_Plot(fig,[6,8]);
print(fig,fullfile(out_dir,'probe_view'),'-dpdf','-bestfit')
%%
    res = MAF_load_sess_adj_matrix(current_path,current_sess,{'SS','PC','MLI','MLI2'},0,1,res_view.jitter_x);
    %save('sample_adj_mat.mat','res');
    
    fig = figure;
    sgtitle([subj_name '-' 'recording session: ' current_sess])
    rectangle('Position',[1,1,size(res.Adj)],'FaceColor',[.7,.7,.7],'EdgeColor',[.7,.7,.7],'LineWidth',0.1);
    hold on
    s = pcolor(blkdiag(res.Adj,0));
    s.EdgeColor = 'none';
    xticks(1:size(res.Adj,1))
    yticks(1:size(res.Adj,1))
    yticklabels({})
    xticklabels(res.cell_names)
    colormap([0.5,0.5,0.5;params.rwb_map]);
    caxis([-20,20]);
    view(90,-90)
    c = colorbar('Location','eastoutside');
    ylabel(c,'rate (Hz)')
    axis equal
    axis tight
    hold on
    rec_size = diff(res.Adj_xline);
    rec_size(end) = rec_size(end);
    all_labels = unique(res.cell_cliques);
    current_labels = ...
        unique(res.cell_cliques);
    ind_ = ismember(all_labels,current_labels);
    color_ = params.clust_color(ind_,:);
    for counter_bound = 1:numel(current_labels)
        st = res.Adj_xline(counter_bound)+.5;
        sz = rec_size(counter_bound);
        pos_ = [st,st,sz,sz];
        rectangle('Position',pos_,...
            'EdgeColor',color_(counter_bound,:),...
            'LineWidth',2);
    end
    
    ESN_Beautify_Plot(fig,[10,10]);
    print(fig,fullfile(out_dir,'adj_matrix_exp'),'-dpdf','-bestfit')
end
