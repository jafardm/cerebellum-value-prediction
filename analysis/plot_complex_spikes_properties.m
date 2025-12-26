function plot_complex_spikes_properties(data_path)

JDM_params_funcs
data_cs  = load(fullfile(data_path, 'population_data','CS_population_clique_prim_sac.mat'),...
    'data');

data_cs  = data_cs.data;

sess_CS    = cell2mat(cellfun(@(x) str2double(x(1:6)),...
    data_cs.cell_ids_tot,'UniformOutput',false));
cliques_CS = sess_CS.*10+data_cs.cliques;

num_CS = numel(cliques_CS);

% plot categorized CSes
load(fullfile(data_path, 'population_data', 'cs_on_rate_spe_rpe_prim_sac'),...
    'cs_on_data','cs_rate','cs_cell_ids');

cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg,cs_on_data,'UniformOutput',false));
cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg,cs_on_data,'UniformOutput',false));

cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg,cs_on_data,'UniformOutput',false));
cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg,cs_on_data,'UniformOutput',false));

num_cs = numel(cs_on_data);

ind_m = cs_on_rho_sac > cs_on_rho_vis;
ind_v = ~ind_m;

cs_on_rho = max(cs_on_rho_vis,cs_on_rho_sac);
cs_on_ang = cs_on_ang_vis/180*pi;
cs_on_ang(ind_m) = cs_on_ang_sac(ind_m)/180*pi;

ang_bins = params.sac.ang_values/180*pi;

[~,cs_on_bin] = min(abs(angdiff(repmat(ang_bins,num_cs,1),cs_on_ang*ones(size(ang_bins)))),[],2);

order_ang = nan(num_cs,8);

for counter_cs = 1:num_cs
    order_ang(counter_cs,:) = circshift(1:8,5-cs_on_bin(counter_cs));
end

cs_on_rate_vis = cell2mat(cellfun(@(x) x.vis.fr_avg,cs_on_data,'UniformOutput',false));
cs_on_rate_sac = cell2mat(cellfun(@(x) x.sac.fr_avg,cs_on_data,'UniformOutput',false));

% organizing the tuning of the cs rates
cs_rate_temp = mean(cs_rate,4,'omitnan');
cs_rate_rot        = nan(num_cs,500,8,3);
order_ang = [order_ang,order_ang(:,1)];
for counter_cs = 1:num_cs
    cs_rate_rot(counter_cs,:,:,:)    = cs_rate_temp(counter_cs,:,order_ang(counter_cs,1:8),:);
end

fig = figure;

n_col = 4;
n_row = 2;


subplot(n_row,n_col,4);
scatter(cs_on_rho_vis,cs_on_rho_sac,5,'k','filled');
hold on
plot([0,.6],[0,.6],'--k');
axis equal
xlabel('visual');
ylabel('motor');
xlim([0,.6]);
ylim([0,.6]);
title(['n = ' num2str(num_cs)])

sess_CS    = cell2mat(cellfun(@(x) str2double(x(1:6)),data_cs.cell_ids_tot,...
    'UniformOutput',false));
cliques_CS = sess_CS.*10+data_cs.cliques;

[cs_cliques_,~,~] = unique(cliques_CS);
num_cliques = numel(cs_cliques_);

ind_cl_cs = ismember(cs_cell_ids,data_cs.cell_ids_tot);

feat_space = [cs_on_rho_vis(ind_cl_cs),cs_on_rho_sac(ind_cl_cs)];

rho_dist = squareform(pdist(feat_space));

num_pc = numel(cliques_CS);

cs_on_ang_cl = cs_on_ang(ind_cl_cs);

ang_dist = nan(num_pc,num_pc);
for counter_=1:num_pc
    ang_dist(counter_,:) = abs(angdiff(cs_on_ang_cl,cs_on_ang_cl(counter_)*ones(num_pc,1)));
end

wit_ind = zeros(num_pc,num_pc);
for counter_clique = 1:num_cliques
    current_clique = cs_cliques_(counter_clique);
    clique_ind = cliques_CS == current_clique;
    wit_ind(clique_ind,clique_ind) = 1;
end
bet_ind = not(wit_ind);
wit_ind = wit_ind & tril(ones(num_pc,'logical'),-1);
bet_ind = bet_ind & tril(ones(num_pc,'logical'),-1);

dist_wit = rho_dist(wit_ind);
dist_bet = rho_dist(bet_ind);

ang_dist_wit = ang_dist(wit_ind)/pi*180;
ang_dist_bet = ang_dist(bet_ind)/pi*180;

num_wit = sum(wit_ind(:));
num_bet = sum(bet_ind(:));

subplot(n_row,n_col+2,n_col+2+1)
violinplot([dist_wit;dist_bet],...
    categorical([repelem({'wit'},num_wit,1);repelem({'bet'},num_bet,1)]));
title(['rank sum p = ', num2str(ranksum(dist_wit,dist_bet),'%.2e')])
ylabel('vis-mot dist (a.u.)')

subplot(n_row,n_col+2,n_col+2+2)
violinplot([ang_dist_wit;ang_dist_bet],...
    categorical([repelem({'wit'},num_wit,1);repelem({'bet'},num_bet,1)]));
title(['rank sum p = ', num2str(ranksum(ang_dist_wit,ang_dist_bet),'%.2e')])
yticks(0:45:180)
ylabel('angle diff (deg)')

[~,sorted_ind] = sort(cs_on_rho_sac-cs_on_rho_vis);
ind_sample1 = sorted_ind(1);
ind_sample2 = sorted_ind(end);

subplot(n_row,n_col,4);
hold on;
scatter(cs_on_rho_vis(ind_sample1),cs_on_rho_sac(ind_sample1),'b');
scatter(cs_on_rho_vis(ind_sample2),cs_on_rho_sac(ind_sample2),'r');

ang_values   = params.sac.ang_values;

subplot(2*n_row,n_col,3)
fr_amplitude = cs_on_rate_sac(ind_sample1,:);
ang_avg = cs_on_ang_sac(ind_sample1);
rho_avg = cs_on_rho_sac(ind_sample1);
plot_data_amp_mean = [fr_amplitude, fr_amplitude(1), nan]';
plot_data_deg_mean = [ang_values, ang_values(1), nan]';
polarplot(deg2rad(plot_data_deg_mean),plot_data_amp_mean, 'color', 'b');
hold on;
polarplot([0 deg2rad(ang_avg)],[0 rho_avg], 'color', 'b');

fr_amplitude = cs_on_rate_sac(ind_sample2,:);
ang_avg = cs_on_ang_sac(ind_sample2);
rho_avg = cs_on_rho_sac(ind_sample2);
plot_data_amp_mean = [fr_amplitude, fr_amplitude(1), nan]';
plot_data_deg_mean = [ang_values, ang_values(1), nan]';
polarplot(deg2rad(plot_data_deg_mean),plot_data_amp_mean, 'color', 'r');
polarplot([0 deg2rad(ang_avg)],[0 rho_avg], 'color', 'r');
% rlim([0,2.6]);
thetaticks(0:45:360);
title('motor')

subplot(2*n_row,n_col,n_col+3)
fr_amplitude = cs_on_rate_vis(ind_sample1,:);
ang_avg = cs_on_ang_vis(ind_sample1);
rho_avg = cs_on_rho_vis(ind_sample1);
plot_data_amp_mean = [fr_amplitude, fr_amplitude(1), nan]';
plot_data_deg_mean = [ang_values, ang_values(1), nan]';
polarplot(deg2rad(plot_data_deg_mean),plot_data_amp_mean, 'color', 'b');
hold on;
polarplot([0 deg2rad(ang_avg)],[0 rho_avg], 'color', 'b');

fr_amplitude = cs_on_rate_vis(ind_sample2,:);
ang_avg = cs_on_ang_vis(ind_sample2);
rho_avg = cs_on_rho_vis(ind_sample2);
plot_data_amp_mean = [fr_amplitude, fr_amplitude(1), nan]';
plot_data_deg_mean = [ang_values, ang_values(1), nan]';
polarplot(deg2rad(plot_data_deg_mean),plot_data_amp_mean, 'color', 'r');
polarplot([0 deg2rad(ang_avg)],[0 rho_avg], 'color', 'r');
% rlim([0,2.6]);
thetaticks(0:45:360);
title('visual')

time_ind = (-150:100)+250;
subplot(n_row,n_col,1)
cs_rate_ = squeeze(cs_rate_rot(sorted_ind,time_ind,5,2));
imagesc(time_ind-250,1:num_cs,cs_rate_);
colormap(params.rwb_map)
caxis([0,4]) 
hold on;
plot([-70,30],[895,895],'k','LineWidth',2)
title('Motor response')
ylabel('# complex spike')
xlabel('time from saccade onset (ms)')
xline(0,'--w','LineWidth',1.5)
s = colorbar;
ylabel(s,'rate (Hz)')

time_ind = (-100:240) + 250;
subplot(n_row,n_col,2)
cs_rate_ = squeeze(cs_rate_rot(sorted_ind,time_ind,5,1));
imagesc(time_ind-250,1:num_cs,cs_rate_);
colormap(params.rwb_map)
caxis([0,4]) 
hold on;
plot([40,85],[895,895],'k','LineWidth',2)
title('Visual response')
ylabel('# complex spike')
xlabel('time from visual onset (ms)')
xline(0,'--w','LineWidth',1.5)
s = colorbar;
ylabel(s,'rate (Hz)')

% similarity as a function of tuning
sim = [];
cs_ang_ = [];
cs_rho_ = [];
cs_rho_diff_norm_ = [];
cs_rho_diff_ = [];
cs_ang_diff_ = [];
clique_cs_on_ang = nan(num_cliques,1);
clique_cs_on_rho = nan(num_cliques,1);
for counter_clique = 1:num_cliques
    current_clique = cs_cliques_(counter_clique);
    current_ind = ismember(cliques_CS,current_clique);
    current_rhos = cs_on_rho_sac(current_ind);
    current_angs = cs_on_ang_sac(current_ind)*pi/180;
    current_vecs = current_rhos.*exp(1i*current_angs);
    vec_net = mean(current_vecs);
    clique_cs_on_ang(counter_clique) = angle(vec_net);
    clique_cs_on_rho(counter_clique) = abs(vec_net);
    sim = [sim;abs(current_vecs-vec_net)/abs(vec_net)];
    cs_rho_ = [cs_rho_;current_rhos];
    cs_ang_ = [cs_ang_;current_angs];
    cs_rho_diff_ = [cs_rho_diff_; (current_rhos-abs(vec_net))];
    cs_rho_diff_norm_ = [cs_rho_diff_norm_;(current_rhos-abs(vec_net))/abs(vec_net)];
    cs_ang_diff_ = [cs_ang_diff_;angdiff(current_angs,angle(vec_net)*ones(size(current_angs)))];
end

subplot(n_row,n_col+2,n_col+2+3)
histogram(cs_ang_*180/pi,50);
xticks(0:90:360)
xlim([0,360])
title('CS ang')

subplot(n_row,n_col+2,n_col+2+4)
histogram(cs_ang_diff_*180/pi,50);
xticks(-180:90:180)
xlim([-180,180])
title('CS ang diff from center')

subplot(n_row,n_col+2,n_col+2+5)
histogram(cs_rho_,30);
title('CS \rho')

subplot(n_row,n_col+2,n_col+2+6)
histogram(cs_rho_diff_,30);
title('CS \rho diff from center');

ESN_Beautify_Plot(fig,[20,8]);
print(fig, fullfile(data_path,'population_figs\ephys',...
    'complex_spikes_properties.pdf'), '-dpdf', '-bestfit');

end