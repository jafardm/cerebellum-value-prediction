function plot_neural_properties(data_path)

JDM_params_funcs

load(fullfile(data_path,'population_data','neural_prop_clique'),'acorr_tot_cs','acorr_tot_ss','cell_types',...
    'wave_tot_cs','wave_tot_ss','xcorr_tot','rate_bl_cs','rate_bl_ss');

ss_ind = ismember(cell_types,{'PC','SS'});
cs_ind = ismember(cell_types,{'PC','CS'});
ml_ind = ismember(cell_types,{'MLI'});
ml2_ind = ismember(cell_types,{'MLI2'});
pc_ind = ismember(cell_types,{'PC'});

fig = figure;

n_row = 5;
n_col = 5;

% auto-corr
subplot(n_row,n_col,1)
sig = acorr_tot_ss(ss_ind,:);
plot(-50:50,sig','LineWidth',.5,'Color',[0,0,0,.1]);
hold on
boundedline(-50:50,mean(sig),2*std(sig)/sqrt(sum(ss_ind)),...
    'LineWidth',1.5,'r','alpha');
xline(0,'--');
title(['SS (n=',num2str(sum(ss_ind)),')']);

subplot(n_row,n_col,2)
sig = acorr_tot_ss(ml_ind,:);
plot(-50:50,sig','LineWidth',.5,'Color',[0,0,0,.1]);
hold on
boundedline(-50:50,mean(sig),2*std(sig)/sqrt(sum(ml_ind)),...
    'LineWidth',1.5,'r','alpha');
xline(0,'--');
title(['MLI (n=',num2str(sum(ml_ind)),')']);


subplot(n_row,n_col,3)
sig = acorr_tot_ss(ml2_ind,:);
plot(-50:50,sig','LineWidth',.5,'Color',[0,0,0,.1]);
hold on
boundedline(-50:50,mean(sig),2*std(sig)/sqrt(sum(ml2_ind)),...
    'LineWidth',1.5,'r','alpha');
xline(0,'--');
title(['MLI2 (n=',num2str(sum(ml2_ind)),')']);

subplot(n_row,n_col,4)
sig = acorr_tot_cs(cs_ind,:);
plot(-500:5:500,sig','LineWidth',.5,'Color',[0,0,0,.1]);
hold on
boundedline(-500:5:500,mean(sig),2*std(sig)/sqrt(sum(cs_ind)),...
    'LineWidth',1.5,'r','alpha');
xline(0,'--');
title(['CS (n=',num2str(sum(cs_ind)),')']);

subplot(n_row,n_col,5)
sig = xcorr_tot(pc_ind,:);
plot(-50:50,sig','LineWidth',.5,'Color',[0,0,0,.1]);
hold on
boundedline(-50:50,mean(sig),2*std(sig)/sqrt(sum(pc_ind)),...
    'LineWidth',1.5,'r','alpha');
xline(0,'--');
title(['SS|CS (n=',num2str(sum(pc_ind)),')']);


% refractory violation
subplot(n_row,n_col,6)
sig = acorr_tot_ss(ss_ind,52);
histogram((sig));
ylabel('# neurons');
title([num2str(median(sig),2) '±' num2str(mad(sig,1),2)])
xlabel('noise rate (Hz)');

subplot(n_row,n_col,7)
sig = acorr_tot_ss(ml_ind,52);
histogram((sig));
title([num2str(median(sig),2) '±' num2str(mad(sig,1),2)])
ylabel('# neurons');
xlabel('noise rate (Hz)');

subplot(n_row,n_col,8)
sig = acorr_tot_ss(ml2_ind,52);
histogram((sig));
title([num2str(median(sig),2) '±' num2str(mad(sig,1),2)])
ylabel('# neurons');
xlabel('noise rate (Hz)');

subplot(n_row,n_col,9)
sig = acorr_tot_cs(cs_ind,102);
histogram((sig));
title([num2str(median(sig),2) '±' num2str(mad(sig,1),2)])
ylabel('# neurons');
xlabel('noise rate (Hz)');

subplot(n_row,n_col,10)
sig = sum(xcorr_tot(pc_ind,51+(1:5)),2);
histogram((sig));
title([num2str(median(sig),2) '±' num2str(mad(sig,1),2)])
ylabel('# neurons');
xlabel('noise rate (Hz)');

% rate dist
subplot(n_row,n_col,11)
ind_selected = ml_ind|ml2_ind|ss_ind;
sig = rate_bl_ss(ind_selected);
violinplot(sig,cell_types(ind_selected));
title('rate dist.')
ylabel('rate (Hz)')

subplot(n_row,n_col,12)
violinplot(rate_bl_cs(cs_ind),cell_types(cs_ind));
ylabel('rate (Hz)')
title('CS rate dist.')
yline(1,'--');

subplot(n_row,n_col,13)
red_dim_wave_ss = run_umap(wave_tot_ss(ind_selected,:),...
    'randomize','false','verbose','none');
gscatter(red_dim_wave_ss(:,1),red_dim_wave_ss(:,2),cell_types(ind_selected));
legend('Box','off')
dn_ind_ss_ = red_dim_wave_ss(:,1) <= 0;
ax_ind_ss_ = ~dn_ind_ss_;

dn_ind_ss = zeros(numel(ind_selected),1,'logical');
ax_ind_ss = zeros(numel(ind_selected),1,'logical');

dn_ind_ss(ind_selected) = dn_ind_ss_;
ax_ind_ss(ind_selected) = ax_ind_ss_;

subplot(n_row,n_col,14)
red_dim_wave_cs = run_umap(wave_tot_cs(cs_ind,:),...
    'randomize','false','verbose','none');
gscatter(red_dim_wave_cs(:,1),red_dim_wave_cs(:,2),cell_types(cs_ind));
legend('Box','off')
rng default
clust_cs = kmeans(red_dim_wave_cs,2,'Replicates',1000);

ax_ind_cs_ = clust_cs == 1;
dn_ind_cs_ = ~ax_ind_cs_;

dn_ind_cs = zeros(numel(cs_ind),1,'logical');
ax_ind_cs = zeros(numel(cs_ind),1,'logical');

dn_ind_cs(cs_ind) = dn_ind_cs_;
ax_ind_cs(cs_ind) = ax_ind_cs_;

% waveform
% ss
len_x = 50;
len_y = 150;
plt_y = [-100,100];
span_x = ((1:len_x)-len_x/3)/30;
span_y = ((1:len_y)-len_y/2)*5;
subplot(n_row,n_col,16)
sig = reshape(wave_tot_ss(ss_ind&ax_ind_ss,:),[],len_y,len_x);

% find size
sig_peak = -squeeze(sig(:,:,17));
num_cell = size(sig_peak,1);

w_r = nan(num_cell,1);
w_l = nan(num_cell,1);
for counter = 1:num_cell
    idx_r = find(sig_peak(counter,77:end)<sig_peak(counter,76)/2,1);
    w_r(counter) = interp1(sig_peak(counter,76+(idx_r-1:idx_r)),...
        [idx_r-1,idx_r]*5,...
        sig_peak(counter,76)/2,'linear');
    
    idx_l = find(sig_peak(counter,75:-1:1)<sig_peak(counter,76)/2,1);
    if ~isempty(idx_l) && idx_l > 1
        x_vals = sig_peak(counter, 76 - (idx_l-1:idx_l));
        t_vals = [idx_l-1, idx_l]*5;
        w_l(counter) = interp1(x_vals, t_vals, sig_peak(counter,76)/2, 'linear');
   end

end

sz_ss_ax = max(w_r,w_l);

plt_ind = span_y <= plt_y(2) & span_y >= plt_y(1);
sig__ = sig(:,plt_ind,:);
img_axes = []; 
ax = gca;
imagesc(span_x,span_y(plt_ind),squeeze(mean(sig__)));
caxis([-1,1])
colormap(params.rwb_map);
img_axes(end+1) = ax;
hold on;
plot([0,0],mean(sz_ss_ax)*[-1,1],'LineWidth',2,'Color','y')
yyaxis right
sig_ = squeeze(sig(:,end/2,:));
boundedline(span_x,mean(sig_),2*std(sig_)/sqrt(sum(ss_ind&ax_ind_ss)),'k','alpha');
ylim([-1,1])
title(['SS (n=',num2str(sum(ss_ind&ax_ind_ss)),')']);

subplot(n_row,n_col,21)
sig = reshape(wave_tot_ss(ss_ind&dn_ind_ss,:),[],len_y,len_x);

% find size
sig_peak = squeeze(sig(:,:,17));
num_cell = size(sig_peak,1);

w_r = nan(num_cell,1);
w_l = nan(num_cell,1);
for counter = 1:num_cell
    idx_r = find(sig_peak(counter,77:end)<sig_peak(counter,76)/2,1);
   if ~isempty(idx_r) && idx_r > 1
    x_vals = sig_peak(counter, 76 + (idx_r-1:idx_r));
    t_vals = [idx_r-1, idx_r]*5;
    w_r(counter) = interp1(x_vals, t_vals, sig_peak(counter,76)/2, 'linear');
   end

    if ~isempty(idx_l) && idx_l > 1
        x_vals = sig_peak(counter, 76 - (idx_l-1:idx_l));
        t_vals = [idx_l-1, idx_l]*5;
        w_l(counter) = interp1(x_vals, t_vals, sig_peak(counter,76)/2, 'linear');
    end

end

sz_ss_dn = max(w_r,w_l);

sig__ = sig(:,plt_ind,:);
ax = gca;
imagesc(span_x,span_y(plt_ind),squeeze(mean(sig__)));
caxis([-1,1])
colormap(params.rwb_map);
img_axes(end+1) = ax;

yyaxis right
sig_ = squeeze(sig(:,end/2,:));
boundedline(((1:len_x)-len_x/3)/30,mean(sig_),2*std(sig_)/sqrt(sum(ss_ind&dn_ind_ss)),'k','alpha');
ylim([-1,1])
title(['SS (n=',num2str(sum(ss_ind&dn_ind_ss)),')']);

% mli
subplot(n_row,n_col,17)
sig = reshape(wave_tot_ss(ml_ind&ax_ind_ss,:),[],len_y,len_x);

% find size
sig_peak = -squeeze(sig(:,:,17));
num_cell = size(sig_peak,1);

w_r = nan(num_cell,1);
w_l = nan(num_cell,1);
for counter = 1:num_cell
    idx_r = find(sig_peak(counter,77:end)<sig_peak(counter,76)/2,1);
   if ~isempty(idx_r) && idx_r > 1
    x_vals = sig_peak(counter, 76 + (idx_r-1:idx_r));
    t_vals = [idx_r-1, idx_r]*5;
    w_r(counter) = interp1(x_vals, t_vals, sig_peak(counter,76)/2, 'linear');
  end

    idx_l = find(sig_peak(counter,75:-1:1)<sig_peak(counter,76)/2,1);
    if ~isempty(idx_l) && idx_l > 1
        x_vals = sig_peak(counter, 76 - (idx_l-1:idx_l));
        t_vals = [idx_l-1, idx_l]*5;
        w_l(counter) = interp1(x_vals, t_vals, sig_peak(counter,76)/2, 'linear');
   end

end

sz_ml_ax = max(w_r,w_l);

sig__ = sig(:,plt_ind,:);
ax = gca;
imagesc(span_x,span_y(plt_ind),squeeze(mean(sig__)));
caxis([-1,1])
colormap(params.rwb_map);
img_axes(end+1) = ax;

hold on;
plot([0,0],mean(sz_ml_ax)*[-1,1],'LineWidth',2,'Color','y')
yyaxis right
sig_ = squeeze(sig(:,end/2,:));
boundedline(((1:len_x)-len_x/3)/30,mean(sig_),2*std(sig_)/sqrt(sum(ml_ind&ax_ind_ss)),'k','alpha');
ylim([-1,1])
title(['MLI (n=',num2str(sum(ml_ind&ax_ind_ss)),')']);

subplot(n_row,n_col,22)
sig = reshape(wave_tot_ss(ml_ind&dn_ind_ss,:),[],len_y,len_x);

% find size
sig_peak = squeeze(sig(:,:,17));
num_cell = size(sig_peak,1);

w_r = nan(num_cell,1);
w_l = nan(num_cell,1);
for counter = 1:num_cell
    idx_r = find(sig_peak(counter,77:end)<sig_peak(counter,76)/2,1);
    if ~isempty(idx_r) && idx_r > 1
        x_vals = sig_peak(counter, 76 + (idx_r-1:idx_r));
        t_vals = [idx_r-1, idx_r]*5;
        w_r(counter) = interp1(x_vals, t_vals, sig_peak(counter,76)/2, 'linear');
   end
   
    idx_l = find(sig_peak(counter,75:-1:1)<sig_peak(counter,76)/2,1);
  if ~isempty(idx_l) && idx_l > 1
    x_vals = sig_peak(counter, 76 - (idx_l-1:idx_l));
    t_vals = [idx_l-1, idx_l]*5;
    w_l(counter) = interp1(x_vals, t_vals, sig_peak(counter,76)/2, 'linear');
  end

end

sz_ml_dn = max(w_r,w_l);

sig__ = sig(:,plt_ind,:);
ax = gca;
imagesc(span_x,span_y(plt_ind),squeeze(mean(sig__)));
caxis([-1,1])
colormap(params.rwb_map);
img_axes(end+1) = ax;

yyaxis right
sig_ = squeeze(sig(:,end/2,:));
boundedline(((1:len_x)-len_x/3)/30,mean(sig_),2*std(sig_)/sqrt(sum(ml_ind&dn_ind_ss)),'k','alpha');
ylim([-1,1])
title(['MLI (n=',num2str(sum(ml_ind&dn_ind_ss)),')']);

% mli2
subplot(n_row,n_col,18)
sig = reshape(wave_tot_ss(ml2_ind&ax_ind_ss,:),[],len_y,len_x);

% find size
sig_peak = -squeeze(sig(:,:,17));
num_cell = size(sig_peak,1);

w_r = nan(num_cell,1);
w_l = nan(num_cell,1);
for counter = 1:num_cell
    idx_r = find(sig_peak(counter,77:end)<sig_peak(counter,76)/2,1);
    if ~isempty(idx_r) && idx_r > 1
        x_vals = sig_peak(counter, 76 + (idx_r-1:idx_r));
        t_vals = [idx_r-1, idx_r]*5;
        w_r(counter) = interp1(x_vals, t_vals, sig_peak(counter,76)/2, 'linear');
   end
   
    idx_l = find(sig_peak(counter,75:-1:1)<sig_peak(counter,76)/2,1);
  if ~isempty(idx_l) && idx_l > 1
    x_vals = sig_peak(counter, 76 - (idx_l-1:idx_l));
    t_vals = [idx_l-1, idx_l]*5;
    w_l(counter) = interp1(x_vals, t_vals, sig_peak(counter,76)/2, 'linear');
  end

end

sz_ml2_ax = max(w_r,w_l);

sig__ = sig(:,plt_ind,:);
ax = gca;
imagesc(span_x,span_y(plt_ind),squeeze(mean(sig__)));
caxis([-1,1])
colormap(params.rwb_map);
img_axes(end+1) = ax;

hold on;
plot([0,0],mean(sz_ml2_ax)*[-1,1],'LineWidth',2,'Color','y')
yyaxis right
sig_ = squeeze(sig(:,end/2,:));
boundedline(((1:len_x)-len_x/3)/30,mean(sig_),2*std(sig_)/sqrt(sum(ml2_ind&ax_ind_ss)),'k','alpha');
ylim([-1,1])
title(['MLI2 (n=',num2str(sum(ml2_ind&ax_ind_ss)),')']);

subplot(n_row,n_col,23)
sig = reshape(wave_tot_ss(ml2_ind&dn_ind_ss,:),[],len_y,len_x);

% find size
sig_peak = squeeze(sig(:,:,17));
num_cell = size(sig_peak,1);

w_r = nan(num_cell,1);
w_l = nan(num_cell,1);
for counter = 1:num_cell
    idx_r = find(sig_peak(counter,77:end)<sig_peak(counter,76)/2,1);
    if ~isempty(idx_r) && idx_r > 1
        x_vals = sig_peak(counter, 76 + (idx_r-1:idx_r));
        t_vals = [idx_r-1, idx_r]*5;
        w_r(counter) = interp1(x_vals, t_vals, sig_peak(counter,76)/2, 'linear');
   end
   
    idx_l = find(sig_peak(counter,75:-1:1)<sig_peak(counter,76)/2,1);
  if ~isempty(idx_l) && idx_l > 1
    x_vals = sig_peak(counter, 76 - (idx_l-1:idx_l));
    t_vals = [idx_l-1, idx_l]*5;
    w_l(counter) = interp1(x_vals, t_vals, sig_peak(counter,76)/2, 'linear');
  end

end

sz_ml2_dn = max(w_r,w_l);

sig__ = sig(:,plt_ind,:);
ax = gca;
imagesc(span_x,span_y(plt_ind),squeeze(mean(sig__)));
caxis([-1,1])
colormap(params.rwb_map);
img_axes(end+1) = ax;
yyaxis right
sig_ = squeeze(sig(:,end/2,:));
boundedline(((1:len_x)-len_x/3)/30,mean(sig_),2*std(sig_)/sqrt(sum(ml2_ind&dn_ind_ss)),'k','alpha');
ylim([-1,1])
title(['MLI2 (n=',num2str(sum(ml2_ind&dn_ind_ss)),')']);

% cs
subplot(n_row,n_col,19)
sig = reshape(wave_tot_cs(cs_ind&ax_ind_cs,:),[],len_y,len_x);
sig__ = sig(:,plt_ind,:);
ax = gca;
imagesc(span_x,span_y(plt_ind),squeeze(mean(sig__)));
caxis([-1,1])
colormap(params.rwb_map);
img_axes(end+1) = ax;
yyaxis right
sig_ = squeeze(sig(:,end/2,:));
boundedline(((1:len_x)-len_x/3)/30,mean(sig_),2*std(sig_)/sqrt(sum(cs_ind&ax_ind_cs)),'k','alpha');
ylim([-1,1])
title(['CS (n=',num2str(sum(cs_ind&ax_ind_cs)),')']);
subplot(n_row,n_col,24)
sig = reshape(wave_tot_cs(cs_ind&dn_ind_cs,:),[],len_y,len_x);
sig__ = sig(:,plt_ind,:);
imagesc(span_x,span_y(plt_ind),squeeze(mean(sig__)));
caxis([-1,1])
colormap(params.rwb_map);
yyaxis right
sig_ = squeeze(sig(:,end/2,:));
boundedline(((1:len_x)-len_x/3)/30,mean(sig_),2*std(sig_)/sqrt(sum(cs_ind&dn_ind_cs)),'k','alpha');
ylim([-1,1])
title(['CS (n=',num2str(sum(cs_ind&dn_ind_cs)),')']);
cb = colorbar('Position', [0.92 0.15 0.02 0.3]);
colormap(cb, params.rwb_map);
set(img_axes, 'CLim', [-1 1]);

% subplot(n_row,n_col,15);
% violinplot(2*[sz_ss_ax;sz_ml_ax;sz_ml2_ax],...
%     repelem({'SS','MLI','MLI2'},...
%     [length(sz_ss_ax),length(sz_ml_ax),length(sz_ml2_ax)]));
% ylabel('soma size (um)')
% title('FWHM')
% 
% subplot(n_row,n_col,20);
% title(['SS:',num2str(median(sz_ss_ax)),'±', num2str(mad(sz_ss_ax,1)),...
%     'MLI:',num2str(median(sz_ml_ax)),'±', num2str(mad(sz_ml_ax,1)),...
%     'MLI2:',num2str(median(sz_ml2_ax)),'±', num2str(mad(sz_ml2_ax,1))])

ESN_Beautify_Plot(fig,[16,9])
img_save_path = fullfile('C:\Users\Jafar\Documents\reward','population_figs\ephys');
print(fig,fullfile(img_save_path,'neural_properties_clique.pdf'),'-dpdf','-bestfit')

end