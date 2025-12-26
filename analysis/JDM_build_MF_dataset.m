function JDM_build_MF_dataset

JDM_params_funcs
load(fullfile(data_path,'population_data','behave_data'));
user = 'JDM'; loc = 'ctx';

path_data_monkey_sorted = params.path_data_monkey_sorted;
animal_list = params.animal_list;
session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);

num_animal = numel(session_list);

% last dim 1:tgt, 2:spont
data_out.MF_rate_tot_sac = zeros(1,500,8,2);
data_out.MF_rate_tot_vis = zeros(1,500,8);
% eye traces
data_out.MF_vm_tot_sac   = zeros(1,500,8,2);

data_out.MF_rate_tot_fix = cell(1);
data_out.MF_cell_ids_tot = cell(1);
data_out.MF_ss_on_data   = cell(1);

for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess = numel(current_sess_list);

    current_sess_cr_ma = data_behav.sess_cr_ma{counter_animal};
    current_sess_cr_ms = data_behav.sess_cr_ms{counter_animal};

    for counter_sess = 1:num_sess
        current_sess = current_sess_list{counter_sess};
        disp(['     ' current_sess]);

        current_cr_ma = current_sess_cr_ma{counter_sess};
        current_cr_ms = current_sess_cr_ms{counter_sess};
        sac_crit.ma   = current_cr_ma;
        sac_crit.ms   = current_cr_ms;

        % load eye_traces
        data_eye = MAF_load_eye_traces(current_path,current_sess);

        % load units info
        units_info = MAF_extract_cell_metadata(current_path,current_sess);

        cell_ids     = units_info.cell_list;
        cell_types   = units_info.cell_type;

        cells_ind  = ismember(cell_types,{'axonal'});

        % mossy fiber
        num_cell           = sum(cells_ind);
        current_cell_ids   = cell_ids(cells_ind);
        current_cell_types = cell_types(cells_ind);
        cell_rate_tot_sac_ = nan(num_cell,500,8,2);
        cell_rate_tot_vis_ = nan(num_cell,500,8);
        cell_vm_tot_sac_   = nan(num_cell,500,8,2);
        cell_rate_tot_fix_ = cell(num_cell,1);
        cell_ids_tot_      = cell(num_cell,1);
        cell_ss_on_data_   = cell(num_cell,1);

        parfor counter_cell = 1:num_cell
            current_cell = current_cell_ids{counter_cell};
            current_cell_type = current_cell_types{counter_cell};
            disp([num2str(counter_cell),'/',num2str(num_cell),': ',...
                current_cell_type, '    ', current_cell]);

            res = load_mf_data_tuned(current_path,...
                current_sess,current_cell,...
                sac_crit,data_eye);

            cell_rate_tot_sac_(counter_cell,:,:,:) = cat(4,res.tgt_rate,res.spnt_rate);
            cell_rate_tot_vis_(counter_cell,:,:,:) = res.vis_rate;

            cell_rate_tot_fix_{counter_cell}       = res.fix;
            cell_ids_tot_{counter_cell}            = current_cell;
            cell_ss_on_data_{counter_cell}         = res.ss_on_data;

            cell_vm_tot_sac_(counter_cell,:,:,:)   = cat(4,res.tgt_vm,res.spnt_vm);
        end

        data_out.MF_rate_tot_sac = cat(1,data_out.MF_rate_tot_sac,cell_rate_tot_sac_);
        data_out.MF_rate_tot_vis = cat(1,data_out.MF_rate_tot_vis,cell_rate_tot_vis_);
        
        data_out.MF_rate_tot_fix = cat(1,data_out.MF_rate_tot_fix,cell_rate_tot_fix_);
        data_out.MF_cell_ids_tot = cat(1,data_out.MF_cell_ids_tot,cell_ids_tot_);
        data_out.MF_ss_on_data   = cat(1,data_out.MF_ss_on_data,cell_ss_on_data_);

        data_out.MF_vm_tot_sac   = cat(1,data_out.MF_vm_tot_sac,cell_vm_tot_sac_);
    end
end

% remove the first empty element
data_out.MF_rate_tot_sac(1,:,:,:) = [];
data_out.MF_rate_tot_vis(1,:,:,:) = [];

data_out.MF_rate_tot_fix(1)       = [];
data_out.MF_ss_on_data(1)         = [];
data_out.MF_cell_ids_tot(1)       = [];

data_out.MF_vm_tot_sac(1,:,:,:)   = [];
out_path = fullfile(data_path,'MF_population_clique');
save(fullfile(out_path, sprintf('%s.mat',file_name)), data_out,'-v7.3');   

end

% MF
function res = load_mf_data_tuned(current_path,current_sess,...
    current_cell,sac_crit,data_eye)

MAF_params_funcs;

tag_id = 1;
time_span = 250;
alignment = 'vmax';

data = MAF_load_cell(current_path,current_sess,current_cell);
eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);

data_recs = [data.data_recordings];
neural_sac_data = MAF_combine_dataset_sac(data_recs,...
    [data.rec_info],{alignment},time_span,tag_id);
eye_vm = MAF_combine_vm_traces(eye_traces,{alignment},time_span,tag_id);

neur_data.eye_vm    = eye_vm.(alignment).eye_vm;
neur_data.SS        = neural_sac_data.(alignment).SS;
sac_data.sac        = neural_sac_data.(alignment).sac;
sac_data.fix        = neural_sac_data.(alignment).fix;
sac_data.sac.has_ss = neural_sac_data.(alignment).has_ss;
sac_data_dir        = MAF_bin_by_ang(neur_data,sac_data,alignment);

neural_vis_data = MAF_combine_dataset_sac(data_recs,...
    [data.rec_info],{'visual'},time_span,1);
neur_data_vis.SS    = neural_vis_data.visual.SS;
vis_data.sac        = neural_vis_data.visual.sac;
vis_data.fix        = neural_vis_data.visual.fix;
vis_data.sac.has_ss = neural_vis_data.visual.has_ss;
vis_data_dir        = MAF_bin_by_ang(neur_data_vis,vis_data,'visual');

res.ss_on_data = MAF_SS_on_analysis_fr([data_recs.eye],...
    [data_recs.Neural_Data],1:10);

on_ang = res.ss_on_data.ang_avg * pi/180;

ang_bins = params.sac.ang_values/180*pi;

[~,on_bin] = min(abs(angdiff(ang_bins,on_ang*ones(size(ang_bins)))));

order_ang = circshift(1:8,1-on_bin);

SS_rate_tot_sac_tgt  = zeros(500,8);
SS_rate_tot_sac_spnt = zeros(500,8);
SS_rate_tot_vis      = zeros(500,8);

eye_vm_tot_sac_tgt   = zeros(500,8);
eye_vm_tot_sac_spnt  = zeros(500,8);

parfor counter_dir = 1:8
    % shifting space wrt SS on
    has_ss = sac_data_dir(order_ang(counter_dir)).sac.has_ss;
    current_amp = sac_data_dir(order_ang(counter_dir)).sac.eye_amp(has_ss);
    current_vel = sac_data_dir(order_ang(counter_dir)).sac.eye_vm_max(has_ss);
    current_dect = (sac_data_dir(order_ang(counter_dir)).sac.time_offset(has_ss)-...
        sac_data_dir(order_ang(counter_dir)).sac.time_vmax(has_ss))*1e3;
    current_acct = (sac_data_dir(order_ang(counter_dir)).sac.time_vmax(has_ss)-...
        sac_data_dir(order_ang(counter_dir)).sac.time_onset(has_ss))*1e3;
    current_raster = sac_data_dir(order_ang(counter_dir)).SS(:,has_ss);
    current_tag    = sac_data_dir(order_ang(counter_dir)).sac.tag(has_ss);
    current_eye_vm = sac_data_dir(order_ang(counter_dir)).eye_vm(:,has_ss);

    % remove invalid sacs
    poly_x = sac_crit.ms(1,:);
    poly_y = sac_crit.ms(2,:);
    [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);

    acc_dec_ = current_dect./current_acct;
    poly_x = sac_crit.ma(1,:);
    poly_y = sac_crit.ma(2,:);
    [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);
    
    % negative acc, dec, outliers in main-seq and amp-acc/dec plot 
    rm_ind = ~(in_ms | on_ms) |...
        ~(in_acc | on_acc) | (current_acct<=0) | (current_dect<=0);

    current_raster(:,rm_ind) = [];
    current_eye_vm(:,rm_ind) = [];
    current_tag(:,rm_ind)    = [];

    ind_tag_tgt  = ismember(current_tag,[1,4,6,7]);
    ind_tag_spnt = ismember(current_tag,10);

    % tgt
    rate_tgt   = ESN_smooth(mean(current_raster(:,ind_tag_tgt),2))*1e3;
    eye_vm_tgt = mean(current_eye_vm(:,ind_tag_tgt),2);

    SS_rate_tot_sac_tgt(:,counter_dir) = rate_tgt;
    eye_vm_tot_sac_tgt(:,counter_dir)  = eye_vm_tgt;
    
    % spont
    rate_spnt   = ESN_smooth(mean(current_raster(:,ind_tag_spnt),2))*1e3;
    eye_vm_spnt = mean(current_eye_vm(:,ind_tag_spnt),2);

    SS_rate_tot_sac_spnt(:,counter_dir) = rate_spnt;
    eye_vm_tot_sac_spnt(:,counter_dir)  = eye_vm_spnt;

    % visual
    % shifting space wrt SS on
    has_ss = vis_data_dir(order_ang(counter_dir)).sac.has_ss;
    current_amp = vis_data_dir(order_ang(counter_dir)).sac.eye_amp(has_ss);
    current_vel = vis_data_dir(order_ang(counter_dir)).sac.eye_vm_max(has_ss);
    current_dect = (vis_data_dir(order_ang(counter_dir)).sac.time_offset(has_ss)-...
        vis_data_dir(order_ang(counter_dir)).sac.time_vmax(has_ss))*1e3;
    current_acct = (vis_data_dir(order_ang(counter_dir)).sac.time_vmax(has_ss)-...
        vis_data_dir(order_ang(counter_dir)).sac.time_onset(has_ss))*1e3;
    current_raster = vis_data_dir(order_ang(counter_dir)).SS(:,has_ss);
    current_tag    = vis_data_dir(order_ang(counter_dir)).sac.tag(has_ss);

    % remove invalid sacs
    poly_x = sac_crit.ms(1,:);
    poly_y = sac_crit.ms(2,:);
    [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);

    acc_dec_ = current_dect./current_acct;
    poly_x = sac_crit.ma(1,:);
    poly_y = sac_crit.ma(2,:);
    [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);
    
    % negative acc, dec, outliers in main-seq and amp-acc/dec plot 
    rm_ind = ~(in_ms | on_ms) |...
        ~(in_acc | on_acc) | (current_acct<=0) | (current_dect<=0);

    current_raster(:,rm_ind) = [];
    current_tag(:,rm_ind)    = [];

    SS_rate_tot_vis(:,counter_dir) = ESN_smooth(mean(current_raster,2))*1e3;

end

res.tgt_rate  = SS_rate_tot_sac_tgt;
res.spnt_rate = SS_rate_tot_sac_spnt;
res.vis_rate  = SS_rate_tot_vis;
res.tgt_vm    = eye_vm_tot_sac_tgt;
res.spnt_vm   = eye_vm_tot_sac_spnt;

% fix
data_fix = MAF_combine_dataset_sac(data.data_recordings,...
        data.rec_info,{'offset'},250,1:10);

fix_dur_after  = ...
    (data_fix.offset.fix.time_fix_offset -...
    data_fix.offset.sac.time_offset)*1e3;

ind_after  = ...
    (data_fix.offset.fix.fix_validity_after == 1) &...
    (fix_dur_after > 200);

raster_offset  = data_fix.offset.SS;

raster_after = raster_offset(:,ind_after);

window_after  = (100:150) + 250;

fr_fix_after  = mean(raster_after(window_after,:))*1e3;

x_fix_after  = data_fix.offset.fix.x_fix_after(ind_after);
y_fix_after  = data_fix.offset.fix.y_fix_after(ind_after);

x_proj = x_fix_after*cos(on_ang) + y_fix_after*sin(on_ang);

res.fix.mdl = fitlm(fr_fix_after,x_proj);
res.fix.x_after = x_fix_after;
res.fix.y_after = y_fix_after;
res.fix.rate_after = fr_fix_after;

edges = params.pop.sac.fix_edges;
n_bin = numel(edges) - 1;
rate_time = nan(500,n_bin);
bins = discretize(x_proj,edges);
parfor counter_bin = 1:n_bin
    rate_time(:,counter_bin) = ESN_smooth(mean(raster_after(:,bins == counter_bin),2))*1e3;
end

res.fix.rate_time = rate_time;
res.fix.window_after = window_after - 250;

end