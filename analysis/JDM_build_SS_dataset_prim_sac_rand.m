function JDM_build_SS_dataset_prim_sac_rand(data_path)

JDM_params_funcs
load(fullfile(data_path,'population_data','behave_data'));

user = 'JDM';
loc  = 'ctx';

path_data_monkey_sorted = params.path_data_monkey_sorted(1:2);
animal_list = params.animal_list(1:2);
session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);

num_animal = numel(session_list);

cell_var = {'SS','MLI'};
cell_var_types = {{'SS','PC'},'MLI'};

%initialize
for counter_var = 1:numel(cell_var)
    current_var = cell_var{counter_var};
    data_out.(current_var).rate_tot_sac  = nan(1,500,8);
    data_out.(current_var).cell_ids_tot  = cell(1);
    data_out.(current_var).cs_on_ang_tot = nan(1);
    data_out.(current_var).cs_on_rho_tot = nan(1);
    data_out.(current_var).cliques       = nan(1);

    % eye traces
    data_out.(current_var).vm_tot_sac    = zeros(1,500,8);
end

% CS
% last dim 1:motor, 2:vis
data_out.CS.rate_tot_sac  = nan(1,500,8,2);
data_out.CS.cell_ids_tot  = cell(1);
data_out.CS.cs_on_ang_tot = nan(1);
data_out.CS.cs_on_rho_tot = nan(1);
data_out.CS.cliques       = nan(1);

% eye traces
data_out.CS.vm_tot_sac    = nan(1,500,8,2);

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
        data_eye     = MAF_load_eye_traces(current_path,current_sess);

        % load units info
        units_info   = MAF_extract_cell_metadata(current_path,current_sess);
        
        cell_ids     = units_info.cell_list;
        cell_types   = units_info.cell_type;
        cell_cliques = units_info.cell_clique;

        clusters     = unique(cell_cliques);
        num_clust    = numel(clusters);

        for counter_clust = 1:num_clust
            current_clust = clusters(counter_clust);
            if current_clust <= 0
                continue;
            end
            disp(['         ' num2str(current_clust)]);
            current_ind_clust = cell_cliques == current_clust;
            cell_ids_clust    = cell_ids(current_ind_clust);
            cell_types_clust  = cell_types(current_ind_clust);
            ind_cs            = ismember(cell_types_clust,{'PC'});

            num_cs_clust = sum(ind_cs);
            % skip if there is no CS in the clique
            if num_cs_clust == 0
                disp('no CS!')
                continue;
            end
            
            % find cs on average and std of the cluster
            cs_ids_clust  = cell_ids_clust(ind_cs);
            cs_on_rho_tot = nan(num_cs_clust,1);
            cs_on_ang_tot = nan(num_cs_clust,1);
            parfor counter_cs = 1:num_cs_clust
                current_cs_ids_clust = cs_ids_clust{counter_cs};
                data = MAF_load_cell(current_path,current_sess,current_cs_ids_clust);
                data_rec = data.data_recordings;
                on_data  = MAF_CS_on_analysis_fr([data_rec.eye],[data_rec.Neural_Data]);

                cs_on_rho_tot(counter_cs) = on_data.vis.rho_avg;
                cs_on_ang_tot(counter_cs) = on_data.vis.ang_avg/180*pi;
            end
            
            vec_net = mean(cs_on_rho_tot.*exp(1i*cs_on_ang_tot));
            cs_on_avg = angle(vec_net);
            cs_on_rho_avg = abs(vec_net);

            ang_bins = params.sac.ang_values/180*pi;

%             [~,clique_cs_on_bin] = min(abs(angdiff(ang_bins,cs_on_avg*ones(size(ang_bins)))));

            cs_on_avg = cs_on_avg*180/pi;

%             order_ang = circshift(1:8,1-clique_cs_on_bin);
            % randomize the cs-on direction
            order_ang = 1:8;


            % loop over SS,MLI1,MLI2
            for counter_var = 1:numel(cell_var)
                current_var = cell_var{counter_var};
                cells_ind  = ismember(cell_types_clust,cell_var_types{counter_var});

                current_cell_ids_clust   = cell_ids_clust(cells_ind);
                current_cell_types_clust = cell_types_clust(cells_ind);
                num_cell_clust           = numel(current_cell_ids_clust);
                cell_rate_tot_sac_       = nan(num_cell_clust,500,8);
                cell_cell_ids_tot_       = cell(num_cell_clust,1);
                cell_cs_on_ang_tot_      = nan(num_cell_clust,1);
                cell_cs_on_rho_tot_      = nan(num_cell_clust,1);
                cell_cliques_            = nan(num_cell_clust,1);
                cell_vm_tot_sac_         = nan(num_cell_clust,500,8);

                parfor counter_cell = 1:num_cell_clust
                    current_cell = current_cell_ids_clust{counter_cell};
                    current_cell_type = current_cell_types_clust{counter_cell};
                    disp([num2str(counter_cell),'/',num2str(num_cell_clust),': ',...
                        current_cell_type, '    ', current_cell]);
                    
                    res = load_ss_data_tuned(current_path,current_sess,current_cell,order_ang,sac_crit,data_eye);
                    cell_rate_tot_sac_(counter_cell,:,:) = res.rate_tot;
                    cell_vm_tot_sac_(counter_cell,:,:)   = res.vm_tot;
                    cell_cell_ids_tot_{counter_cell} = current_cell;
                    cell_cs_on_ang_tot_(counter_cell) = cs_on_avg;
                    cell_cs_on_rho_tot_(counter_cell) = cs_on_rho_avg;
                    cell_cliques_(counter_cell) = current_clust;
                end
                data_out.(current_var).rate_tot_sac = cat(1,data_out.(current_var).rate_tot_sac,cell_rate_tot_sac_);
                data_out.(current_var).vm_tot_sac = cat(1,data_out.(current_var).vm_tot_sac,cell_vm_tot_sac_);
                data_out.(current_var).cell_ids_tot = cat(1,data_out.(current_var).cell_ids_tot,cell_cell_ids_tot_);
                data_out.(current_var).cs_on_ang_tot = cat(1,data_out.(current_var).cs_on_ang_tot,cell_cs_on_ang_tot_);
                data_out.(current_var).cs_on_rho_tot = cat(1,data_out.(current_var).cs_on_rho_tot,cell_cs_on_rho_tot_);
                data_out.(current_var).cliques = cat(1,data_out.(current_var).cliques,cell_cliques_);
            end
            
            % complex spike
            CS_cells_ind        = ismember(cell_types_clust,{'PC'});
            cs_cell_ids_clust   = cell_ids_clust(CS_cells_ind);
            cs_cell_types_clust = cell_types_clust(CS_cells_ind);
            num_cs_clust = numel(cs_cell_ids_clust);
            CS_rate_tot_sac_ = nan(num_cs_clust,500,8,2);
            CS_vm_tot_sac_ = nan(num_cs_clust,500,8,2);
            CS_cell_ids_tot_ = cell(num_cs_clust,1);
            CS_cs_on_ang_tot_ = nan(num_cs_clust,1);
            CS_cs_on_rho_tot_ = nan(num_cs_clust,1);
            CS_cliques_       = nan(num_cs_clust,1);

            parfor counter_cs = 1:num_cs_clust
                current_cs = cs_cell_ids_clust{counter_cs};
                current_cs_type = cs_cell_types_clust{counter_cs};
                disp([num2str(counter_cs),'/',num2str(num_cs_clust),': ',...
                    current_cs_type, '    ', current_cs]);

                res = load_cs_data_tuned(current_path,current_sess,current_cs,order_ang,sac_crit,data_eye);
                CS_rate_tot_sac_(counter_cs,:,:,:) = res.rate_tot;
                CS_vm_tot_sac_(counter_cs,:,:,:) = res.vm_tot;
                CS_cell_ids_tot_{counter_cs} = current_cs;
                CS_cs_on_ang_tot_(counter_cs) = cs_on_avg;
                CS_cs_on_rho_tot_(counter_cs) = cs_on_rho_avg;
                CS_cliques_(counter_cs) = current_clust;
            end

            data_out.CS.rate_tot_sac = cat(1,data_out.CS.rate_tot_sac,CS_rate_tot_sac_);
            data_out.CS.vm_tot_sac = cat(1,data_out.CS.vm_tot_sac,CS_vm_tot_sac_);
            data_out.CS.cell_ids_tot = cat(1,data_out.CS.cell_ids_tot,CS_cell_ids_tot_);
            data_out.CS.cs_on_ang_tot = cat(1,data_out.CS.cs_on_ang_tot,CS_cs_on_ang_tot_);
            data_out.CS.cs_on_rho_tot = cat(1,data_out.CS.cs_on_rho_tot,CS_cs_on_rho_tot_);
            data_out.CS.cliques = cat(1,data_out.CS.cliques,CS_cliques_);
        end
    end
end

% remove the first empty element and save data
for counter_var = 1:numel(cell_var)
    current_var = cell_var{counter_var};
    data_out.(current_var).cell_ids_tot(1)         = [];
    data_out.(current_var).rate_tot_sac(1,:,:)     = [];
    data_out.(current_var).cs_on_ang_tot(1)        = [];
    data_out.(current_var).cs_on_rho_tot(1)        = [];
    data_out.(current_var).cliques(1)              = [];
    data_out.(current_var).vm_tot_sac(1,:,:)   = [];

    data = data_out.(current_var);
    file_name = '_rate_prim_sac_rand';
    save(fullfile(data_path,'population_data',[current_var file_name '.mat']), ...
        'data','-v7.3');
end

% data_out.CS.cell_ids_tot(1) = [];
% data_out.CS.rate_tot_sac(1,:,:,:) = [];
% data_out.CS.cs_on_ang_tot(1) = [];
% data_out.CS.cs_on_rho_tot(1) = [];
% data_out.CS.cliques(1) = [];
% 
% data = data_out.CS;

% save([pop_path,'CS_population_clique_rand.mat'],'data','-v7.3');
    
  
end

%%

function res = load_ss_data_tuned(current_path,current_sess,current_cell,order_ang,sac_crit,data_eye)

tag_id = [1];
time_span = 250;
alignment = 'vmax';

data = MAF_load_cell(current_path,current_sess,current_cell);
eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);

neural_sac_data = MAF_combine_dataset_sac([data.data_recordings],[data.rec_info],{alignment},time_span,tag_id);
eye_vm = MAF_combine_vm_traces(eye_traces,{alignment},time_span,tag_id);

neur_data.eye_vm    = eye_vm.(alignment).eye_vm;
neur_data.SS        = neural_sac_data.(alignment).SS;
sac_data.sac        = neural_sac_data.(alignment).sac;
sac_data.fix        = neural_sac_data.(alignment).fix;
sac_data.sac.has_ss = neural_sac_data.(alignment).has_ss;
sac_data_dir        = MAF_bin_by_ang(neur_data,sac_data,alignment);

ss_rate_tot_sac = zeros(500,8);
ss_vm_tot_sac   = zeros(500,8);

parfor counter_dir = 1:8
    % shifting space wrt CS on
    has_ss = sac_data_dir(order_ang(counter_dir)).sac.has_ss;
    current_amp = sac_data_dir(order_ang(counter_dir)).sac.eye_amp(has_ss);
    current_vel = sac_data_dir(order_ang(counter_dir)).sac.eye_vm_max(has_ss);
    current_dec = (sac_data_dir(order_ang(counter_dir)).sac.time_offset(has_ss) -...
        sac_data_dir(order_ang(counter_dir)).sac.time_vmax(has_ss))*1e3;
    current_acc = (sac_data_dir(order_ang(counter_dir)).sac.time_vmax(has_ss) -...
        sac_data_dir(order_ang(counter_dir)).sac.time_onset(has_ss))*1e3;
    current_raster = sac_data_dir(order_ang(counter_dir)).SS(:,has_ss);
    current_eye_vm = sac_data_dir(order_ang(counter_dir)).eye_vm(:,has_ss);

    % remove trials without spike if ss also remove invalid sacs
    poly_x = sac_crit.ms(1,:);
    poly_y = sac_crit.ms(2,:);
    [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);

    acc_dec_ = current_dec./current_acc;
    poly_x = sac_crit.ma(1,:);
    poly_y = sac_crit.ma(2,:);
    [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);
    
    % no spike, negative acc, dec, outliers in main-seq and amp-acc/dec plot
    rm_ind = (sum(current_raster)==0) | ~(in_ms | on_ms) | ~(in_acc | on_acc) | (current_acc<=0) | (current_dec<=0);

    current_raster(:,rm_ind) = [];
    current_eye_vm(:,rm_ind) = [];

    ss_rate_tot_sac(:,counter_dir) = ESN_smooth(mean(current_raster,2))*1e3;
    ss_vm_tot_sac(:,counter_dir)   = mean(current_eye_vm,2);
end

res.rate_tot = ss_rate_tot_sac;
res.vm_tot   = ss_vm_tot_sac;

end
%%
function res = load_cs_data_tuned(current_path,current_sess,current_cell,order_ang,sac_crit,data_eye)

tags = 1;

data = MAF_load_cell(current_path,current_sess,current_cell);
eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);

neural_sac_data = MAF_combine_dataset_sac([data.data_recordings],[data.rec_info],{'visual','onset'},250,tags);
eye_vm = MAF_combine_vm_traces(eye_traces,{'visual','onset'},250,tags);

% vis
neur_data_vis.eye_vm    = eye_vm.visual.eye_vm;
neur_data_vis.CS        = neural_sac_data.visual.CS;
sac_data_vis.sac        = neural_sac_data.visual.sac;
sac_data_vis.fix        = neural_sac_data.visual.fix;
sac_data_vis.sac.has_ss = neural_sac_data.visual.has_ss;
sac_data_vis.sac.has_cs = neural_sac_data.visual.has_cs;
sac_data_dir_vis        = MAF_bin_by_ang(neur_data_vis,sac_data_vis,'visual');

% sac
neur_data_sac.eye_vm    = eye_vm.onset.eye_vm;
neur_data_sac.CS        = neural_sac_data.onset.CS;
sac_data_sac.sac        = neural_sac_data.onset.sac;
sac_data_sac.fix        = neural_sac_data.onset.fix;
sac_data_sac.sac.has_cs = neural_sac_data.onset.has_cs;
sac_data_dir_sac   = MAF_bin_by_ang(neur_data_sac,sac_data_sac,'onset');

cs_rate_tot_sac = nan(500,8);
cs_rate_tot_vis = nan(500,8);

cs_vm_tot_sac = nan(500,8);
cs_vm_tot_vis = nan(500,8);

parfor counter_dir = 1:8
    % shifting space wrt CS on
    % vis
    has_cs_vis = sac_data_dir_vis(order_ang(counter_dir)).sac.has_cs;
    current_amp = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_amp(has_cs_vis);
    current_vel = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_vm_max(has_cs_vis);
    current_dec = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_offset(has_cs_vis) -...
        sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(has_cs_vis))*1e3;
    current_acc = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(has_cs_vis) -...
        sac_data_dir_vis(order_ang(counter_dir)).sac.time_onset(has_cs_vis))*1e3;
    current_raster_vis = sac_data_dir_vis(order_ang(counter_dir)).CS(:,has_cs_vis);
    current_vm_vis = sac_data_dir_vis(order_ang(counter_dir)).eye_vm(:,has_cs_vis);

    % remove trials without spike if ss also remove invalid sacs
    poly_x = sac_crit.ms(1,:);
    poly_y = sac_crit.ms(2,:);
    [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);

    acc_dec_ = current_dec./current_acc;
    poly_x = sac_crit.ma(1,:);
    poly_y = sac_crit.ma(2,:);
    [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);
    
    % negative acc, dec, outliers in main-seq and amp-acc/dec plot
    rm_ind = ~(in_ms | on_ms) | ~(in_acc | on_acc) | (current_acc<=0) | (current_dec<=0);

    current_raster_vis(:,rm_ind) = [];
    current_vm_vis(:,rm_ind) = [];

    cs_rate_tot_vis(:,counter_dir) = ESN_smooth(mean(current_raster_vis,2))*1e3;
    cs_vm_tot_vis(:,counter_dir) = mean(current_vm_vis,2);

    % sac
    has_cs_sac = sac_data_dir_sac(order_ang(counter_dir)).sac.has_cs;
    current_amp = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_amp(has_cs_sac);
    current_vel = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_vm_max(has_cs_sac);
    current_dec = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_offset(has_cs_sac) -...
        sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(has_cs_sac))*1e3;
    current_acc = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(has_cs_sac) -...
        sac_data_dir_sac(order_ang(counter_dir)).sac.time_onset(has_cs_sac))*1e3;
    current_raster_sac = sac_data_dir_sac(order_ang(counter_dir)).CS(:,has_cs_sac);
    current_vm_sac = sac_data_dir_sac(order_ang(counter_dir)).eye_vm(:,has_cs_sac);

    % remove trials without spike if ss also remove invalid sacs
    poly_x = sac_crit.ms(1,:);
    poly_y = sac_crit.ms(2,:);
    [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);

    acc_dec_ = current_dec./current_acc;
    poly_x = sac_crit.ma(1,:);
    poly_y = sac_crit.ma(2,:);
    [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);
    
    % no spike, negative acc, dec, outliers in main-seq and amp-acc/dec plot
    rm_ind = ~(in_ms | on_ms) | ~(in_acc | on_acc) | (current_acc<=0) | (current_dec<=0);

    current_raster_sac(:,rm_ind) = [];
    current_vm_sac(:,rm_ind) = [];

    cs_rate_tot_sac(:,counter_dir) = ESN_smooth(mean(current_raster_sac,2))*1e3;
    cs_vm_tot_sac(:,counter_dir) = mean(current_vm_sac,2);

end

res.rate_tot = cat(3,cs_rate_tot_vis,cs_rate_tot_sac);
res.vm_tot   = cat(3,cs_vm_tot_vis,cs_vm_tot_sac);

end