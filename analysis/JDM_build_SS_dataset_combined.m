function JDM_build_SS_dataset_combined(data_path)
% If Vel_edges is omitted or empty, we will use 3 bins (tertiles) per animal.

    % | `order_ang(i)` | Represents  |
    % | -------------- | ----------- |
    % | i = 1          | CS-on (0°)  |
    % | i = 2          | CS-on +45°  |
    % | i = 3          | CS-on +90°  |
    % | i = 4          | CS-on +135° |
    % | i = 5          | CS-on +180° |
    % | i = 6          | CS-on -135° |
    % | i = 7          | CS-on -90°  |
    % | i = 8          | CS-on -45°  |

    JDM_params_funcs
    load(fullfile(data_path,'population_data','behave_data'));
    max_vlocities = load(fullfile(data_path, 'population_data', 'sac_kinematic_prim_sac.mat'),'tag1_max_vel');
    user = 'JDM';loc  = 'ctx';

    % reward conditions
    forced_conds = [params.reward_tag1_conds;...
                    params.reward_tag4_conds];

    n_conds = size(forced_conds,1);

    % amplitude bins
    Amp_edges = [1 2 3 4 5 100];
    n_amp_bins = numel(Amp_edges)-1;
    Vel_edges = [100 200 300 400 10000];   % => Low:100–200, Mid:200–300, High:300–400, Very High:≥400 °/s

    Vel_edges = Vel_edges(:)';
    n_vel_bins = numel(Vel_edges)-1;
  
    % tag id
    tag_id = [1 4];
    path_data_monkey_sorted = params.path_data_monkey_sorted;
    animal_list  = params.animal_list;
    session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);
    num_animal   = numel(session_list);

    % containers
    cell_var       = {'SS','MLI','MLI2'};
    cell_var_types = {{'SS','PC'},'MLI','MLI2'};

    for counter_var = 1:numel(cell_var)
        current_var = cell_var{counter_var};
        % legacy (no amp/vel split): [cell x 500 x 8 x n_conds]
        data_out.(current_var).rate_tot_sac  = [];
        data_out.(current_var).vm_tot_sac    = [];
        % new amp-binned: [cell x 500 x 8 x n_conds x n_amp_bins]
        data_out.(current_var).rate_tot_sac_amp = [];
        data_out.(current_var).vm_tot_sac_amp   = [];
        data_out.(current_var).n_trials_sac_amp = []; % [cell x 8 x n_conds x n_amp_bins]
        % new vel-binned: [cell x 500 x 8 x n_conds x n_vel_bins]
        data_out.(current_var).rate_tot_sac_vel = [];
        data_out.(current_var).vm_tot_sac_vel   = [];
        data_out.(current_var).n_trials_sac_vel = []; % [cell x 8 x n_conds x n_vel_bins]
        % trial counts (overall)
        data_out.(current_var).n_trials_tot     = []; % [cell x 8 x n_conds]
        % per-cell meta
        data_out.(current_var).cell_ids_tot  = {};
        data_out.(current_var).cs_on_ang_tot = [];
        data_out.(current_var).cs_on_rho_tot = [];
        data_out.(current_var).cliques       = [];
    end

    % CS container (unchanged)
    data_out.CS.rate_tot_sac  = []; % [cell x 500 x 8 x n_conds x 2]
    data_out.CS.vm_tot_sac    = [];
    data_out.CS.cell_ids_tot  = {};
    data_out.CS.cs_on_ang_tot = [];
    data_out.CS.cs_on_rho_tot = [];
    data_out.CS.cliques       = [];

    clique_info = struct('session',{}, 'clique_id',{}, 'cs_on_avg_deg',{}, ...
                         'cs_on_rho',{}, 'cs_on_bin',{}, 'order_ang',{});

    % keep track of edges actually used per animal (helpful if using fallback)
    vel_edges_by_animal = cell(num_animal,1);

   fprintf('\n================  Building datasets for TAG [%s]  ================\n', ...
    sprintf('%d ', tag_id));


    for counter_animal = 1:num_animal
        disp(animal_list{counter_animal});
        current_path      = path_data_monkey_sorted{counter_animal};
        current_sess_list = session_list{counter_animal};
        num_sess          = numel(current_sess_list);

        current_sess_cr_ma = data_behav.sess_cr_ma{counter_animal};
        current_sess_cr_ms = data_behav.sess_cr_ms{counter_animal};

        % decide animal-level velocity edges
        if nargin >= 2 && ~isempty(Vel_edges)
            Vel_edges_local = Vel_edges;               % user-provided
        else
            % tertiles over this animal's tag1 maxima
            max_vel_sub = max_vlocities.tag1_max_vel{counter_animal,1};  
            all_vel = vertcat(max_vel_sub{:});
            q = quantile(all_vel, [1/3, 2/3]);
            if ~issorted(q) || any(~isfinite(q)), q = quantile(all_vel, [0.3334, 0.6667]); end
            if q(2) <= q(1), q(2) = q(1) + eps(q(1)); end
            Vel_edges_local = [0, q(1), q(2), inf];   % 3 bins
        end
        vel_edges_by_animal{counter_animal} = Vel_edges_local;

        for counter_sess = 1:num_sess
            current_sess = current_sess_list{counter_sess};
            disp(['     ' current_sess]);

            sac_crit.ma = current_sess_cr_ma{counter_sess};
            sac_crit.ms = current_sess_cr_ms{counter_sess};

            data_eye   = MAF_load_eye_traces(current_path,current_sess);
            units_info = MAF_extract_cell_metadata(current_path,current_sess);
            cell_ids     = units_info.cell_list;
            cell_types   = units_info.cell_type;
            cell_cliques = units_info.cell_clique;

            clusters  = unique(cell_cliques);
            num_clust = numel(clusters);

            for counter_clust = 1:num_clust
                current_clust = clusters(counter_clust);
                if current_clust <= 0, continue; end
                disp(['         Clique ' num2str(current_clust)]);

                current_ind_clust = cell_cliques == current_clust;
                cell_ids_clust    = cell_ids(current_ind_clust);
                cell_types_clust  = cell_types(current_ind_clust);

                ind_cs       = ismember(cell_types_clust,{'PC'});
                num_cs_clust = sum(ind_cs);
                if num_cs_clust == 0
                    disp('No CS, skipping clique!')
                    continue;
                end

                % compute clique CS-on average
                cs_ids_clust  = cell_ids_clust(ind_cs);
                cs_on_rho_tot = nan(num_cs_clust,1);
                cs_on_ang_tot = nan(num_cs_clust,1);
                for counter_cs = 1:num_cs_clust
                    current_cs_ids_clust = cs_ids_clust{counter_cs};
                    cell_data = MAF_load_cell(current_path,current_sess,current_cs_ids_clust);
                    data_rec  = cell_data.data_recordings;
                    on_data = MAF_CS_on_analysis_fr([data_rec.eye], [data_rec.Neural_Data],[1,4,6,7],[40,85],[-70,30]);
                    cs_on_rho_tot(counter_cs) = on_data.vis.rho_avg;
                    cs_on_ang_tot(counter_cs) = on_data.vis.ang_avg/180*pi;
                end

                vec_net        = mean(cs_on_rho_tot.*exp(1i*cs_on_ang_tot));
                cs_on_avg      = angle(vec_net);          % rad
                cs_on_rho_avg  = abs(vec_net);
                ang_bins       = params.sac.ang_values/180*pi;
                [~,clique_cs_on_bin] = min(abs(angdiff(ang_bins,cs_on_avg*ones(size(ang_bins)))));
                cs_on_avg_deg  = cs_on_avg*180/pi;
                order_ang      = circshift(1:8,1-clique_cs_on_bin);

                ci = struct;
                ci.session       = current_sess;
                ci.clique_id     = current_clust;
                ci.cs_on_avg_deg = cs_on_avg_deg;
                ci.cs_on_rho     = cs_on_rho_avg;
                ci.cs_on_bin     = clique_cs_on_bin;
                ci.order_ang     = order_ang;
                clique_info(end+1) = ci; %#ok<AGROW>

                % ---------- loop over SS / MLI / MLI2 ----------
                for counter_var = 1:numel(cell_var)
                    current_var  = cell_var{counter_var};
                    cells_ind    = ismember(cell_types_clust,cell_var_types{counter_var});
                    current_cell_ids_clust   = cell_ids_clust(cells_ind);
                    current_cell_types_clust = cell_types_clust(cells_ind);
                    num_cell_clust           = numel(current_cell_ids_clust);

                    % prealloc (legacy / amp / vel)
                    cell_rate_tot_sac_      = nan(num_cell_clust,500,8,n_conds);
                    cell_vm_tot_sac_        = nan(num_cell_clust,500,8,n_conds);

                    cell_rate_tot_sac_amp_  = nan(num_cell_clust,500,8,n_conds,n_amp_bins);
                    cell_vm_tot_sac_amp_    = nan(num_cell_clust,500,8,n_conds,n_amp_bins);
                    cell_n_trials_sac_amp_  = nan(num_cell_clust,8,n_conds,n_amp_bins);

                    cell_rate_tot_sac_vel_  = nan(num_cell_clust,500,8,n_conds,n_vel_bins);
                    cell_vm_tot_sac_vel_    = nan(num_cell_clust,500,8,n_conds,n_vel_bins);
                    cell_n_trials_sac_vel_  = nan(num_cell_clust,8,n_conds,n_vel_bins);

                    cell_n_trials_tot_      = nan(num_cell_clust,8,n_conds);

                    cell_cell_ids_tot_  = cell(num_cell_clust,1);
                    cell_cs_on_ang_tot_ = nan(num_cell_clust,1);
                    cell_cs_on_rho_tot_ = nan(num_cell_clust,1);
                    cell_cliques_       = nan(num_cell_clust,1);

                    for counter_cell = 1:num_cell_clust
                        current_cell      = current_cell_ids_clust{counter_cell};
                        current_cell_type = current_cell_types_clust{counter_cell};
                        disp([num2str(counter_cell),'/',num2str(num_cell_clust),': ', ...
                              current_cell_type, '    ', current_cell]);

                        res = load_ss_data_tuned(current_path,current_sess,current_cell, ...
                              order_ang,sac_crit,data_eye,forced_conds,Amp_edges,Vel_edges,tag_id);

                        % legacy
                        cell_rate_tot_sac_(counter_cell,:,:,:) = res.rate_tot;
                        cell_vm_tot_sac_(counter_cell,:,:,:)   = res.vm_tot;

                        % amp-binned
                        if isfield(res,'rate_tot_amp')
                            cell_rate_tot_sac_amp_(counter_cell,:,:,:,:) = res.rate_tot_amp;
                        end
                        if isfield(res,'vm_tot_amp')
                            cell_vm_tot_sac_amp_(counter_cell,:,:,:,:)   = res.vm_tot_amp;
                        end
                        if isfield(res,'n_trials_amp')
                            cell_n_trials_sac_amp_(counter_cell,:,:,:)   = res.n_trials_amp;
                        end

                        % vel-binned
                        if isfield(res,'rate_tot_vel')
                            cell_rate_tot_sac_vel_(counter_cell,:,:,:,:) = res.rate_tot_vel;
                        end
                        if isfield(res,'vm_tot_vel')
                            cell_vm_tot_sac_vel_(counter_cell,:,:,:,:)   = res.vm_tot_vel;
                        end
                        if isfield(res,'n_trials_vel')
                            cell_n_trials_sac_vel_(counter_cell,:,:,:)   = res.n_trials_vel;
                        end

                        % overall trials
                        if isfield(res,'n_trials_tot')
                            cell_n_trials_tot_(counter_cell,:,:)         = res.n_trials_tot;
                        end

                        % per-cell meta
                        cell_cell_ids_tot_{counter_cell} = current_cell;
                        cell_cs_on_ang_tot_(counter_cell) = cs_on_avg_deg;
                        cell_cs_on_rho_tot_(counter_cell) = cs_on_rho_avg;
                        cell_cliques_(counter_cell)       = current_clust;
                    end

                    % concat into global containers
                    data_out.(current_var).rate_tot_sac      = cat(1,data_out.(current_var).rate_tot_sac,     cell_rate_tot_sac_);
                    data_out.(current_var).vm_tot_sac        = cat(1,data_out.(current_var).vm_tot_sac,       cell_vm_tot_sac_);

                    data_out.(current_var).rate_tot_sac_amp  = cat(1,data_out.(current_var).rate_tot_sac_amp, cell_rate_tot_sac_amp_);
                    data_out.(current_var).vm_tot_sac_amp    = cat(1,data_out.(current_var).vm_tot_sac_amp,   cell_vm_tot_sac_amp_);
                    data_out.(current_var).n_trials_sac_amp  = cat(1,data_out.(current_var).n_trials_sac_amp, cell_n_trials_sac_amp_);

                    data_out.(current_var).rate_tot_sac_vel  = cat(1,data_out.(current_var).rate_tot_sac_vel, cell_rate_tot_sac_vel_);
                    data_out.(current_var).vm_tot_sac_vel    = cat(1,data_out.(current_var).vm_tot_sac_vel,   cell_vm_tot_sac_vel_);
                    data_out.(current_var).n_trials_sac_vel  = cat(1,data_out.(current_var).n_trials_sac_vel, cell_n_trials_sac_vel_);

                    data_out.(current_var).n_trials_tot      = cat(1,data_out.(current_var).n_trials_tot,     cell_n_trials_tot_);
                    data_out.(current_var).cell_ids_tot      = cat(1,data_out.(current_var).cell_ids_tot,     cell_cell_ids_tot_);
                    data_out.(current_var).cs_on_ang_tot     = cat(1,data_out.(current_var).cs_on_ang_tot,    cell_cs_on_ang_tot_);
                    data_out.(current_var).cs_on_rho_tot     = cat(1,data_out.(current_var).cs_on_rho_tot,    cell_cs_on_rho_tot_);
                    data_out.(current_var).cliques           = cat(1,data_out.(current_var).cliques,          cell_cliques_);
                end

                % ---------- complex spikes (unchanged) ----------
                CS_cells_ind        = ismember(cell_types_clust,{'PC'});
                cs_cell_ids_clust   = cell_ids_clust(CS_cells_ind);
                cs_cell_types_clust = cell_types_clust(CS_cells_ind);
                num_cs_clust        = numel(cs_cell_ids_clust);

                CS_rate_tot_sac_   = nan(num_cs_clust,500,8,n_conds,2);
                CS_vm_tot_sac_     = nan(num_cs_clust,500,8,n_conds,2);
                CS_cell_ids_tot_   = cell(num_cs_clust,1);
                CS_cs_on_ang_tot_  = nan(num_cs_clust,1);
                CS_cs_on_rho_tot_  = nan(num_cs_clust,1);
                CS_cliques_        = nan(num_cs_clust,1);

                for counter_cs = 1:num_cs_clust
                    current_cs      = cs_cell_ids_clust{counter_cs};
                    current_cs_type = cs_cell_types_clust{counter_cs};
                    disp([num2str(counter_cs),'/',num2str(num_cs_clust),': ', ...
                          current_cs_type, '    ', current_cs]);

                    res = load_cs_data_tuned(current_path,current_sess,current_cs, ...
                          order_ang,sac_crit,data_eye,forced_conds,tag_id);

                    CS_rate_tot_sac_(counter_cs,:,:,:,:) = res.rate_tot;
                    CS_vm_tot_sac_(counter_cs,:,:,:,:)   = res.vm_tot;
                    CS_cell_ids_tot_{counter_cs}         = current_cs;
                    CS_cs_on_ang_tot_(counter_cs)        = cs_on_avg_deg;
                    CS_cs_on_rho_tot_(counter_cs)        = cs_on_rho_avg;
                    CS_cliques_(counter_cs)              = current_clust;
                end

                data_out.CS.rate_tot_sac  = cat(1,data_out.CS.rate_tot_sac, CS_rate_tot_sac_);
                data_out.CS.vm_tot_sac    = cat(1,data_out.CS.vm_tot_sac,   CS_vm_tot_sac_);
                data_out.CS.cell_ids_tot  = cat(1,data_out.CS.cell_ids_tot, CS_cell_ids_tot_);
                data_out.CS.cs_on_ang_tot = cat(1,data_out.CS.cs_on_ang_tot,CS_cs_on_ang_tot_);
                data_out.CS.cs_on_rho_tot = cat(1,data_out.CS.cs_on_rho_tot,CS_cs_on_rho_tot_);
                data_out.CS.cliques       = cat(1,data_out.CS.cliques,      CS_cliques_);
            end
        end
    end

    % ---------------- META (and save) ----------------
    meta = struct();
    meta.tag_id         = tag_id;
    meta.forced_conds   = forced_conds;
    meta.ss_alignment   = 'vmax';
    meta.cs_alignments  = {'visual','onset'};
    meta.time_span_ms   = 250;
    meta.ang_edges_deg  = params.sac.ang_edges;
    meta.ang_values_deg = params.sac.ang_values;
    meta.animal_list    = animal_list;
    meta.session_list   = session_list;
    meta.date_built     = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    meta.clique_info    = clique_info;
    meta.amp_edges        = Amp_edges(:)';                                     
    meta.amp_bin_centers  = mean([Amp_edges(1:end-1); Amp_edges(2:end)],1);
    meta.vel_edges       = Vel_edges(:)';
    meta.vel_bin_centers = mean([Vel_edges(1:end-1); Vel_edges(2:end)],1);

    % save SS/MLI/MLI2
    for counter_var = 1:numel(cell_var)
        current_var = cell_var{counter_var};
        data = data_out.(current_var);
        file_name = [current_var '_population_clique_combined_sac'];
        save(fullfile(data_path,'population_data',sprintf('%s.mat',file_name)), ...
             'data','meta','-v7.3');
    end

    % save CS
    data = data_out.CS;
    file_name = 'CS_population_clique_combined_sac.mat';
    save(fullfile(data_path,'population_data',file_name), 'data','meta','-v7.3');
end


%% load_ss_data_tuned
function res = load_ss_data_tuned(current_path,current_sess,current_cell,...
    order_ang,sac_crit,data_eye,forced_conds,Amp_edges,Vel_edges,tag_id)

if nargin < 9 || isempty(Amp_edges)
    Amp_edges = [1 2 3 4 5 100];   % degrees
end
if nargin < 10 || isempty(Vel_edges)
    Vel_edges = [100 200 300 400 10000]; 
end

n_amp_bins = numel(Amp_edges)-1;
n_vel_bins = numel(Vel_edges)-1;

n_forc_cond = size(forced_conds,1);
n_conds     = n_forc_cond ;
time_span   = 250;
alignment   = 'vmax';

data       = MAF_load_cell(current_path,current_sess,current_cell);
eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);

neural_sac_data = JDM_combine_dataset_sac([data.data_recordings],[data.rec_info],{alignment},time_span,tag_id);
eye_vm          = MAF_combine_vm_traces(eye_traces,{alignment},time_span,tag_id);

neur_data.eye_vm    = eye_vm.(alignment).eye_vm;
neur_data.SS        = neural_sac_data.(alignment).SS;
sac_data.sac        = neural_sac_data.(alignment).sac;
sac_data.fix        = neural_sac_data.(alignment).fix;
sac_data.sac.has_ss = neural_sac_data.(alignment).has_ss;
sac_data_dir        = MAF_bin_by_ang(neur_data,sac_data,alignment);

% --- original (3D) outputs ---
ss_rate_tot_sac = nan(500,8,n_conds);
ss_vm_tot_sac   = nan(500,8,n_conds);
n_trials_tot    = zeros(8,n_conds);

% --- amplitude-binned (4D) outputs: [time x dir x cond x ampBin] ---
ss_rate_tot_sac_amp = nan(500,8,n_conds,n_amp_bins);
ss_vm_tot_sac_amp   = nan(500,8,n_conds,n_amp_bins);
n_trials_amp        = zeros(8,n_conds,n_amp_bins);

% --- velocity-binned (4D) outputs: [time x dir x cond x velBin] ---
ss_rate_tot_sac_vel = nan(500,8,n_conds,n_vel_bins);
ss_vm_tot_sac_vel   = nan(500,8,n_conds,n_vel_bins);
n_trials_vel        = zeros(8,n_conds,n_vel_bins);

for counter_dir = 1:8
    for cond_idx = 1:n_conds
        task_cond  = sac_data_dir(order_ang(counter_dir)).sac.task_cond;
        task_tag   = sac_data_dir(order_ang(counter_dir)).sac.tag;
        tgt_cond   = sac_data_dir(order_ang(counter_dir)).sac.tgt_cond;
        jump_cond  = sac_data_dir(order_ang(counter_dir)).sac.jump_cond;
        rew_cond   = sac_data_dir(order_ang(counter_dir)).sac.rew_cond;
        has_ss     = sac_data_dir(order_ang(counter_dir)).sac.has_ss;

        istuned_conds = ismember([task_cond' task_tag' tgt_cond' jump_cond' rew_cond'], ...
                                     forced_conds(cond_idx,:),'rows');
       

        idx_tuned_ = istuned_conds & has_ss';
        mask_      = idx_tuned_ > 0;
        if ~any(mask_), continue; end

        % Grab trial-wise features
        current_amp    = sac_data_dir(order_ang(counter_dir)).sac.eye_amp(mask_);
        current_vel    = sac_data_dir(order_ang(counter_dir)).sac.eye_vm_max(mask_);
        current_dec_ms = (sac_data_dir(order_ang(counter_dir)).sac.time_offset(mask_) - ...
                          sac_data_dir(order_ang(counter_dir)).sac.time_vmax(mask_))*1e3;
        current_acc_ms = (sac_data_dir(order_ang(counter_dir)).sac.time_vmax(mask_) - ...
                          sac_data_dir(order_ang(counter_dir)).sac.time_onset(mask_))*1e3;
        current_raster = sac_data_dir(order_ang(counter_dir)).SS(:,mask_);
        current_eye_vm = sac_data_dir(order_ang(counter_dir)).eye_vm(:,mask_);

        % Main-sequence gate
        poly_x = sac_crit.ms(1,:);  poly_y = sac_crit.ms(2,:);
        [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);

        % Acc/dec gate
        acc_dec_ = current_dec_ms ./ max(current_acc_ms,eps);
        poly_x = sac_crit.ma(1,:);  poly_y = sac_crit.ma(2,:);
        [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);

        % Remove: no spikes, bad gates, non-positive acc/dec
        rm_ind = (sum(current_raster)==0) | ~(in_ms | on_ms) | ~(in_acc | on_acc) | ...
                 (current_acc_ms<=0) | (current_dec_ms<=0);

        % Apply removal consistently
        current_raster(:,rm_ind) = [];
        current_eye_vm(:,rm_ind) = [];
        current_amp(rm_ind)      = [];
        current_vel(rm_ind)      = [];  % <-- important for vel bins

        if isempty(current_amp), continue; end

        % --- overall (3D) averages (kept for backward compatibility) ---
      ss_rate_tot_sac(:,counter_dir,cond_idx) = ESN_smooth(mean(current_raster,2,'omitnan'))*1e3;
      ss_vm_tot_sac(:,counter_dir,cond_idx)   = mean(current_eye_vm,2,'omitnan');

        n_trials_tot(counter_dir,cond_idx)      = size(current_raster,2);

        % --- amplitude-binned (4D) ---
        amp_bin_idx = discretize(current_amp, Amp_edges);
        for b = 1:n_amp_bins
            bin_mask = (amp_bin_idx == b);
            if ~any(bin_mask), continue; end
            ss_rate_tot_sac_amp(:,counter_dir,cond_idx,b) = ESN_smooth(mean(current_raster(:,bin_mask),2,'omitnan'))*1e3;
            ss_vm_tot_sac_amp(:,counter_dir,cond_idx,b)   =   mean(current_eye_vm(:,bin_mask),2,'omitnan');
            n_trials_amp(counter_dir,cond_idx,b)          = sum(bin_mask);
        end

        % --- velocity-binned (4D) ---
        vel_bin_idx = discretize(current_vel, Vel_edges);
        for v = 1:n_vel_bins
            vmask = (vel_bin_idx == v);
            if ~any(vmask), continue; end
            ss_rate_tot_sac_vel(:,counter_dir,cond_idx,v) =  ESN_smooth(mean(current_raster(:,vmask),2,'omitnan'))*1e3;
            ss_vm_tot_sac_vel(:,counter_dir,cond_idx,v)   = mean(current_eye_vm(:,vmask),2,'omitnan');
            n_trials_vel(counter_dir,cond_idx,v)          = sum(vmask);
        end
    end
end

% -------- outputs --------
res.rate_tot        = ss_rate_tot_sac;                 % [500 x 8 x n_conds]
res.vm_tot          = ss_vm_tot_sac;                   % [500 x 8 x n_conds]
res.n_trials_tot    = n_trials_tot;                    % [8 x n_conds]

res.rate_tot_amp    = ss_rate_tot_sac_amp;             % [500 x 8 x n_conds x n_amp_bins]
res.vm_tot_amp      = ss_vm_tot_sac_amp;               % [500 x 8 x n_conds x n_amp_bins]
res.n_trials_amp    = n_trials_amp;                    % [8 x n_conds x n_amp_bins]
res.amp_edges       = Amp_edges(:)';                   
res.amp_bin_centers = mean([Amp_edges(1:end-1); Amp_edges(2:end)],1);

res.rate_tot_vel    = ss_rate_tot_sac_vel;             % [500 x 8 x n_conds x n_vel_bins]
res.vm_tot_vel      = ss_vm_tot_sac_vel;               % [500 x 8 x n_conds x n_vel_bins]
res.n_trials_vel    = n_trials_vel;                    % [8 x n_conds x n_vel_bins]
res.vel_edges       = Vel_edges(:)';                   
res.vel_bin_centers = mean([Vel_edges(1:end-1); Vel_edges(2:end)],1);
end

%% load_cs_data_tuned
function res = load_cs_data_tuned(current_path,current_sess,current_cell,...
    order_ang,sac_crit,data_eye,force_conds,tag_id)
    n_timepoints = 500;
    n_dirs       = 8;
    n_forc_cond  = size(force_conds,1);
    n_conds      = n_forc_cond;
    
    cs_rate_tot_sac = nan(n_timepoints, n_dirs, n_conds);
    cs_rate_tot_vis = nan(n_timepoints, n_dirs, n_conds);
    
    cs_vm_tot_sac = nan(n_timepoints, n_dirs, n_conds);
    cs_vm_tot_vis = nan(n_timepoints, n_dirs, n_conds);
    
    data = MAF_load_cell(current_path,current_sess,current_cell);
    eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);
    
    neural_sac_data = JDM_combine_dataset_sac([data.data_recordings],[data.rec_info],{'visual','onset'},250,tag_id);
    eye_vm = MAF_combine_vm_traces(eye_traces,{'visual','onset'},250,tag_id);
    
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


 for counter_dir = 1:n_dirs  
    for cond_idx = 1:n_conds
        task_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.task_cond;
        task_tag = sac_data_dir_vis(order_ang(counter_dir)).sac.tag;
        tgt_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.tgt_cond;
        jmp_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.jump_cond;
        rwd_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.rew_cond;
        has_cs_vis = sac_data_dir_vis(order_ang(counter_dir)).sac.has_cs;

        istuned_vis = ismember([task_cond' task_tag' tgt_cond' jmp_cond' rwd_cond'],...
            force_conds(cond_idx,:),'rows');

        idx_tuned_vis = istuned_vis & has_cs_vis';
        mask_vis   = idx_tuned_vis > 0;
        if ~any(mask_vis), continue; end
      
        current_amp = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_amp(mask_vis);
        current_vel = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_vm_max(mask_vis);
        current_dec = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_offset(mask_vis) -...
            sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(mask_vis))*1e3;
        current_acc = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(mask_vis) -...
            sac_data_dir_vis(order_ang(counter_dir)).sac.time_onset(mask_vis))*1e3;
        current_raster_vis = sac_data_dir_vis(order_ang(counter_dir)).CS(:,mask_vis);
        current_vm_vis = sac_data_dir_vis(order_ang(counter_dir)).eye_vm(:,mask_vis);

        % remove trials without spike if ss also remove invalid sacs
        poly_x = sac_crit.ms(1,:);
        poly_y = sac_crit.ms(2,:);
        [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);
    
        acc_dec_ = current_dec./current_acc;
        poly_x = sac_crit.ma(1,:);
        poly_y = sac_crit.ma(2,:);
        [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);
        
        % negative acc, dec, outliers in main-seq and amp-acc/dec plot
        rm_ind = ~(in_ms | on_ms) | ~(in_acc | on_acc) | (current_acc<=0) | (current_dec<=0) ;
    
        current_raster_vis(:,rm_ind) = [];
        current_vm_vis(:,rm_ind) = [];
    
        cs_rate_tot_vis(:, counter_dir, cond_idx) = ESN_smooth(mean(current_raster_vis,2,'omitnan'))*1e3;
        cs_vm_tot_vis(:, counter_dir, cond_idx) =   mean(current_vm_vis,2,'omitnan');
    
        % sac
        task_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.task_cond;
        tgt_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.tgt_cond;
        task_tag = sac_data_dir_sac(order_ang(counter_dir)).sac.tag;
        jmp_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.jump_cond;
        rwd_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.rew_cond;
        has_cs_sac = sac_data_dir_sac(order_ang(counter_dir)).sac.has_cs;

       
       istuned_sac = ismember([task_cond' task_tag' tgt_cond' jmp_cond' rwd_cond'],...
            force_conds(cond_idx,:),'rows');
  
        idx_tuned_sac = istuned_sac & has_cs_sac';
        mask_sac   = idx_tuned_sac > 0;
        if ~any(mask_sac), continue; end
       
        current_amp = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_amp(mask_sac);
        current_vel = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_vm_max(mask_sac);
        current_dec = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_offset(mask_sac) -...
            sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(mask_sac))*1e3;
        current_acc = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(mask_sac) -...
            sac_data_dir_sac(order_ang(counter_dir)).sac.time_onset(mask_sac))*1e3;
        current_raster_sac = sac_data_dir_sac(order_ang(counter_dir)).CS(:,mask_sac);
        current_vm_sac = sac_data_dir_sac(order_ang(counter_dir)).eye_vm(:,mask_sac);
    
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
    
       cs_rate_tot_sac(:,counter_dir,cond_idx) = ESN_smooth(mean(current_raster_sac,2,'omitnan'))*1e3;
       cs_vm_tot_sac(:,  counter_dir,cond_idx) = mean(current_vm_sac,2,'omitnan');
    end
 end
    res.rate_tot = cat(4, cs_rate_tot_vis, cs_rate_tot_sac); % [500 x 8 x n_conds  x 2]
    res.vm_tot   = cat(4, cs_vm_tot_vis, cs_vm_tot_sac);

end

