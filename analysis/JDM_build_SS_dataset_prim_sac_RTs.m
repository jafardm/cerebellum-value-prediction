function JDM_build_SS_dataset_prim_sac_RTs(data_path)

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
    user = 'JDM'; loc  = 'ctx';

    % reward conditions
    force_conds  =  ...
         [1, 1, 1;  % High
          1, 1, 0]; % Low
    n_conds = size(force_conds,1);

    % RT bins
    RT_edges  = [60:10:200 600];
    n_RT_bins = numel(RT_edges)-1;
   
    % tag id
    tag_id = 1;
    path_data_monkey_sorted = params.path_data_monkey_sorted;
    animal_list  = params.animal_list;
    session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);
    num_animal   = numel(session_list);

    % containers
    cell_var       = {'SS','MLI','MLI2'};
    cell_var_types = {{'SS','PC'},'MLI','MLI2'};

    for counter_var = 1:numel(cell_var)
        current_var = cell_var{counter_var};

        % legacy vmax-only: [cell x 500 x 8 x n_conds]
        data_out.(current_var).rate_tot_sac       = [];
        data_out.(current_var).vm_tot_sac         = [];

        % NEW: onset-aligned (same dims)
        % %%% added onset-aligned containers %%%
        data_out.(current_var).rate_tot_sac_onset = [];
        data_out.(current_var).vm_tot_sac_onset   = [];

        % vmax RT-binned: [cell x 500 x 8 x n_conds x n_RT_bins]
        data_out.(current_var).rate_tot_sac_RT       = [];
        data_out.(current_var).vm_tot_sac_RT         = [];
        data_out.(current_var).n_trials_sac_RT       = []; % [cell x 8 x n_conds x n_RT_bins]

        % NEW: onset RT-binned
        % %%% added onset RT containers %%%
        data_out.(current_var).rate_tot_sac_RT_onset = [];
        data_out.(current_var).vm_tot_sac_RT_onset   = [];

        % trial counts (overall)
        data_out.(current_var).n_trials_tot     = []; % [cell x 8 x n_conds]

        % per-cell meta
        data_out.(current_var).cell_ids_tot  = {};
        data_out.(current_var).cs_on_ang_tot = [];
        data_out.(current_var).cs_on_rho_tot = [];
        data_out.(current_var).cliques       = [];
    end

    % CS container (unchanged)
    data_out.CS.rate_tot_sac  = []; % [cell x 500 x 8 x n_conds x n_RT_bins x 2]
    data_out.CS.vm_tot_sac    = [];
    data_out.CS.cell_ids_tot  = {};
    data_out.CS.cs_on_ang_tot = [];
    data_out.CS.cs_on_rho_tot = [];
    data_out.CS.cliques       = [];

    clique_info = struct('session',{}, 'clique_id',{}, 'cs_on_avg_deg',{}, ...
                         'cs_on_rho',{}, 'cs_on_bin',{}, 'order_ang',{});


    fprintf('\n================  Building datasets for TAG %d  ================\n', tag_id);

    for counter_animal = 1:num_animal
        disp(animal_list{counter_animal});
        current_path      = path_data_monkey_sorted{counter_animal};
        current_sess_list = session_list{counter_animal};
        num_sess          = numel(current_sess_list);

        current_sess_cr_ma = data_behav.sess_cr_ma{counter_animal};
        current_sess_cr_ms = data_behav.sess_cr_ms{counter_animal};

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

                    % prealloc (vmax & onset)
                    cell_rate_tot_sac_        = nan(num_cell_clust,500,8,n_conds);
                    cell_vm_tot_sac_          = nan(num_cell_clust,500,8,n_conds);
                    % %%% onset-aligned prealloc %%%
                    cell_rate_tot_sac_onset_  = nan(num_cell_clust,500,8,n_conds);
                    cell_vm_tot_sac_onset_    = nan(num_cell_clust,500,8,n_conds);

                    % RT-binned (vmax & onset)
                    cell_rate_tot_sac_RT_        = nan(num_cell_clust,500,8,n_conds,n_RT_bins);
                    cell_vm_tot_sac_RT_          = nan(num_cell_clust,500,8,n_conds,n_RT_bins);
                    % %%% onset RT prealloc %%%
                    cell_rate_tot_sac_RT_onset_  = nan(num_cell_clust,500,8,n_conds,n_RT_bins);
                    cell_vm_tot_sac_RT_onset_    = nan(num_cell_clust,500,8,n_conds,n_RT_bins);

                    cell_n_trials_sac_RT_  = nan(num_cell_clust,8,n_conds,n_RT_bins);
                    cell_n_trials_tot_      = nan(num_cell_clust,8,n_conds);

                    cell_cell_ids_tot_  = cell(num_cell_clust,1);
                    cell_cs_on_ang_tot_ = nan(num_cell_clust,1);
                    cell_cs_on_rho_tot_ = nan(num_cell_clust,1);
                    cell_cliques_       = nan(num_cell_clust,1);

                    parfor counter_cell = 1:num_cell_clust
                        current_cell      = current_cell_ids_clust{counter_cell};
                        current_cell_type = current_cell_types_clust{counter_cell};
                        disp([num2str(counter_cell),'/',num2str(num_cell_clust),': ', ...
                              current_cell_type, '    ', current_cell]);

                        res = load_ss_data_tuned(current_path,current_sess,current_cell, ...
                              order_ang,sac_crit,data_eye,force_conds,RT_edges,tag_id);

                        % ------------ vmax-aligned (legacy) ------------
                        cell_rate_tot_sac_(counter_cell,:,:,:) = res.rate_tot;
                        cell_vm_tot_sac_(counter_cell,:,:,:)   = res.vm_tot;

                        if isfield(res,'rate_tot_RT')
                            cell_rate_tot_sac_RT_(counter_cell,:,:,:,:) = res.rate_tot_RT;
                        end
                        if isfield(res,'vm_tot_RT')
                            cell_vm_tot_sac_RT_(counter_cell,:,:,:,:)   = res.vm_tot_RT;
                        end

                        % ------------ onset-aligned (new) ------------
                        if isfield(res,'rate_tot_onset')
                            cell_rate_tot_sac_onset_(counter_cell,:,:,:) = res.rate_tot_onset;
                        end
                        if isfield(res,'vm_tot_onset')
                            cell_vm_tot_sac_onset_(counter_cell,:,:,:)   = res.vm_tot_onset;
                        end
                        if isfield(res,'rate_tot_RT_onset')
                            cell_rate_tot_sac_RT_onset_(counter_cell,:,:,:,:) = res.rate_tot_RT_onset;
                        end
                        if isfield(res,'vm_tot_RT_onset')
                            cell_vm_tot_sac_RT_onset_(counter_cell,:,:,:,:)   = res.vm_tot_RT_onset;
                        end

                        % ------------ trial counts ------------
                        if isfield(res,'n_trials_tot')
                            cell_n_trials_tot_(counter_cell,:,:) = res.n_trials_tot;  % [8 x n_conds]
                        end
                        if isfield(res,'n_trials_RT')
                            cell_n_trials_sac_RT_(counter_cell,:,:,:) = res.n_trials_RT; % [8 x n_conds x n_RT_bins]
                        end

                        % per-cell meta
                        cell_cell_ids_tot_{counter_cell}  = current_cell;
                        cell_cs_on_ang_tot_(counter_cell) = cs_on_avg_deg;
                        cell_cs_on_rho_tot_(counter_cell) = cs_on_rho_avg;
                        cell_cliques_(counter_cell)       = current_clust;
                    end

                    % concat into global containers
                    data_out.(current_var).rate_tot_sac       = cat(1,data_out.(current_var).rate_tot_sac,      cell_rate_tot_sac_);
                    data_out.(current_var).vm_tot_sac         = cat(1,data_out.(current_var).vm_tot_sac,        cell_vm_tot_sac_);

                    % %%% save onset-aligned into global struct %%%
                    data_out.(current_var).rate_tot_sac_onset = cat(1,data_out.(current_var).rate_tot_sac_onset, cell_rate_tot_sac_onset_);
                    data_out.(current_var).vm_tot_sac_onset   = cat(1,data_out.(current_var).vm_tot_sac_onset,   cell_vm_tot_sac_onset_);

                    data_out.(current_var).rate_tot_sac_RT       = cat(1,data_out.(current_var).rate_tot_sac_RT,      cell_rate_tot_sac_RT_);
                    data_out.(current_var).vm_tot_sac_RT         = cat(1,data_out.(current_var).vm_tot_sac_RT,        cell_vm_tot_sac_RT_);
                    data_out.(current_var).rate_tot_sac_RT_onset = cat(1,data_out.(current_var).rate_tot_sac_RT_onset,cell_rate_tot_sac_RT_onset_);
                    data_out.(current_var).vm_tot_sac_RT_onset   = cat(1,data_out.(current_var).vm_tot_sac_RT_onset,  cell_vm_tot_sac_RT_onset_);

                    data_out.(current_var).n_trials_sac_RT  = cat(1,data_out.(current_var).n_trials_sac_RT, cell_n_trials_sac_RT_);
                    data_out.(current_var).n_trials_tot     = cat(1,data_out.(current_var).n_trials_tot,    cell_n_trials_tot_);

                    data_out.(current_var).cell_ids_tot     = cat(1,data_out.(current_var).cell_ids_tot,    cell_cell_ids_tot_);
                    data_out.(current_var).cs_on_ang_tot    = cat(1,data_out.(current_var).cs_on_ang_tot,   cell_cs_on_ang_tot_);
                    data_out.(current_var).cs_on_rho_tot    = cat(1,data_out.(current_var).cs_on_rho_tot,   cell_cs_on_rho_tot_);
                    data_out.(current_var).cliques          = cat(1,data_out.(current_var).cliques,         cell_cliques_);
                end

                % ---------- complex spikes   ----------
                CS_cells_ind        = ismember(cell_types_clust,{'PC'});
                cs_cell_ids_clust   = cell_ids_clust(CS_cells_ind);
                cs_cell_types_clust = cell_types_clust(CS_cells_ind);
                num_cs_clust        = numel(cs_cell_ids_clust);

                CS_rate_tot_sac_   = nan(num_cs_clust,500,8,n_conds,n_RT_bins,2);
                CS_vm_tot_sac_     = nan(num_cs_clust,500,8,n_conds,n_RT_bins,2);
                CS_cell_ids_tot_   = cell(num_cs_clust,1);
                CS_cs_on_ang_tot_  = nan(num_cs_clust,1);
                CS_cs_on_rho_tot_  = nan(num_cs_clust,1);
                CS_cliques_        = nan(num_cs_clust,1);

                for counter_cs = 1:num_cs_clust
                    current_cs      = cs_cell_ids_clust{counter_cs};
                    current_cs_type = cs_cell_types_clust{counter_cs};
                    disp([num2str(counter_cs),'/',num2str(num_cs_clust),': ', ...
                          current_cs_type, '    ', current_cs]);

                    res = load_cs_data_tuned_RT(current_path,current_sess,current_cs, ...
                          order_ang,sac_crit,data_eye,force_conds,RT_edges,tag_id);

                    CS_rate_tot_sac_(counter_cs,:,:,:,:,:) = res.rate_tot;
                    CS_vm_tot_sac_(counter_cs,:,:,:,:,:)   = res.vm_tot;
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
    meta.forced_conds   = force_conds;
    % %%% updated to reflect both SS alignments %%%
    meta.ss_alignments  = {'vmax','onset'};
    meta.cs_alignments  = {'visual','onset'};
    meta.time_span_ms   = 250;
    meta.ang_edges_deg  = params.sac.ang_edges;
    meta.ang_values_deg = params.sac.ang_values;
    meta.animal_list    = animal_list;
    meta.session_list   = session_list;
    meta.date_built     = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    meta.clique_info    = clique_info;
    meta.RT_edges       = RT_edges';                                     

    % save SS/MLI/MLI2
    for counter_var = 1:numel(cell_var)
        current_var = cell_var{counter_var};
        data = data_out.(current_var);
        file_name = [current_var '_population_clique_prim_sac_RT'];
        save(fullfile(data_path,'population_data',sprintf('%s.mat',file_name)), ...
             'data','meta','-v7.3');
    end

    % save CS
    data = data_out.CS;
    file_name = 'CS_population_clique_prim_sac_RT.mat';
    save(fullfile(data_path,'population_data',file_name), 'data','meta','-v7.3');
end


%% load_ss_data_tuned
function res = load_ss_data_tuned(current_path,current_sess,current_cell,...
    order_ang,sac_crit,data_eye,force_conds,RT_edges,tag_id)

    n_RT_bins   = numel(RT_edges)-1;
    n_forc_cond = size(force_conds,1);
    n_conds     = n_forc_cond;
    time_span   = 250;
 
    % -----------------------------
    % Load data & build sac datasets
    % -----------------------------
    data       = MAF_load_cell(current_path,current_sess,current_cell);
    eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);

    % Both vmax- and onset-aligned data
    neural_sac_data = JDM_combine_dataset_sac([data.data_recordings], ...
                                              [data.rec_info], ...
                                              {'vmax','onset'}, ...
                                              time_span, tag_id);
    eye_vm          = MAF_combine_vm_traces(eye_traces, ...
                                            {'vmax','onset'}, ...
                                            time_span, tag_id);

    % ---------- vmax-aligned ----------
    neur_data_vmax.eye_vm    = eye_vm.vmax.eye_vm;
    neur_data_vmax.SS        = neural_sac_data.vmax.SS;
    sac_data_vmax.sac        = neural_sac_data.vmax.sac;
    sac_data_vmax.fix        = neural_sac_data.vmax.fix;
    sac_data_vmax.sac.has_ss = neural_sac_data.vmax.has_ss;
    sac_data_vmax_dir        = MAF_bin_by_ang(neur_data_vmax, sac_data_vmax, 'vmax');

    % ---------- onset-aligned ----------
    neur_data_onset.eye_vm    = eye_vm.onset.eye_vm;
    neur_data_onset.SS        = neural_sac_data.onset.SS;
    sac_data_onset.sac        = neural_sac_data.onset.sac;
    sac_data_onset.fix        = neural_sac_data.onset.fix;
    sac_data_onset.sac.has_ss = neural_sac_data.onset.has_ss;
    sac_data_onset_dir        = MAF_bin_by_ang(neur_data_onset, sac_data_onset, 'onset');

    % ==========================================
    % Preallocate containers
    % time x dir x cond (and x RTbin)
    % ==========================================
    % --- vmax-aligned (original) ---
    ss_rate_tot_sac_vmax   = nan(500,8,n_conds);
    ss_vm_tot_sac_vmax     = nan(500,8,n_conds);
    ss_rate_tot_sac_RT_vmax = nan(500,8,n_conds,n_RT_bins);
    ss_vm_tot_sac_RT_vmax   = nan(500,8,n_conds,n_RT_bins);

    % --- onset-aligned (new) ---
    ss_rate_tot_sac_onset   = nan(500,8,n_conds);
    ss_vm_tot_sac_onset     = nan(500,8,n_conds);
    ss_rate_tot_sac_RT_onset = nan(500,8,n_conds,n_RT_bins);
    ss_vm_tot_sac_RT_onset   = nan(500,8,n_conds,n_RT_bins);

    % --- trial counts (common) ---
    n_trials_tot = zeros(8,n_conds);
    n_trials_RT  = zeros(8,n_conds,n_RT_bins);

    % ==========================================
    % Loop over direction and condition
    % ==========================================
    for counter_dir = 1:8
        for cond_idx = 1:n_conds

            % For convenience
            sac_dir_vmax   = sac_data_vmax_dir(order_ang(counter_dir));
            sac_dir_onset  = sac_data_onset_dir(order_ang(counter_dir));

            task_cond  = sac_dir_vmax.sac.task_cond;
            task_tag   = sac_dir_vmax.sac.tag;
            tgt_cond   = sac_dir_vmax.sac.tgt_cond;
            has_ss     = sac_dir_vmax.sac.has_ss;

            istuned_conds = ismember([task_cond' task_tag' tgt_cond'], ...
                                      force_conds(cond_idx,:),'rows');

            idx_tuned_ = istuned_conds & has_ss';
            mask_      = idx_tuned_ > 0;
            if ~any(mask_), continue; end

            % -------- trial-wise features (vmax-aligned sac fields) --------
            current_amp    = sac_dir_vmax.sac.eye_amp(mask_);
            current_vel    = sac_dir_vmax.sac.eye_vm_max(mask_);
            current_dec_ms = (sac_dir_vmax.sac.time_offset(mask_) - ...
                              sac_dir_vmax.sac.time_vmax(mask_))*1e3;
            current_acc_ms = (sac_dir_vmax.sac.time_vmax(mask_) - ...
                              sac_dir_vmax.sac.time_onset(mask_))*1e3;

            % Rasters & vm for BOTH alignments
            current_raster_vmax  = sac_dir_vmax.SS(:,mask_);
            current_eye_vm_vmax  = sac_dir_vmax.eye_vm(:,mask_);
            current_raster_onset = sac_dir_onset.SS(:,mask_);
            current_eye_vm_onset = sac_dir_onset.eye_vm(:,mask_);

            % --- Reaction time (visual → saccade onset), tuned trials only ---
            time_onset  = sac_dir_vmax.sac.time_onset(mask_);
            time_visual = sac_dir_vmax.sac.time_visual(mask_);
            current_rt  = (time_onset - time_visual)*1e3;   % ms

            % ============================
            % Main-sequence gate
            % ============================
            poly_x = sac_crit.ms(1,:);  poly_y = sac_crit.ms(2,:);
            [in_ms,on_ms] = inpolygon(log10(current_amp), ...
                                      log10(current_vel), ...
                                      poly_x, poly_y);

            % Acc/dec gate
            acc_dec_ = current_dec_ms ./ max(current_acc_ms,eps);
            poly_x = sac_crit.ma(1,:);  poly_y = sac_crit.ma(2,:);
            [in_acc,on_acc] = inpolygon(log10(current_amp), ...
                                        log10(acc_dec_), ...
                                        poly_x, poly_y);

            % Remove: no spikes, bad gates, non-positive acc/dec
            rm_ind = (sum(current_raster_vmax)==0) | ...
                     ~(in_ms | on_ms) | ...
                     ~(in_acc | on_acc) | ...
                     (current_acc_ms<=0) | ...
                     (current_dec_ms<=0);

            % Apply removal consistently (both alignments + kinematics/RT)
            current_raster_vmax(:,rm_ind)  = [];
            current_eye_vm_vmax(:,rm_ind)  = [];
            current_raster_onset(:,rm_ind) = [];
            current_eye_vm_onset(:,rm_ind) = [];

            current_amp(rm_ind)      = [];
            current_vel(rm_ind)      = [];
            current_dec_ms(rm_ind)   = [];
            current_acc_ms(rm_ind)   = [];
            current_rt(rm_ind)       = [];

            if isempty(current_amp), continue; end

            % ============================
            % Overall averages (3D)
            % ============================
            % vmax-aligned
            ss_rate_tot_sac_vmax(:,counter_dir,cond_idx) = ...
                ESN_smooth(mean(current_raster_vmax,2,'omitnan'))*1e3;
            ss_vm_tot_sac_vmax(:,counter_dir,cond_idx)   = ...
                mean(current_eye_vm_vmax,2,'omitnan');

            % onset-aligned
            ss_rate_tot_sac_onset(:,counter_dir,cond_idx) = ...
                ESN_smooth(mean(current_raster_onset,2,'omitnan'))*1e3;
            ss_vm_tot_sac_onset(:,counter_dir,cond_idx)   = ...
                mean(current_eye_vm_onset,2,'omitnan');

            % trial count (same for both alignments)
            n_trials_tot(counter_dir,cond_idx) = size(current_raster_vmax,2);

            % ============================
            % RT-binned averages (4D)
            % ============================
            rt_bin_idx = discretize(current_rt, RT_edges);   % 1..n_RT_bins

            for iRT = 1:n_RT_bins
                bin_mask = (rt_bin_idx == iRT);
                if ~any(bin_mask), continue; end

                % vmax-aligned
                ss_rate_tot_sac_RT_vmax(:,counter_dir,cond_idx,iRT) = ...
                    ESN_smooth(mean(current_raster_vmax(:,bin_mask),2,'omitnan'))*1e3;
                ss_vm_tot_sac_RT_vmax(:,counter_dir,cond_idx,iRT)   = ...
                    mean(current_eye_vm_vmax(:,bin_mask),2,'omitnan');

                % onset-aligned
                ss_rate_tot_sac_RT_onset(:,counter_dir,cond_idx,iRT) = ...
                    ESN_smooth(mean(current_raster_onset(:,bin_mask),2,'omitnan'))*1e3;
                ss_vm_tot_sac_RT_onset(:,counter_dir,cond_idx,iRT)   = ...
                    mean(current_eye_vm_onset(:,bin_mask),2,'omitnan');

                % trial counts
                n_trials_RT(counter_dir,cond_idx,iRT) = sum(bin_mask);
            end
        end
    end

    % ============================
    % Outputs
    % ============================
    % --- vmax-aligned (backward compatible names) ---
    res.rate_tot     = ss_rate_tot_sac_vmax;            % [500 x 8 x n_conds]
    res.vm_tot       = ss_vm_tot_sac_vmax;              % [500 x 8 x n_conds]
    res.n_trials_tot = n_trials_tot;                    % [8 x n_conds]

    res.rate_tot_RT  = ss_rate_tot_sac_RT_vmax;         % [500 x 8 x n_conds x n_RT_bins]
    res.vm_tot_RT    = ss_vm_tot_sac_RT_vmax;           % [500 x 8 x n_conds x n_RT_bins]
    res.n_trials_RT  = n_trials_RT;                     % [8 x n_conds x n_RT_bins]

    % --- onset-aligned (new) ---
    res.rate_tot_onset    = ss_rate_tot_sac_onset;      % [500 x 8 x n_conds]
    res.vm_tot_onset      = ss_vm_tot_sac_onset;        % [500 x 8 x n_conds]

    res.rate_tot_RT_onset = ss_rate_tot_sac_RT_onset;   % [500 x 8 x n_conds x n_RT_bins]
    res.vm_tot_RT_onset   = ss_vm_tot_sac_RT_onset;     % [500 x 8 x n_conds x n_RT_bins]

end


%% load_cs_data_tuned
function res = load_cs_data_tuned_RT(current_path,current_sess,current_cell,...
    order_ang,sac_crit,data_eye,force_conds,rt_edges,tags)
    JDM_params_funcs;
  
    
    n_timepoints = 500;
    n_dirs       = 8;
    n_forc_cond  = size(force_conds,1);
    n_conds      = n_forc_cond ;
    n_rt_bins = numel(rt_edges) - 1;
    
    cs_rate_tot_sac    = nan(n_timepoints, n_dirs, n_conds,n_rt_bins);
    cs_rate_tot_vis    = nan(n_timepoints, n_dirs, n_conds,n_rt_bins);
    
    cs_vm_tot_sac    = nan(n_timepoints, n_dirs, n_conds,n_rt_bins);
    cs_vm_tot_vis    = nan(n_timepoints, n_dirs, n_conds,n_rt_bins);
    
    data       = MAF_load_cell(current_path,current_sess,current_cell);
    eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);
    
    %neural_sac_data2 = MAF_combine_dataset_sac([data.data_recordings],[data.rec_info],{'visual','onset'},250,tags);
    neural_sac_data = JDM_combine_dataset_sac([data.data_recordings],[data.rec_info],{'visual','onset'},250,tags);
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
    
    %parfor counter_dir = 1:8
     for counter_dir = 1:n_dirs  
        for cond_idx = 1:n_conds
           task_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.task_cond;
           task_tag = sac_data_dir_vis(order_ang(counter_dir)).sac.tag;
           tgt_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.tgt_cond;
           has_cs_vis = sac_data_dir_vis(order_ang(counter_dir)).sac.has_cs;
           istuned_vis = ismember([task_cond' task_tag' tgt_cond'],force_conds(cond_idx,:),'rows');
           % compute RT
           time_onset = sac_data_dir_vis(order_ang(counter_dir)).sac.time_onset;
           time_visual = sac_data_dir_vis(order_ang(counter_dir)).sac.time_visual; 
           current_rt = (time_onset - time_visual)*1e3; % ms
           rt_bin_vis = discretize(current_rt, rt_edges);
    
           for iRT = 1:n_rt_bins
            trials_in_bin = (rt_bin_vis == iRT);
            % vis
            idx_tuned_vis = istuned_vis(:) & has_cs_vis(:) & trials_in_bin(:);
    
            if ~any(idx_tuned_vis)
                cs_rate_tot_vis(:, counter_dir, cond_idx, iRT) = nan;
                cs_vm_tot_vis(:, counter_dir, cond_idx, iRT)   = nan;
                continue
            end
    
            current_amp = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_amp(idx_tuned_vis);
            current_vel = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_vm_max(idx_tuned_vis);
            current_dec = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_offset(idx_tuned_vis) -...
                sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(idx_tuned_vis))*1e3;
            current_acc = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(idx_tuned_vis) -...
                sac_data_dir_vis(order_ang(counter_dir)).sac.time_onset(idx_tuned_vis))*1e3;
            current_raster_vis = sac_data_dir_vis(order_ang(counter_dir)).CS(:,idx_tuned_vis);
            current_vm_vis = sac_data_dir_vis(order_ang(counter_dir)).eye_vm(:,idx_tuned_vis);
    
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
    
            cs_rate_tot_vis(:, counter_dir, cond_idx,iRT) = ESN_smooth(mean(current_raster_vis,2))*1e3;
            cs_vm_tot_vis(:, counter_dir, cond_idx,iRT) = mean(current_vm_vis,2);
        
           % sac
           task_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.task_cond;
           tgt_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.tgt_cond;
           task_tag = sac_data_dir_sac(order_ang(counter_dir)).sac.tag;
           has_cs_sac = sac_data_dir_sac(order_ang(counter_dir)).sac.has_cs;
           istuned_sac = ismember([task_cond' task_tag' tgt_cond'],force_conds(cond_idx,:),'rows');
    
           time_onset_sac  = sac_data_dir_sac(order_ang(counter_dir)).sac.time_onset;
           time_visual_sac = sac_data_dir_sac(order_ang(counter_dir)).sac.time_visual;
           current_rt_sac  = (time_onset_sac - time_visual_sac) * 1e3;
           rt_bin_sac      = discretize(current_rt_sac, rt_edges);
           trials_in_bin_sac = (rt_bin_sac == iRT);
           idx_tuned_sac = istuned_sac(:) & has_cs_sac(:) & trials_in_bin_sac(:);
    
           if ~any(idx_tuned_sac)
            cs_rate_tot_sac(:, counter_dir, cond_idx, iRT) = nan;
            cs_vm_tot_sac(:, counter_dir, cond_idx, iRT)   = nan;
            continue
           end
    
            current_amp = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_amp(idx_tuned_sac);
            current_vel = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_vm_max(idx_tuned_sac);
            current_dec = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_offset(idx_tuned_sac) -...
                sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(idx_tuned_sac))*1e3;
            current_acc = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(idx_tuned_sac) -...
                sac_data_dir_sac(order_ang(counter_dir)).sac.time_onset(idx_tuned_sac))*1e3;
            current_raster_sac = sac_data_dir_sac(order_ang(counter_dir)).CS(:,idx_tuned_sac);
            current_vm_sac = sac_data_dir_sac(order_ang(counter_dir)).eye_vm(:,idx_tuned_sac);
        
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
      
            cs_rate_tot_sac(:, counter_dir, cond_idx,iRT) = ESN_smooth(mean(current_raster_sac,2))*1e3;
            cs_vm_tot_sac(:, counter_dir, cond_idx,iRT)   = mean(current_vm_sac,2);
          end
        end % cond end
     end % dir end
    res.rate_tot = cat(5, cs_rate_tot_vis, cs_rate_tot_sac); % [500 x 8 x n_conds x nRTbins x 2]
    res.vm_tot   = cat(5, cs_vm_tot_vis, cs_vm_tot_sac);

end% function end