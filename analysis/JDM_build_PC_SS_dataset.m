function JDM_build_PC_SS_dataset(data_path)
    JDM_params_funcs
    S = load(fullfile(data_path,'population_data','behave_data'));
    user = 'JDM'; loc = 'ctx';

    % --- tags we want ---
    tag_list = [1,4];  

    path_data_monkey_sorted = params.path_data_monkey_sorted;
    animal_list = params.animal_list;
    session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);
    num_animal = numel(session_list);

    % =========================================================
    % Loop over tags separately
    % =========================================================
    for tag_id = tag_list
        fprintf('\n================  Building dataset for TAG %d (saccades)  ================\n', tag_id);

        % initialize container
        cell_var = {'PC','SS'};
        data_out = struct();
        for counter_var = 1:numel(cell_var)
            current_var = cell_var{counter_var};
            data_out.(current_var).rate_tot_sac = [];   % [cell × 500 × 8 × 4]
            data_out.(current_var).vm_tot_sac   = [];   % [cell × 500 × 8 × 4]
            data_out.(current_var).cell_ids_tot = {};
        end

        % =====================================================
        % Iterate over animals & sessions
        % =====================================================
        for counter_animal = 1:num_animal
            disp(animal_list{counter_animal});
            current_path = path_data_monkey_sorted{counter_animal};
            current_sess_list = session_list{counter_animal};
            num_sess = numel(current_sess_list);

            current_sess_cr_ma = S.data_behav.sess_cr_ma{counter_animal};
            current_sess_cr_ms = S.data_behav.sess_cr_ms{counter_animal};

            for counter_sess = 1:num_sess
                current_sess = current_sess_list{counter_sess};
                disp(['     ' current_sess]);

                sac_crit.ma = current_sess_cr_ma{counter_sess};
                sac_crit.ms = current_sess_cr_ms{counter_sess};

                % load eye traces
                data_eye = MAF_load_eye_traces(current_path,current_sess);

                % load units
                units_info = MAF_extract_cell_metadata(current_path,current_sess);
                cell_ids   = units_info.cell_list;
                cell_types = units_info.cell_type;

                % find PC or SS cells
                keep_idx = find(ismember(cell_types,cell_var));
                if isempty(keep_idx)
                    fprintf('---- No SS/PC in session %s-%s\n', ...
                        animal_list{counter_animal}, current_sess);
                    continue
                end

                num_cells = numel(keep_idx);

                % preallocate temp storage for this session
                tmp_rate = nan(num_cells,500,8,4);
                tmp_vm   = nan(num_cells,500,8,4);
                tmp_ids  = cell(num_cells,1);
                tmp_type = cell(num_cells,1);

                % loop cells
                parfor counter_cell = 1:num_cells
                    current_cell = cell_ids{keep_idx(counter_cell)};
                    current_type = cell_types{keep_idx(counter_cell)};

                    disp([num2str(counter_cell),'/',num2str(num_cells),': ',...
                        current_type,' ',current_cell]);

                    res = load_ss_data_tuned(current_path,current_sess,current_cell, ...
                                              1:8,sac_crit,data_eye,tag_id);

                    tmp_rate(counter_cell,:,:,:) = res.rate_tot;
                    tmp_vm(counter_cell,:,:,:)   = res.vm_tot;
                    tmp_ids{counter_cell}        = current_cell;
                    tmp_type{counter_cell}       = current_type;
                end

                % append to output
                for counter_var = 1:numel(cell_var)
                    current_var = cell_var{counter_var};
                    mask = strcmp(tmp_type,current_var);
                    if any(mask)
                        data_out.(current_var).rate_tot_sac = cat(1,data_out.(current_var).rate_tot_sac,tmp_rate(mask,:,:,:));
                        data_out.(current_var).vm_tot_sac   = cat(1,data_out.(current_var).vm_tot_sac,tmp_vm(mask,:,:,:));
                        data_out.(current_var).cell_ids_tot = cat(1,data_out.(current_var).cell_ids_tot,tmp_ids(mask));
                    end
                end
            end
        end

        % ---------------- META ----------------
        meta = struct();
        meta.tag_id         = tag_id;
        meta.ss_alignment   = 'vmax';
        meta.time_span_ms   = 250;
        meta.animal_list    = animal_list;
        meta.session_list   = session_list;
        meta.sac_crit_ma    = S.data_behav.sess_cr_ma;
        meta.sac_crit_ms    = S.data_behav.sess_cr_ms;
        meta.date_built     = datestr(now,'yyyy-mm-dd HH:MM:SS');
        meta.reward_labels  = {'HH','HL','LL','LH'};

        % save separate file for each tag
        save_name = sprintf('PC_SS_dataset_tag%d.mat',tag_id);
        save(fullfile(data_path,'population_data',save_name),'data_out','meta','-v7.3');
    end
end


%% load_ss_data_tuned
function res = load_ss_data_tuned(current_path,current_sess,current_cell,order_ang,sac_crit,data_eye,tag_id)

time_span = 250;
alignment = 'vmax';

% -------------------------------
% Select condition matrix by tag
% -------------------------------
switch tag_id
    case 1
        conds = [ ...
            1, 1, 1, 0, 1;  % HH
            1, 1, 1, 0, 0;  % HL
            1, 1, 0, 0, 0;  % LL
            1, 1, 0, 0, 1]; % LH
    case 4
        conds = [ ...
            1, 4, 1, 1, 1;  % HH
            1, 4, 1, 1, 0;  % HL
            1, 4, 0, 1, 0;  % LL
            1, 4, 0, 1, 1]; % LH
    otherwise
        error('Tag %d not implemented in load_ss_data_tuned', tag_id);
end

% -------------------------------
% Load neural & eye data
% -------------------------------
data = MAF_load_cell(current_path,current_sess,current_cell);
eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);

neural_sac_data = MAF_combine_dataset_sac([data.data_recordings], ...
                                          [data.rec_info],{alignment}, ...
                                          time_span,tag_id);
eye_vm = MAF_combine_vm_traces(eye_traces,{alignment},time_span,tag_id);

neur_data.eye_vm    = eye_vm.(alignment).eye_vm;
neur_data.SS        = neural_sac_data.(alignment).SS;
sac_data.sac        = neural_sac_data.(alignment).sac;
sac_data.fix        = neural_sac_data.(alignment).fix;
sac_data.sac.has_ss = neural_sac_data.(alignment).has_ss;
sac_data_dir        = MAF_bin_by_ang(neur_data,sac_data,alignment);

% -------------------------------
% Preallocate [time × dir × cond]
% -------------------------------
nTime = 500; nDir = 8; nCond = size(conds,1);
ss_rate_tot_sac = zeros(nTime,nDir,nCond);
ss_vm_tot_sac   = zeros(nTime,nDir,nCond);

% -------------------------------
% Loop over directions and reward conds
% -------------------------------
for counter_dir = 1:nDir
    for r = 1:nCond
        cond_mask = ...
            sac_data_dir(order_ang(counter_dir)).sac.task_cond == conds(r,1) & ...
            sac_data_dir(order_ang(counter_dir)).sac.tag       == conds(r,2) & ...
            sac_data_dir(order_ang(counter_dir)).sac.tgt_cond  == conds(r,3) & ...
            sac_data_dir(order_ang(counter_dir)).sac.jump_cond == conds(r,4) & ...
            sac_data_dir(order_ang(counter_dir)).sac.rew_cond  == conds(r,5);

        has_ss       = sac_data_dir(order_ang(counter_dir)).sac.has_ss & cond_mask;
        current_amp  = sac_data_dir(order_ang(counter_dir)).sac.eye_amp(has_ss);
        current_vel  = sac_data_dir(order_ang(counter_dir)).sac.eye_vm_max(has_ss);
        current_dec  = (sac_data_dir(order_ang(counter_dir)).sac.time_offset(has_ss) - ...
                        sac_data_dir(order_ang(counter_dir)).sac.time_vmax(has_ss))*1e3;
        current_acc  = (sac_data_dir(order_ang(counter_dir)).sac.time_vmax(has_ss) - ...
                        sac_data_dir(order_ang(counter_dir)).sac.time_onset(has_ss))*1e3;
        current_raster = sac_data_dir(order_ang(counter_dir)).SS(:,has_ss);
        current_eye_vm = sac_data_dir(order_ang(counter_dir)).eye_vm(:,has_ss);

        % --- apply saccade criteria ---
        poly_x = sac_crit.ms(1,:); poly_y = sac_crit.ms(2,:);
        [in_ms,on_ms] = inpolygon(log10(current_amp),log10(current_vel),poly_x,poly_y);

        acc_dec_ = current_dec./current_acc;
        poly_x = sac_crit.ma(1,:); poly_y = sac_crit.ma(2,:);
        [in_acc,on_acc] = inpolygon(log10(current_amp),log10(acc_dec_),poly_x,poly_y);

        % remove bad trials
        rm_ind = (sum(current_raster)==0) | ...
                 ~(in_ms | on_ms) | ...
                 ~(in_acc | on_acc) | ...
                 (current_acc<=0) | (current_dec<=0);

        current_raster(:,rm_ind) = [];
        current_eye_vm(:,rm_ind) = [];

        % --- average spike rate & velocity ---
        if ~isempty(current_raster)
            ss_rate_tot_sac(:,counter_dir,r) = ESN_smooth(mean(current_raster,2))*1e3;
            ss_vm_tot_sac(:,counter_dir,r)   = mean(current_eye_vm,2);
        else
            ss_rate_tot_sac(:,counter_dir,r) = nan;
            ss_vm_tot_sac(:,counter_dir,r)   = nan;
        end
    end
end

% -------------------------------
% Pack results
% -------------------------------
res.rate_tot = ss_rate_tot_sac;   % [500 × 8 × 4] (time × direction × cond)
res.vm_tot   = ss_vm_tot_sac;     % [500 × 8 × 4]
res.conds    = conds;             % save condition matrix for reference

end
