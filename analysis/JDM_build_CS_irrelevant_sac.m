function JDM_build_CS_irrelevant_sac(data_path)
JDM_params_funcs
load(fullfile(data_path,'population_data','behave_data'));
user = 'JDM';
loc  = 'ctx';

path_data_monkey_sorted = params.path_data_monkey_sorted;
animal_list = params.animal_list;
session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);

num_animal = numel(session_list);

cs_on_data  = cell(0);
cs_cell_ids = cell(0);
cs_rate     = [];
cs_vm       = [];

for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess = numel(current_sess_list);

    cs_on_data_sess  = cell(num_sess,1);
    cs_cell_ids_sess = cell(num_sess,1);
    cs_rate_sess     = cell(num_sess,1);
    cs_vm_sess       = cell(num_sess,1);

    current_sess_cr_ma = data_behav.sess_cr_ma{counter_animal};
    current_sess_cr_ms = data_behav.sess_cr_ms{counter_animal};
    %parfor
    parfor counter_sess = 1:num_sess
        current_sess = current_sess_list{counter_sess};
        disp(['     ' current_sess]);
        current_cr_ma = current_sess_cr_ma{counter_sess};
        current_cr_ms = current_sess_cr_ms{counter_sess};
        % load eye_traces
        data_eye = MAF_load_eye_traces(current_path,current_sess);
        units_info = MAF_extract_cell_metadata(current_path,current_sess);
        cell_types = units_info.cell_type;
        cell_ids   = units_info.cell_list;
        cs_cell_ids_ = cell_ids(ismember(cell_types,{'CS','PC'}));

        [cs_on_data_,~,cs_rate_,cs_vm_] = par_loop_cs_on_data...
            (current_path,current_sess,data_eye,cs_cell_ids_,...
            current_cr_ma,current_cr_ms);
        cs_cell_ids_sess{counter_sess} = cs_cell_ids_;
        cs_on_data_sess{counter_sess}  = cs_on_data_;
        cs_rate_sess{counter_sess}     = cs_rate_;
        cs_vm_sess{counter_sess}       = cs_vm_;
    end
    cs_on_data  = [cs_on_data; vertcat(cs_on_data_sess{:})];
    cs_cell_ids = [cs_cell_ids; vertcat(cs_cell_ids_sess{:})];
    if counter_animal == 1
        cs_rate = cat(1,cs_rate_sess{:});
        cs_vm   = cat(1,cs_vm_sess{:});
    else
        cs_rate = cat(1,cs_rate, cat(1,cs_rate_sess{:}));
        cs_vm   = cat(1,cs_vm, cat(1,cs_vm_sess{:}));
    end
end
file_name = 'cs_on_rate_irrelevant_onset_sac';
save(fullfile(data_path,'population_data',sprintf('%s.mat',file_name)),...
   'cs_on_data','cs_cell_ids','cs_rate','cs_vm','-v7.3');
end

%% par_loop_cs_on_data
function [cs_on_data,perturb,cs_rate,cs_vm] = par_loop_cs_on_data...
    (current_path,current_sess,data_eye,cell_ids,...
    current_cr_ma,current_cr_ms)
    num_cs = numel(cell_ids);
    cs_on_data = cell(num_cs,1);
    perturb = cell(num_cs,1);

    cs_rate = nan(num_cs,500,8);
    cs_vm  = nan(num_cs,500,8);
 
    sac_crit.ma   = current_cr_ma;
    sac_crit.ms   = current_cr_ms;
   %parfor
   parfor counter_cs = 1:num_cs
        current_cs = cell_ids{counter_cs};
    
        data = MAF_load_cell(current_path,current_sess,current_cs);
        data_recordings = data.data_recordings;
        cs_on_data{counter_cs} = MAF_CS_on_analysis_fr([data_recordings.eye],...
            [data_recordings.Neural_Data],[1,4,6,7],[40,85],[-70,30]);

        if cs_on_data{counter_cs}.vis.rho_avg > cs_on_data{counter_cs}.sac.rho_avg
            cs_on_ang = cs_on_data{counter_cs}.vis.ang_avg*pi/180;
        else
            cs_on_ang = cs_on_data{counter_cs}.sac.ang_avg*pi/180;
        end

        res_rate = load_cs_data_tuned(current_path,...
            current_sess,current_cs,1:8,sac_crit,data_eye);
  
        cs_rate(counter_cs,:,:) = res_rate.rate_tot;
        cs_vm(counter_cs,:,:)  = res_rate.vm_tot;
    end
end
%% load_cs_data_tuned
function res = load_cs_data_tuned(current_path,current_sess,current_cell,order_ang,sac_crit,data_eye)
    JDM_params_funcs;
    tags         = 10;
    n_timepoints = 500;
    n_dirs       = 8;

    % Preallocate outputs
    cs_rate_tot_sac = nan(n_timepoints, n_dirs);
    cs_vm_tot_sac   = nan(n_timepoints, n_dirs);

    % Load data
    data       = MAF_load_cell(current_path,current_sess,current_cell);
    eye_traces = data_eye.eye_traces(data.rec_info.rec_flag);

    neural_sac_data = JDM_combine_dataset_sac([data.data_recordings],[data.rec_info],{'onset'},250,tags);
    eye_vm          = MAF_combine_vm_traces(eye_traces,{'onset'},250,tags);

    % Package for angular binning
    neur_data_sac.eye_vm = eye_vm.onset.eye_vm;
    neur_data_sac.CS     = neural_sac_data.onset.CS;

    sac_data_sac.sac        = neural_sac_data.onset.sac;
    sac_data_sac.fix        = neural_sac_data.onset.fix;
    sac_data_sac.sac.has_cs = neural_sac_data.onset.has_cs; % kept (not used to gate)

    sac_data_dir_sac = MAF_bin_by_ang(neur_data_sac, sac_data_sac, 'onset');

    % -------- per-direction processing --------
    for counter_dir = 1:n_dirs
        d = order_ang(counter_dir);

        % Trial-level vectors
        amp_all   = sac_data_dir_sac(d).sac.eye_amp(:);        % deg
        vel_all   = sac_data_dir_sac(d).sac.eye_vm_max(:);     % peak vel (deg/s)
        t_on      = sac_data_dir_sac(d).sac.time_onset(:);     % s
        t_vmax    = sac_data_dir_sac(d).sac.time_vmax(:);      % s
        t_off     = sac_data_dir_sac(d).sac.time_offset(:);    % s

        % Fixation validity
        fix_before = sac_data_dir_sac(d).fix.fix_validity_before(:);
        fix_after  = sac_data_dir_sac(d).fix.fix_validity_after(:);
        idx_fix    = fix_before & fix_after;

        % Base inclusion: finite, 1–10 deg amplitude, valid fixation
        idx_base = isfinite(amp_all) & amp_all >= 1 & amp_all <= 10 & idx_fix;

        % Index into current trials
        current_amp = amp_all(idx_base);
        current_vel = vel_all(idx_base);
        current_on  = t_on(idx_base);
        current_vm  = sac_data_dir_sac(d).eye_vm(:, idx_base);
        current_cs  = sac_data_dir_sac(d).CS(:,     idx_base);

        % Kinematic times (ms)
        current_t_vmax = t_vmax(idx_base);
        current_t_off  = t_off(idx_base);
        acc_ms = (current_t_vmax - current_on) * 1e3;
        dec_ms = (current_t_off  - current_t_vmax) * 1e3;

        % Polygon filters (main sequence & amp–acc/dec)
        [in_ms, on_ms] = inpolygon(log10(current_amp), log10(current_vel), ...
                                   sac_crit.ms(1,:), sac_crit.ms(2,:));

        acc_dec_ratio = dec_ms ./ acc_ms; % can be Inf/NaN if acc_ms<=0 (removed below)
        [in_ad, on_ad] = inpolygon(log10(current_amp), log10(acc_dec_ratio), ...
                                   sac_crit.ma(1,:), sac_crit.ma(2,:));

        % Remove: outside polygons, non-positive acc/dec, non-finite velocity
        rm = ~(in_ms | on_ms) | ~(in_ad | on_ad) | (acc_ms <= 0) | (dec_ms <= 0) | ~isfinite(current_vel);

        current_cs(:, rm) = [];
        current_vm(:, rm) = [];

        if isempty(current_cs)
            cs_rate_tot_sac(:, counter_dir) = nan(n_timepoints,1);
            cs_vm_tot_sac(:,   counter_dir) = nan(n_timepoints,1);
        else
            cs_rate_tot_sac(:, counter_dir) = ESN_smooth(mean(current_cs, 2)) * 1e3; % Hz
            cs_vm_tot_sac(:,   counter_dir) = mean(current_vm, 2);
        end
    end

    % Outputs
    res.rate_tot = cs_rate_tot_sac;   % [500 × 8]
    res.vm_tot   = cs_vm_tot_sac;     % [500 × 8]
end
