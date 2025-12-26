function JDM_build_CS_rate_dir_corr_sac(data_path)
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

for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess = numel(current_sess_list);

    cs_on_data_sess  = cell(num_sess,1);
    cs_cell_ids_sess = cell(num_sess,1);
    cs_rate_sess     = cell(num_sess,1);


    current_sess_cr_ma = data_behav.sess_cr_ma{counter_animal};
    current_sess_cr_ms = data_behav.sess_cr_ms{counter_animal};
    %parfor
    parfor counter_sess = 1:num_sess
        current_sess = current_sess_list{counter_sess};
        disp(['     ' current_sess]);
        current_cr_ma = current_sess_cr_ma{counter_sess};
        current_cr_ms = current_sess_cr_ms{counter_sess};
        % load eye_traces
       
        units_info = MAF_extract_cell_metadata(current_path,current_sess);
        cell_types = units_info.cell_type;
        cell_ids   = units_info.cell_list;
        cs_cell_ids_ = cell_ids(ismember(cell_types,{'CS','PC'}));

        [cs_on_data_,cs_rate_] = par_loop_cs_on_data...
            (current_path,current_sess,cs_cell_ids_,...
            current_cr_ma,current_cr_ms);
        cs_cell_ids_sess{counter_sess} = cs_cell_ids_;
        cs_on_data_sess{counter_sess}  = cs_on_data_;
        cs_rate_sess{counter_sess}     = cs_rate_;
    end
    cs_on_data  = [cs_on_data; vertcat(cs_on_data_sess{:})];
    cs_cell_ids = [cs_cell_ids; vertcat(cs_cell_ids_sess{:})];
    if counter_animal == 1
        cs_rate = cat(1,cs_rate_sess{:});
    else
        cs_rate = cat(1,cs_rate, cat(1,cs_rate_sess{:}));
    end
end

file_name = 'CS_rate_dir_corr_sac_50-150';
save(fullfile(data_path,'population_data',sprintf('%s.mat',file_name)),...
   'cs_on_data','cs_cell_ids','cs_rate');
end
%% par_loop_cs_on_data
function [cs_on_data,cs_rate] = par_loop_cs_on_data(current_path,current_sess,cell_ids,...
  current_cr_ma,current_cr_ms)
  num_cs = numel(cell_ids);
  cs_on_data = cell(num_cs,1);
  cs_rate = nan(num_cs,8,4,2);

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

        res_rate = compute_cs_dir_tuned(current_path,...
            current_sess,current_cs,1:8,sac_crit);
  
        cs_rate(counter_cs,:,:,:) = res_rate.rate_tot;

    end
end
%% load_cs_data_tuned
function res = compute_cs_dir_tuned(current_path,current_sess,current_cell,order_ang,sac_crit)
MAF_params_funcs;

tags = 4;
win_vis = [50,150];
win_sac = [-70,30];

force_conds  = params.reward_tag4_conds;
n_dirs = 8;
n_conds = size(force_conds,1);

cs_rate_tot_sac = nan(1,n_dirs, n_conds);
cs_rate_tot_vis = nan(1, n_dirs, n_conds);

cell_data = MAF_load_cell(current_path,current_sess,current_cell);
neural_sac_vis =  JDM_combine_dataset_cs([cell_data.data_recordings],[cell_data.rec_info],{'visual'},win_vis,tags);
neural_sac_sac =  JDM_combine_dataset_cs([cell_data.data_recordings],[cell_data.rec_info],{'onset'},win_sac,tags);

% vis
neur_data_vis.CS        = neural_sac_vis.visual.CS;
sac_data_vis.sac        = neural_sac_vis.visual.sac;
sac_data_vis.fix        = neural_sac_vis.visual.fix;
sac_data_vis.sac.has_ss = neural_sac_vis.visual.has_ss;
sac_data_vis.sac.has_cs = neural_sac_vis.visual.has_cs;
sac_data_dir_vis        = MAF_bin_by_ang(neur_data_vis,sac_data_vis,'visual');

% sac
neur_data_sac.CS        = neural_sac_sac.onset.CS;
sac_data_sac.sac        = neural_sac_sac.onset.sac;
sac_data_sac.fix        = neural_sac_sac.onset.fix;
sac_data_sac.sac.has_cs = neural_sac_sac.onset.has_cs;
sac_data_dir_sac   = MAF_bin_by_ang(neur_data_sac,sac_data_sac,'onset');


 for counter_dir = 1:n_dirs  
    for cond_idx = 1:n_conds
       task_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.task_cond;
       task_tag = sac_data_dir_vis(order_ang(counter_dir)).sac.tag;
       tgt_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.tgt_cond;
       jmp_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.jump_cond;
       rwd_cond = sac_data_dir_vis(order_ang(counter_dir)).sac.rew_cond;

       istuned_vis = ismember([task_cond' task_tag' tgt_cond' jmp_cond' rwd_cond'],...
        force_conds(cond_idx,:),'rows');
       has_cs_vis = sac_data_dir_vis(order_ang(counter_dir)).sac.has_cs;
        idx_tuned_vis = istuned_vis & has_cs_vis';
        current_amp = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_amp(idx_tuned_vis);
        current_vel = sac_data_dir_vis(order_ang(counter_dir)).sac.eye_vm_max(idx_tuned_vis);
        current_dec = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_offset(idx_tuned_vis) -...
            sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(idx_tuned_vis))*1e3;
        current_acc = (sac_data_dir_vis(order_ang(counter_dir)).sac.time_vmax(idx_tuned_vis) -...
            sac_data_dir_vis(order_ang(counter_dir)).sac.time_onset(idx_tuned_vis))*1e3;
        current_raster_vis = sac_data_dir_vis(order_ang(counter_dir)).CS(:,idx_tuned_vis);
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
        cs_rate_tot_vis(:, counter_dir, cond_idx) = mean(ESN_smooth(mean(current_raster_vis,2))*1e3,'omitnan');

       % sac
       task_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.task_cond;
       tgt_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.tgt_cond;
       task_tag = sac_data_dir_sac(order_ang(counter_dir)).sac.tag;
       jmp_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.jump_cond;
       rwd_cond = sac_data_dir_sac(order_ang(counter_dir)).sac.rew_cond;

       istuned_sac = ismember([task_cond' task_tag' tgt_cond' jmp_cond' rwd_cond'],...
        force_conds(cond_idx,:),'rows');
      

        has_cs_sac = sac_data_dir_sac(order_ang(counter_dir)).sac.has_cs;
        idx_tuned_sac = istuned_sac & has_cs_sac';
        current_amp = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_amp(idx_tuned_sac);
        current_vel = sac_data_dir_sac(order_ang(counter_dir)).sac.eye_vm_max(idx_tuned_sac);
        current_dec = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_offset(idx_tuned_sac) -...
            sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(idx_tuned_sac))*1e3;
        current_acc = (sac_data_dir_sac(order_ang(counter_dir)).sac.time_vmax(idx_tuned_sac) -...
            sac_data_dir_sac(order_ang(counter_dir)).sac.time_onset(idx_tuned_sac))*1e3;
        current_raster_sac = sac_data_dir_sac(order_ang(counter_dir)).CS(:,idx_tuned_sac);
    
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
        cs_rate_tot_sac(:, counter_dir, cond_idx) = mean(ESN_smooth(mean(current_raster_sac,2))*1e3,'omitnan');    

    end
 end
res.rate_tot = cat(4, cs_rate_tot_vis, cs_rate_tot_sac); % [1 x 8 x 2 x 2] high\low

end
%% JDM_combine_dataset_cs
function neural_sac_data = JDM_combine_dataset_cs(data_recs, rec_info, alignment_vec, time_span, tags_of_interest)
    
    MAF_params_funcs;
    
    if nargin == 4
        tags_of_interest = 1:numel(params.sac.tag_name_list);
    end
    
    recs = find(rec_info.rec_flag);
    
    for counter_alignment = 1:numel(alignment_vec)
        current_alignment = alignment_vec{counter_alignment};
%         mod_span = [-time_span(1),time_span(2)] + [1,0];
        neural_sac_data.([current_alignment '_time']) = time_span(1)+1:time_span(2);

        num_rec = numel(data_recs);
    
        for counter_rec = 1:num_rec
            current_rec = data_recs(counter_rec);
            Neural_Data = current_rec.Neural_Data;
            neural_sac_data_ = JDM_buildSacData(current_rec.eye, Neural_Data, current_alignment, time_span, tags_of_interest);
            
            % Use filtered saccade times from neural_sac_data_ instead of raw data
            filtered_sac_times = neural_sac_data_.sac.(['time_' current_alignment]); % Key fix: Use filtered saccade times
            
            if counter_rec == 1
                neural_sac_data.(current_alignment) = neural_sac_data_;
                neural_sac_data.(current_alignment).rec_num = recs(1) .* ones(size(neural_sac_data_.sac.tag));
                
                % Initialize has_ss and has_cs using filtered saccades
                has_ss = ~isempty(Neural_Data.SS_time) .* ones(size(neural_sac_data_.sac.tag));
                if isfield(Neural_Data, 'SS_invalids') && ~isempty(Neural_Data.SS_invalids)
                    starts = Neural_Data.SS_invalids(:,1);
                    ends = Neural_Data.SS_invalids(:,2);
        
                    in_interval = false(size(filtered_sac_times));
                    for i = 1:size(Neural_Data.SS_invalids, 1)
                          in_interval = in_interval | ...
                          (filtered_sac_times >= starts(i) & filtered_sac_times <= ends(i));
                    end
                    has_ss(in_interval) = 0;
                end
                neural_sac_data.(current_alignment).has_ss = has_ss;
                
                % Repeat for has_cs
                has_cs = ~isempty(Neural_Data.CS_time) .* ones(size(neural_sac_data_.sac.tag));
                if isfield(Neural_Data, 'CS_invalids') && ~isempty(Neural_Data.CS_invalids)
                    starts = Neural_Data.CS_invalids(:,1);
                    ends = Neural_Data.CS_invalids(:,2);
        
                    in_interval = false(size(filtered_sac_times));
                    for i = 1:size(Neural_Data.CS_invalids, 1)
                          in_interval = in_interval | ...
                          (filtered_sac_times >= starts(i) & filtered_sac_times <= ends(i));
                    end
                    has_cs(in_interval) = 0;
                end
                neural_sac_data.(current_alignment).has_cs = has_cs;
                
            else
                % Concatenate data fields
                neural_sac_data.(current_alignment).SS = [neural_sac_data.(current_alignment).SS, neural_sac_data_.SS];
                neural_sac_data.(current_alignment).CS = [neural_sac_data.(current_alignment).CS, neural_sac_data_.CS];
                neural_sac_data.(current_alignment).sac = MAF_CatStructFields(...
                    neural_sac_data.(current_alignment).sac, neural_sac_data_.sac, 2);
                neural_sac_data.(current_alignment).fix = MAF_CatStructFields(...
                    neural_sac_data.(current_alignment).fix, neural_sac_data_.fix, 2);
                neural_sac_data.(current_alignment).rec_num = [neural_sac_data.(current_alignment).rec_num, recs(counter_rec) .* ones(size(neural_sac_data_.sac.tag))];
                
                % Compute has_ss for current recording using filtered saccades
                current_has_ss = ~isempty(Neural_Data.SS_time) .* ones(size(neural_sac_data_.sac.tag));
                if isfield(Neural_Data, 'SS_invalids') && ~isempty(Neural_Data.SS_invalids)
                    starts = Neural_Data.SS_invalids(:,1);
                    ends = Neural_Data.SS_invalids(:,2);
                    filtered_sac_times = filtered_sac_times(:);
                    % Vectorized comparison using implicit expansion (Mx1 vs 1xN)
                    in_interval = (filtered_sac_times >= starts') & (filtered_sac_times <= ends');
                    ss_invalid = any(in_interval, 2);
                    current_has_ss(ss_invalid) = 0;
                end
                neural_sac_data.(current_alignment).has_ss = [neural_sac_data.(current_alignment).has_ss, current_has_ss];
                
                % Repeat for has_cs
                current_has_cs = ~isempty(Neural_Data.CS_time) .* ones(size(neural_sac_data_.sac.tag));
                if isfield(Neural_Data, 'CS_invalids') && ~isempty(Neural_Data.CS_invalids)
                    starts = Neural_Data.CS_invalids(:,1);
                    ends = Neural_Data.CS_invalids(:,2);
                    filtered_sac_times = filtered_sac_times(:);
                    in_interval = (filtered_sac_times >= starts') & (filtered_sac_times <= ends');
                    cs_invalid = any(in_interval, 2);
                    current_has_cs(cs_invalid) = 0;
                end
                neural_sac_data.(current_alignment).has_cs = [neural_sac_data.(current_alignment).has_cs, current_has_cs];
            end
        end
    end
end

%% JDM_buildSacData
function neural_sac_data = JDM_buildSacData(sac_data,spikeData,alignment,...
    window,tags_of_interest)

if nargin == 5
    ind_tag = ismember(sac_data.sac.tag,tags_of_interest);
elseif nargin == 4
    ind_tag = ones('like',sac_data.sac.tag);
end

MAF_params_funcs;

ss_time = spikeData.SS_time;
cs_time = spikeData.CS_time;

event_time  = sac_data.sac.(['time_',alignment]);
ind_valid = not(isnan(event_time)) & ind_tag;

event_time = event_time(ind_valid);
neural_sac_data.SS = JDM_sp_time_to_raster(ss_time,event_time,window);
neural_sac_data.CS = JDM_sp_time_to_raster(cs_time,event_time,window);
neural_sac_data.sac = structfun(@(x) x(:,ind_valid),sac_data.sac,...
    'UniformOutput',false);
neural_sac_data.fix = structfun(@(x) x(:,ind_valid),sac_data.fix,...
    'UniformOutput',false);

end
%% JDM_sp_time_to_raster
function raster = JDM_sp_time_to_raster(sp_time, event_time, window)
sp_time = sp_time * 1e3;
event_time = event_time * 1e3;

num_event = length(event_time);

window_size = window(2) - window(1) + 1;
raster_ind = cell(num_event,1);
raster = zeros(window_size,num_event,'logical');
if isempty(sp_time)
    return;
end
parfor counter_event = 1:num_event
    current_event = event_time(counter_event);
    spike_time_rel = sp_time - current_event;
    spike_time_rel(spike_time_rel<window(1)) = [];
    spike_time_rel(spike_time_rel>window(2)) = [];
    raster_ind{counter_event} = ...
        int32(round(spike_time_rel - window(1))+1)+...
        (counter_event-1)*window_size;
end
raster_ind(cell2mat(cellfun(@isempty,...
    raster_ind,'UniformOutput',false))) = [];
raster_ind = cell2mat(raster_ind);
raster(raster_ind) = 1;

end
