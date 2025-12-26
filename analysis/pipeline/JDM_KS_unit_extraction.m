%% function cambridge unit_extraction
function JDM_KS_unit_extraction(sess_path)

pat_SS = 'SS' + digitsPattern(1,2);
pat_CS = 'CS' + digitsPattern(1,2);
if ~strcmp(sess_path(end), filesep)
    sess_path = [sess_path filesep];
end

folders_ = strsplit(sess_path, filesep);
current_sess = folders_{end-1};
sess_name = regexprep(current_sess(3:end),'-','');
sess_meta_data = readtable([sess_path sess_name '.xls']);

rec_list = sess_meta_data.folder_name(logical...
    (sess_meta_data.ephys .* sess_meta_data.eye));

concat_recs_path = fullfile(sess_path, 'cat_recs');

 if isfolder(concat_recs_path)
    % Unit extraction for concatenated data (Cambridge)
    
    path_to_sorted = fullfile(sess_path, 'cat_recs', filesep, 'analyzed_data', filesep, 'sorted_data', filesep);

      % specify that the data sorted
    params_path = fullfile(path_to_sorted, 'params.py');
    if exist(params_path,'file') ~=2 
       error("The DATA has not been sorted, First run Kilosort!")
    end

    offset_dir = dir(fullfile(path_to_sorted, '*offsets.txt'));
    offsets = read_offsets_file(fullfile(offset_dir.folder, offset_dir.name));
    
    % Start and stop indices and times
    start_indices = offsets.smp_imap0;
    num_inds = round(diff([offsets.sec_imap0  offsets.fileTimeSecs]) * offsets.srate_imap0(1));
    if length(start_indices) ~= length(num_inds)
        error("Length of the offsets should be matched!")
    end
    stop_indices = offsets.smp_imap0 + num_inds - 1;
    onset_times = offsets.sec_imap0;
    offset_times = offsets.sec_imap0 + diff([offsets.sec_imap0, offsets.fileTimeSecs]);
    durations = offset_times - onset_times;
    % load the spikes indices
    sp = loadKSdir(path_to_sorted);


    num_rec = length(rec_list);
    for counter_rec = 1 : num_rec
        rec_name = rec_list{counter_rec};
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec)...
            ' Analyzing rec ', current_rec '\n']);

        units_dir = [sess_path, current_rec, filesep, ...
            'analyzed_data', filesep, 'units' filesep];
        if exist(units_dir, 'dir')
            rmdir(units_dir, 's');
        end
        path_to_data = [sess_path, current_rec, filesep, ...
            'raw_data', filesep];
    
        ch_file = dir([path_to_data,'*_CH1.continuous']);
        if isempty(ch_file)
            ch_file = dir([path_to_data,'*_1.continuous']);
            if isempty(ch_file)
                fprintf('\n No ephys!\n')
                continue;
            end
        end
        [~, ch_time, ~] = load_open_ephys_data(...
            [path_to_data,ch_file(1).name]);
        t_start = ch_time(1);

        %t_start = onset_times(1);
        ind_st  = start_indices(counter_rec);
        ind_end = stop_indices(counter_rec);
        n_samples = ind_end - ind_st + 1;
       
        ind_rec = sp.sp_info.st >= ind_st & sp.sp_info.st <= ind_end;
        sp_inds = sp.sp_info.st(ind_rec) - ind_st;
        sp_clust = sp.sp_info.clu(ind_rec);
   
        sample_rate = sp.sample_rate;
        cl_waveform = sp.waveform;
        ch_num = sp.clu_info.ch;
        id = sp.clu_info.id;
        [~, ind_good] = ismember(sp.clu_info.group, 'good', 'rows');
        good_ids = id(ind_good == 1);
        neurontype = sp.clu_info.neurontype;

        UNIT_NAME = cell(numel(good_ids), 1);
        cell_idx = cell(numel(good_ids), 1);
        neuron_type = cell(numel(good_ids), 1);

        for counter_gid = 1:numel(good_ids)
            current_id = good_ids(counter_gid);
            label = deblank(neurontype(id == current_id, :));
            if isnan(label)
                label = '';
            end

            if startsWith(label, pat_SS)
                ss_id = current_id;
                pc_num = str2double(erase(label, 'SS'));
                ind_cs = contains(string(neurontype), ['CS' num2str(pc_num)]) .* ind_good;
                % Handle cases of label overlap, i.e., CS1/CS10/CS11
                if sum(ind_cs) > 1
                    ind_cs = matches(deblank(string(neurontype)), string(['CS' num2str(pc_num)])) .* ind_good;
                end
                cs_id = id(ind_cs == 1);
                ch_id = cs_id;
                type = 'PC';

            elseif startsWith(label, pat_CS)
                continue;

            elseif strcmp(label, 'CS')
                ss_id = nan;
                cs_id = current_id;
                ch_id = cs_id;
                type = 'CS';

            else
                ss_id = current_id;
                cs_id = nan;
                ch_id = ss_id;
                type = label;
            end

            ss_index_int = sp_inds(sp_clust == ss_id);
            cs_index_int = sp_inds(sp_clust == cs_id);
            ch_ = ch_num(id == ch_id);
            time_ind = 1:n_samples;
            ss_index = ismember(time_ind, ss_index_int)';
            cs_index = ismember(time_ind, cs_index_int)';

       
           % end validity
            ss_index = resolve_auto_conflicts(ss_index, sample_rate, .5e-3);

            if sum(cs_index) > 0
                cs_index = resolve_auto_conflicts(cs_index, sample_rate, 5e-3);
                ss_index = resolve_cs_ss_conflicts(ss_index, cs_index, sample_rate, .5e-3, .5e-3);
            end

            waveform.ch_map.x = sp.xcoords;
            waveform.ch_map.y = sp.ycoords;
            waveform.ch_map.map = sp.ch_map;

            SS_rate = sum(ss_index) / (n_samples/sample_rate);
            CS_rate = sum(cs_index) / (n_samples/sample_rate);
            
            % Skip neurons with no spikes, but preserve alignment
            if SS_rate < 0.1 && CS_rate < 0.1
                fprintf('Neuron ID %d skipped: No SS or CS spikes\n', current_id);
               ss_invalid_inds = [];
               cs_invalid_inds = [];
               UNIT_NAME{counter_gid} = [num2str(ch_+1, '%.2i') '_' num2str(current_id, '%.3i')];
               cell_idx{counter_gid} = -1;
               neuron_type{counter_gid} = type;
               sess_name = regexprep(rec_name(3:end),'-','');
               unit_name = [sess_name '_' num2str(ch_+1,'%.2i') '_' ...
                             num2str(current_id,'%.3i')];
               out_dir = [units_dir unit_name filesep];
               mkdir(out_dir);
                 save([out_dir unit_name '_sorted_KSA'],'waveform','t_start','sample_rate',...
                'ss_invalid_inds','cs_invalid_inds','ss_index','cs_index','type','-v7.3');

                continue; 
            end
          
           %  add validity
           SS_id_total  = sp.sp_info.st(sp.sp_info.clu == ss_id);
           SS_amp_total = sp.sp_info.amp(sp.sp_info.clu == ss_id);
           [st_ss_ind,end_ss_inds] = JDM_spike_validity(SS_amp_total, SS_id_total);   
           ss_st  = SS_id_total(st_ss_ind) - ind_st;
           ss_end = SS_id_total(end_ss_inds) - ind_st;
           ss_invalid_inds = invalidity_intervals(ss_st,ss_end,ss_index_int);

  
           
           CS_id_total  = sp.sp_info.st(sp.sp_info.clu == cs_id);
           CS_amp_total = sp.sp_info.amp(sp.sp_info.clu == cs_id);
           [st_cs_ind,end_cs_inds] = JDM_spike_validity(CS_amp_total, CS_id_total);   
           cs_st  = CS_id_total(st_cs_ind) - ind_st;
           cs_end = CS_id_total(end_cs_inds) - ind_st;
           cs_invalid_inds = invalidity_intervals(cs_st,cs_end,cs_index_int);

            if sum(ss_index) > 0
                ss_wave = squeeze(cl_waveform.mean_wf(id == ss_id, :, :))';
                ss_wave_ch = squeeze(cl_waveform.ch_num(id == ss_id, :))';
                ss_peak = sp.sp_info.amp(sp_clust == ss_id);
            end

            if sum(cs_index) > 0
                cs_wave = squeeze(cl_waveform.mean_wf(id == cs_id, :, :))';
                cs_wave_ch = squeeze(cl_waveform.ch_num(id == cs_id, :))';
                cs_peak = sp.sp_info.amp(sp_clust == cs_id);
            end

            if sum(ss_index) > 0
                waveform.ss_wave    = ss_wave;
                waveform.ss_wave_ch = ss_wave_ch;
                waveform.ss_peak    = ss_peak;
            else
                waveform.ss_wave    = nan(size(cs_wave));
                waveform.ss_wave_ch = cs_wave_ch;
                waveform.ss_peak    = [];
            end

            if sum(cs_index) > 0
                waveform.cs_wave    = cs_wave;
                waveform.cs_wave_ch = cs_wave_ch;
                waveform.cs_peak    = cs_peak;
            else
                waveform.cs_wave    = nan(size(ss_wave));
                waveform.cs_wave_ch = ss_wave_ch;
                waveform.cs_peak    = [];
            end

            sess_name = regexprep(rec_name(3:end),'-','');
            unit_name = [sess_name '_' num2str(ch_+1,'%.2i') '_' ...
                num2str(current_id,'%.3i')];
            out_dir = [units_dir unit_name filesep];
            mkdir(out_dir);
            save([out_dir unit_name '_sorted_KSA'],'waveform','t_start','sample_rate',...
                'ss_invalid_inds','cs_invalid_inds','ss_index','cs_index','type','-v7.3');
            close all;
            MAF_plot_cell_summary([out_dir unit_name '_sorted_KSA']);

            UNIT_NAME{counter_gid} = [num2str(ch_+1, '%.2i') '_' num2str(current_id, '%.3i')];
            cell_idx{counter_gid} = current_id;
            neuron_type{counter_gid} = type;
        end

        emptyRows = all(cellfun(@isempty, UNIT_NAME), 2);
        meta_data_temp = [UNIT_NAME cell_idx neuron_type];
        meta_data_temp(emptyRows, :) = [];
        [~,idx_r] = sort(meta_data_temp(:,1));
        unitName = meta_data_temp(idx_r,1);
        cell_idx = meta_data_temp(idx_r,2);
        neuron_type = meta_data_temp(idx_r,3);
      
        fullFilePath = fullfile(sess_path, rec_name, sess_name);
        rec_meta_data = table(unitName, cell_idx, neuron_type);
        rec_meta_data.Properties.VariableNames = {'unit_name', 'cell', 'type'};
        writetable(rec_meta_data, fullFilePath, 'FileType', 'spreadsheet', 'WriteMode', 'overwritesheet');

    end

 else
    num_rec = length(rec_list);
    
    for counter_rec = 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec)...
            ' Analyzing rec ', current_rec '\n']);
    
        units_dir = [sess_path, current_rec, filesep, ...
            'analyzed_data', filesep, 'units' filesep];
        if exist(units_dir,'dir')
            rmdir(units_dir,'s');
        end
    
        path_to_sort = [sess_path, current_rec, filesep, ...
            'analyzed_data', filesep, 'sorted_data', filesep];
        path_to_data = [sess_path, current_rec, filesep, ...
            'raw_data', filesep];
    
        ch_file = dir([path_to_data,'*_CH1.continuous']);
        if isempty(ch_file)
            ch_file = dir([path_to_data,'*_1.continuous']);
            if isempty(ch_file)
                fprintf('\n No ephys!\n')
                continue;
            end
        end
        [~, ch_time, ~] = load_open_ephys_data(...
            [path_to_data,ch_file(1).name]);
        t_start = ch_time(1);
        n_samples = length(ch_time);
        sp = loadKSdir(path_to_sort);
        sp_inds = sp.sp_info.st;
        sp_clust = sp.sp_info.clu;
        sample_rate = sp.sample_rate;
        cl_waveform = sp.waveform;
    
        ch_num = sp.clu_info.ch;
        id = sp.clu_info.id;
        [~, ind_good] = ismember(sp.clu_info.group, 'good', 'rows');
        good_ids = id(ind_good == 1);
        neurontype = sp.clu_info.neurontype;
        rec_name = rec_list{counter_rec};
    
        for counter_gid = 1:numel(good_ids)
            current_id = good_ids(counter_gid);
            label = deblank(neurontype(id == current_id,:));
            if isnan(label)
                label = '';
            end
            if startsWith(label,pat_SS)
                ss_id = current_id;
                pc_num = str2double(erase(label,'SS'));
                ind_cs = contains(string(neurontype),...
                    ['CS' num2str(pc_num)]) .* ind_good;
                % handle cases of label overlap, i.e. CS1/CS10/CS11
                if sum(ind_cs) > 1
                    ind_cs = matches(deblank(string(neurontype)),...
                        string(['CS' num2str(pc_num)])) .* ind_good;
                end
                cs_id = id(ind_cs == 1);
                ch_id = cs_id;
                type = 'PC';
            elseif startsWith(label,pat_CS)
                continue;
            elseif strcmp(label,'CS')
                ss_id = nan;
                cs_id = current_id;
                ch_id = cs_id;
                type = 'CS';
            else
                ss_id = current_id;
                cs_id = nan;
                ch_id = ss_id;
                type = label;
            end
    
            ss_index_int = sp_inds(sp_clust == ss_id);
            cs_index_int = sp_inds(sp_clust == cs_id);
            ch_ = ch_num(id == ch_id);
    
            time_ind = 1:n_samples;
            ss_index = ismember(time_ind,ss_index_int)';
            cs_index = ismember(time_ind,cs_index_int)';
    
            ss_index = resolve_auto_conflicts(ss_index,sample_rate,.5e-3);
    
            if sum(cs_index) > 0
                cs_index = resolve_auto_conflicts(cs_index,sample_rate,5e-3);
                ss_index = resolve_cs_ss_conflicts(ss_index,cs_index,sample_rate,.5e-3,.5e-3);
            end
    
            waveform.ch_map.x = sp.xcoords;
            waveform.ch_map.y = sp.ycoords;
            waveform.ch_map.map = sp.ch_map;
    
            if sum(ss_index) > 0
                ss_wave = squeeze(cl_waveform.mean_wf(id == ss_id,:,:))';
                ss_wave_ch = squeeze(cl_waveform.ch_num(id == ss_id,:))';
                ss_peak = sp.sp_info.amp(sp_clust == ss_id);
            end
    
            if sum(cs_index) > 0
                cs_wave = squeeze(cl_waveform.mean_wf(id == cs_id,:,:))';
                cs_wave_ch = squeeze(cl_waveform.ch_num(id == cs_id,:))';
                cs_peak = sp.sp_info.amp(sp_clust == cs_id);
            end
    
            if sum(ss_index) > 0
                waveform.ss_wave    = ss_wave;
                waveform.ss_wave_ch = ss_wave_ch;
                waveform.ss_peak    = ss_peak;
            else
                waveform.ss_wave    = nan(size(cs_wave));
                waveform.ss_wave_ch = cs_wave_ch;
                waveform.ss_peak    = [];
            end
    
            if sum(cs_index) > 0
                waveform.cs_wave    = cs_wave;
                waveform.cs_wave_ch = cs_wave_ch;
                waveform.cs_peak    = cs_peak;
            else
                waveform.cs_wave    = nan(size(ss_wave));
                waveform.cs_wave_ch = ss_wave_ch;
                waveform.cs_peak    = [];
            end
    
            sess_name = regexprep(rec_name(3:end),'-','');
            unit_name = [sess_name '_' num2str(ch_+1,'%.2i') '_' ...
                num2str(current_id,'%.3i')];
            out_dir = [units_dir unit_name filesep];
            mkdir(out_dir);
            save([out_dir unit_name '_sorted_KSA'],'waveform','t_start',...
                'sample_rate','ss_index','cs_index','type','-v7.3');
            close all;
            MAF_plot_cell_summary([out_dir unit_name '_sorted_KSA']);
    
        end
    
     end
  end
end

%% function solve SS-SS Conflicts
function ss_index_new = resolve_auto_conflicts(ss_index,sample_rate,win_look_around)
window_len = floor(win_look_around * sample_rate);
ss_index_ = ss_index;
ss_index_int_ = find(ss_index_);
for counter_ss = 1:length(ss_index_int_)
    ss_index_local_ = ss_index_int_(counter_ss);
    if ss_index_local_ < window_len
        ss_index_(ss_index_local_) = 0;
        continue
    end
    if ss_index_local_ > (length(ss_index_) - window_len)
        ss_index_(ss_index_local_) = 0;
        continue
    end
    search_win_inds = ss_index_local_-window_len:ss_index_local_+window_len;
    ss_search_win_bool = ss_index_(search_win_inds);
    ss_search_win_int  = find(ss_search_win_bool);
    if length(ss_search_win_int) < 2
        continue
    end
    if length(ss_search_win_int) > 1
        valid_ind = ss_search_win_int(1);
        ss_search_win_bool = zeros(size(search_win_inds));
        ss_search_win_bool(valid_ind) = 1;
        ss_index_(search_win_inds) = ss_search_win_bool;
    end
end
ss_index_new = ss_index_;
end

%% function solve SS-CS Conflicts
function ss_index_new = resolve_cs_ss_conflicts(ss_index,cs_index,sample_rate, win_look_before, win_look_after)
window_len_back = floor(win_look_before * sample_rate);
window_len_front = floor(win_look_after * sample_rate);
cs_index_int = find(cs_index);
ss_index_ = ss_index;
for counter_cs = 1:length(cs_index_int)
    cs_index_local = cs_index_int(counter_cs);
    search_win_inds = cs_index_local-window_len_back:cs_index_local+window_len_front;
    ss_search_win_bool = ss_index_(search_win_inds);
    ss_search_win_int  = find(ss_search_win_bool);
    if ~isempty(ss_search_win_int)
        ss_ind_invalid = ss_search_win_int + cs_index_local - window_len_back;
        ss_index_(ss_ind_invalid-1) = 0;
    end
end
ss_index_new = ss_index_;
end
%% Read supercat offset
function offsets = read_offsets_file(filename)
    % Open the file for reading
    fid = fopen(filename, 'r');
    
    % Initialize the output structure
    offsets = struct('smp_imap0', [], 'sec_imap0', [], 'srate_imap0', []);
    
    % Read through each line of the file
    while ~feof(fid)
        % Read the current line
        line = fgetl(fid);
        
        % Parse based on the line prefix (smp_imap0, sec_imap0, or srate_imap0)
        if startsWith(line, 'smp_imap0:')
            % Remove 'smp_imap0:' and extract all the numeric values
            line = strrep(line, 'smp_imap0:', '');
            offsets.smp_imap0 = sscanf(line, '%f')';
        elseif startsWith(line, 'sec_imap0:')
            % Remove 'sec_imap0:' and extract all the numeric values
            line = strrep(line, 'sec_imap0:', '');
            offsets.sec_imap0 = sscanf(line, '%f')';
        elseif startsWith(line, 'srate_imap0:')
            % Remove 'srate_imap0:' and extract all the numeric values
            line = strrep(line, 'srate_imap0:', '');
            offsets.srate_imap0 = sscanf(line, '%f')';
        elseif startsWith(line, 'fileTimeSecs:')
            % Remove 'fileTimeSecs:' and extract all the numeric values
            line = strrep(line, 'fileTimeSecs:', '');
            offsets.fileTimeSecs = sscanf(line, '%f')';
        end
    end
    
    % Close the file
    fclose(fid);
end
%%
function invalidity_inds = invalidity_intervals(statrt_sp_inds,end_sp_inds,spike_indices)

if numel(statrt_sp_inds) ~= numel(end_sp_inds)
    error('Number of starts and ends times must be the same!')
end

if isempty(statrt_sp_inds) && isempty(end_sp_inds)
    invalidity_inds = [];
else
    invalidity_inds = nan(numel(statrt_sp_inds),2);
     for ii = 1:numel(statrt_sp_inds)
         in_rec_spikes = spike_indices(spike_indices >= statrt_sp_inds(ii) & spike_indices <= end_sp_inds(ii));
           if ~isempty(in_rec_spikes)
                invalidity_inds(ii,:) = [in_rec_spikes(1) in_rec_spikes(end)];
           end
     end
end

end