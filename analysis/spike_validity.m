
clear
kilosort_output_path = 'Z:\Ephys\data_132F\2024-03\2024-03-11\cat_recs\analyzed_data\sorted_data';
sp = loadKSdir(kilosort_output_path);

%%
% Change this to the ID of the unit you want to plot
id_n            = sp.clu_info.id;
[~, idx_group]  = ismember(sp.clu_info.group,'good','row');
idx_goods       = id_n(idx_group == 1);
neurontype      = sp.clu_info.neurontype;
labels_cells    = neurontype(idx_group==1,:);
% wave_forms      = squeeze(sp.waveform.mean_wf(idx_group==1,:,:));
% cl_waveform     = sp.waveform;
% wave_chs        = squeeze(sp.waveform.ch_num(idx_group==1,:,:));


chunk_duration  = 5*60;  % Duration of each chunk in seconds
responses       = nan(numel(idx_goods),1);
id_cell         = nan(numel(idx_goods),1);

check_length = zeros(numel(idx_goods),2);

for ii = 1:numel(idx_goods)

    unit_id = idx_goods(ii);

    % Spike times of the unit
    unit_spike_inds = double(sp.sp_info.st(sp.sp_info.clu == unit_id)); 

    unit_spike_times_sec = unit_spike_inds/3e4;
    % Amplitudes of the unit's spikes
    unit_amplitudes      = double(sp.sp_info.amp(sp.sp_info.clu == unit_id));  
    
    id_label = deblank(labels_cells(ii,:));   % Neurontye of the cell
  
    % Get the total recording duration in seconds
    timeChunks = [min(unit_spike_times_sec):chunk_duration:max(unit_spike_times_sec), max(unit_spike_times_sec)];
    
    [perc_missing, ~, ~, ~, ~, ~,spike_times_presence] = ...
        percSpikesMissing(unit_amplitudes, unit_spike_times_sec, timeChunks,id_label,unit_id,1);

    [st_ind,end_ind] = JDM_spike_validity(unit_amplitudes, unit_spike_inds);

    % Loop until a valid key ('k' or 'l') is pressed
    validInput = false;
    while ~validInput
        % Ask for input and capture the first character
        keyPressed = input('Press "k" or "l" to continue: ', 's');
        
        % Check if the input is either 'k' or 'l'
        if ismember(keyPressed, ['k', 'l'])
            validInput = true;
            responses(ii,1) = keyPressed; % Save response
            disp(['You pressed: ', keyPressed]);
            id_cell(ii) = unit_id;
        else
            disp('Invalid input. Please press "k" or "l".');
        end
    end
    % Close the figure window
     close(gcf);

end 
disp('Curation of this session completed');
%%
figure;
subplot(411)
plot(rec1.ss_index), hold on
xline(rec1.ss_invalid_inds,'r','linewidth',2)

subplot(412)
plot(rec2.ss_index), hold on
xline(rec2.ss_invalid_inds,'r','linewidth',2)

subplot(413)
plot(rec3.ss_index), hold on
xline(rec3.ss_invalid_inds,'r','linewidth',2)

subplot(414)
plot(rec4.ss_index), hold on
xline(rec4.ss_invalid_inds,'r','linewidth',2)

