function [startIndices, endIndices] = JDM_spike_validity(spikes_amp, spike_indices,params)
% Estimate the amount of spikes missing (below the detection threshold)
% by fitting a gaussian to the amplitude distribution and spike rate 
% for each timeChunk
% compute aplitude of missing spikes adapted from Bombcell toolbox and
% modified.(https://github.com/Julie-Fabre/bombcell)
% -------------------------------------------------------------------------

if nargin  < 3
    params.chunkLength = 10;
    params.winSize = 0.1;
    params.LowThr = 0.8;
    params.highThr = 0.5;
    params.missingAmpThr= 25;
end

% convert spike indices to second
spike_times = double(spike_indices) / 3e4;
spikes_amp  = double(spikes_amp);

chunk_duration  = params.chunkLength*60;  % Duration of each chunk in seconds(minutes*60)

% Get the total recording duration in seconds
timeChunks = [min(spike_times):chunk_duration:max(spike_times), max(spike_times)];
percentMissing_gaussian = nan(numel(timeChunks)-1, 1);
percentMissing_symmetric = nan(numel(timeChunks)-1, 1);
mean_firing   = nan(numel(timeChunks)-1, 1);


for iTimeChunk = 1:numel(timeChunks) - 1

    % get firing rate
    spk_rate = compute_firingRate (spike_times(spike_times >= timeChunks(iTimeChunk) & spike_times < timeChunks(iTimeChunk+1)),params.winSize);
    mean_firing(iTimeChunk) = mean(spk_rate);
    
    % amplitude histogram
    nBins = 50;
    % get amplitude histogram 
    [spikeCountsPerAmpliBin, bins] = histcounts(spikes_amp(spike_times >= timeChunks(iTimeChunk) & ...
        spike_times < timeChunks(iTimeChunk+1)), nBins);
    if sum(spikeCountsPerAmpliBin) > 5 % if there is at least 5 spikes in this time bin 
        
        maxAmpli_val = find(spikeCountsPerAmpliBin == max(spikeCountsPerAmpliBin));
        mode_seed = bins(maxAmpli_val); %guess mode - this is only accurate if mode is present in histogram
        if numel(mode_seed) > 1 %if two or more modes, take the mean
            mode_seed = mean(mode_seed);
        end
        bin_steps = diff(bins(1:2)); %size of a bin

        %% symmetric 
        % mirror the unit's amplitudes 
        % median smooth ti get more accurate peak ]
        spikeCountsPerAmpliBin_smooth = smoothdata(spikeCountsPerAmpliBin, 'movmedian', 5);
        maxAmpli_val_smooth = find(spikeCountsPerAmpliBin_smooth == max(spikeCountsPerAmpliBin_smooth),1);
       
        surrogate_amplitudes = [spikeCountsPerAmpliBin_smooth(end:-1:maxAmpli_val_smooth), fliplr(spikeCountsPerAmpliBin_smooth(end:-1:maxAmpli_val_smooth+1))];
        surrogate_bins = [fliplr(bins(maxAmpli_val_smooth)-bin_steps:-bin_steps:bins(maxAmpli_val_smooth) - bin_steps*floor(size(surrogate_amplitudes,2)/2)),...
            bins(maxAmpli_val:end)];
        
        % remove any bins where the amplitudes would be < 0, that doesn't make any sense 
        surrogate_amplitudes(surrogate_bins<0) = [];
        surrogate_bins(surrogate_bins<0) = [];
        surrogate_area = sum(surrogate_amplitudes)*bin_steps;
        
        % estimate the percentage of missing spikes 
        pMissing = (surrogate_area - sum(spikeCountsPerAmpliBin)*bin_steps)/surrogate_area * 100;
        if pMissing < 0 % if inferior, distribution isn't symmetric
            pMissing = 0;
        end
        percentMissing_symmetric(iTimeChunk) = pMissing; 

        %% evaluate gaussian distribution from symmetric 'surrogate' data 
        [~,ksTest_pValue(iTimeChunk)] = kstest(surrogate_amplitudes);

        %% gaussian
        ampliBin_gaussian = bins(1:end-1) + bin_steps / 2;
        next_low_bin = ampliBin_gaussian(1) - bin_steps;
        add_points = 0:bin_steps:next_low_bin; % add points so amplitude values starts at 0
        ampliBin_gaussian = [add_points, ampliBin_gaussian];
        spikeCountsPerAmpliBin_gaussian = [zeros(size(add_points, 2), 1)', spikeCountsPerAmpliBin];

        p0 = [max(spikeCountsPerAmpliBin_gaussian), mode_seed, 2 * nanstd(spikes_amp), prctile(spikes_amp, 1)]; % seed


        f = @(x, xdata)gaussian_cut(x, xdata); % get anonymous function handle

        options = optimoptions('lsqcurvefit', 'OptimalityTolerance', 1e-32, 'FunctionTolerance', 1e-32, 'Display', 'off'); %,'StepTolerance', 1e-20,...
        %'MaxFunctionEvaluations', 5000);%'MaxFunctionEvaluations', 10000, 'MaxIterations', 1000);
        lb = [];
        ub = [];
        fitOutput = lsqcurvefit(f, p0, ampliBin_gaussian, spikeCountsPerAmpliBin_gaussian, lb, ub, options); %QQ need to look into local minimum error that sometimes happens

        %norm area calculated by fit parameters

        gaussianFit_cutoff = JF_gaussian_cut(ampliBin_gaussian, fitOutput(1), fitOutput(2), fitOutput(3), fitOutput(4));
        %    n_fit_no_cut = JF_gaussian_cut(bin_centers, fitOutput(1), fitOutput(2), fitOutput(3), 0);
        norm_area_ndtr = normcdf((fitOutput(2) - fitOutput(4))/fitOutput(3)); %ndtr((popt[1] - min_amplitude) /popt[2])
        percentMissing_gaussian(iTimeChunk) = 100 * (1 - norm_area_ndtr);
    else
        percentMissing_gaussian(iTimeChunk) = NaN;
        percentMissing_symmetric(iTimeChunk) = NaN;
        ksTest_pValue(iTimeChunk) = NaN;
        ampliBin_gaussian = NaN;
        spikeCountsPerAmpliBin_gaussian = NaN;
        gaussianFit_cutoff = NaN;
    end  
end

%% Set the threshold for firing rate and missing amps
 
temp_amp = percentMissing_gaussian;
temp_amp(temp_amp==100)= 0;

mn_missing =  mean(temp_amp,'omitnan');
if mn_missing > params.missingAmpThr
    temp_amp = temp_amp(temp_amp < mn_missing);
else
    temp_amp = temp_amp(temp_amp<params.missingAmpThr);
end

threshold = mean(temp_amp,'omitnan'); 
% Set values below the threshold to zero
ind_missing = (percentMissing_gaussian < threshold);

spike_times_presence = [];

if all(ind_missing) 
     spike_times_presence = spike_indices;
else
    mean_firing(isnan(mean_firing)) = 0;
    spike_count = compute_firingRate(spike_times, params.winSize);
    avg_ref_rate = mean(mean_firing(ind_missing),'omitnan');
    
    if  avg_ref_rate < 2 % threshold for neuron with low rate
        lower_bound = avg_ref_rate - (avg_ref_rate * params.LowThr);
    else % threshold for neuron with high rate
        lower_bound = avg_ref_rate - (avg_ref_rate * params.highThr);
    end
    
    upper_bound = max(spike_count);

    %% 
    %% 
     if ~isempty(ind_missing)
        
        for iChunk = 1:length(timeChunks) - 1
            % Find start index (first spike in the chunk)
            [~, startIdx] = min(abs(spike_times - timeChunks(iChunk)));
            
            % Find stop index (last spike in the chunk)
            [~, stopIdx] = min(abs(spike_times - timeChunks(iChunk + 1)));
            stopIdx = min(stopIdx, length(spike_times));  % Ensure stopIdx is within bounds
            
            % Check if mean firing rate is within bounds and set presence
            if mean_firing(iChunk) >= lower_bound && mean_firing(iChunk) <= upper_bound
                spike_times_presence = [spike_times_presence,startIdx:stopIdx];
            end
        end
       
    end


end
valids_spikes = false(numel(spike_indices), 1);
valids_spikes(spike_times_presence) = true;

diffVec = diff([1; valids_spikes; 1]); % Pad with 1 to handle edges
startIndices = find(diffVec == -1); % Start of 0 sequences
endIndices = find(diffVec == 1) - 1; % End of 0 sequences

end
%% SUB FUNCTIONS

function F = gaussian_cut(x, bin_centers)
    F = x(1) * exp(-(bin_centers - x(2)).^2/(2 * x(3).^2));
    F(bin_centers < x(4)) = 0;
end

%gaussian cut 

function g = JF_gaussian_cut(x, a, x0, sigma, xcut)
     g = a .* exp(-(x - x0) .^ 2 / (2 * sigma .^ 2));
     g(x < xcut) = 0;
end

% compute spike rate
function spike_rate = compute_firingRate(spike_times, time_window)
    % Computes firing rate in bins of specified time window

    % Check if spike_times is empty
    if isempty(spike_times)
        spike_rate = [];
        return;
    end
    % Ensure spike_times is sorted
    spike_times = sort(spike_times);
    % Calculate the total duration of the recording
    total_duration = spike_times(end) - spike_times(1);
    % Handle edge case where total_duration is less than the time window
    if total_duration < time_window
        spike_rate = [];
        return;
    end
    % Define the edges of the time bins
    edges = 0:time_window:total_duration;
    % Ensure edges has at least two elements
    if numel(edges) < 2
        error('Time window is too large for the duration of spike times.');
    end
    % Count spikes in each bin
    spike_counts = histcounts(spike_times - spike_times(1), edges);
    % Calculate the firing rate (spikes per second)
    spike_rate = spike_counts / time_window;
end

