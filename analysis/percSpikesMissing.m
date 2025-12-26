function [percentMissing_gaussian, percentMissing_symmetric, ksTest_pValue,...
    ampliBin_gaussian, spikeCountsPerAmpliBin_gaussian, gaussianFit_cutoff,...
    spike_times_presence] = percSpikesMissing(theseAmplitudes, theseSpikeTimes, timeChunks,neuron_label,unit_id, plotThis)
% Estimate the amount of spikes missing (below the detection threshold)
% by fitting a gaussian to the amplitude distribution for each timeChunk
% defined in timeChunks
% ------

% Adapted from Bombcell toolbox and modified.

warning off; % prevent warnings from optifit function about possible local minimums (time-consumming)
% initialize variables 
percentMissing_gaussian = nan(numel(timeChunks)-1, 1);
percentMissing_symmetric = nan(numel(timeChunks)-1, 1);
ksTest_pValue = nan(numel(timeChunks)-1, 1);
mean_firing   = nan(numel(timeChunks)-1, 1);

if plotThis
    fig=figure('Color', 'none');
    sgtitle(['Neuron type: ' neuron_label '       ' 'unit_id: ' num2str(unit_id)])
end


for iTimeChunk = 1:numel(timeChunks) - 1

    % get firing rate
    spk_rate = compute_firingRate (theseSpikeTimes(theseSpikeTimes >= timeChunks(iTimeChunk) & theseSpikeTimes < timeChunks(iTimeChunk+1)),0.1);
    mean_firing(iTimeChunk) = mean(spk_rate);

    % amplitude histogram
    nBins = 50;
    % get amplitude histogram 
    [spikeCountsPerAmpliBin, bins] = histcounts(theseAmplitudes(theseSpikeTimes >= timeChunks(iTimeChunk) & ...
        theseSpikeTimes < timeChunks(iTimeChunk+1)), nBins);
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

        p0 = [max(spikeCountsPerAmpliBin_gaussian), mode_seed, 2 * nanstd(theseAmplitudes), prctile(theseAmplitudes, 1)]; % seed


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

    roundedP = num2str(round(percentMissing_gaussian(iTimeChunk), 1));
    if plotThis
        subplot(5, numel(timeChunks)-1, numel(timeChunks)-1+iTimeChunk)
        if sum(spikeCountsPerAmpliBin_gaussian) > 0
            hold on;
            plot(gaussianFit_cutoff, ampliBin_gaussian, 'r');
            barh(ampliBin_gaussian, spikeCountsPerAmpliBin_gaussian, 'FaceColor', [0, 0.35, 0.71], 'EdgeColor', [0, 0.35, 0.71]);

            if iTimeChunk == 1
                xlabel('count')
                ylabel('amplitude')
            end
            if iTimeChunk > 1
                set(gca,'Yticklabel',[]) 
            end
            if iTimeChunk == 1
                title({['% missing spikes: ', newline, roundedP]}, 'Color', [0.7, 0.7, 0.7])
            else
                title(roundedP, 'Color', [0.7, 0.7, 0.7]);
            end
        hold off
        end
             
    end
end

warning on;
if plotThis 
    if exist('prettify_plot', 'file')
        prettify_plot('FigureColor', 'w', 'XLimits', 'keep', 'YLimits', 'all')
    else
        warning('https://github.com/Julie-Fabre/prettify-matlab repo missing - download it and add it to your matlab path to make plots pretty')
        makepretty('none')
    end

    subplot(4, numel(timeChunks)-1, [1:numel(timeChunks) - 1])
    scatter(theseSpikeTimes, theseAmplitudes, 4, [0, 0.35, 0.71], 'filled');
    hold on;
    % chunk lines
    ylims = ylim;
    for iTimeChunk = 1:length(timeChunks)
        line([timeChunks(iTimeChunk), timeChunks(iTimeChunk)], ...
            [ylims(1), ylims(2)], 'Color', [0.7, 0.7, 0.7])
    end
    xlabel('time (s)')
    ylabel(['amplitude scaling', newline, 'factor'])
    if exist('prettify_plot', 'file')
        prettify_plot('FigureColor', 'w', 'XLimits', 'keep', 'YLimits', 'keep')
    else
        warning('https://github.com/Julie-Fabre/prettify-matlab repo missing - download it and add it to your matlab path to make plots pretty')
    end
hold off;
    
   
    bin_size = 10; 
    [spike_rate,time_bins] = JDM_compute_firingRate(theseSpikeTimes, bin_size);
    
    % Plot histogram of spike_rate as a function of edges
    subplot(4, 1, 3)
    bar(time_bins, spike_rate, 'histc'); % 'histc' ensures bars are centered on the edges
    xlabel('Time (s)');
    ylabel('Spike Rate (spikes/s)'); 
    grid on;
    box off
    prettify_plot
   
end



% chunks with the least missing spikes
percentMissing_gaussian = percentMissing_gaussian(~isnan(percentMissing_gaussian));
avg_missing = mean(percentMissing_gaussian);

ind_miss      = percentMissing_gaussian < (mean(percentMissing_gaussian) + avg_missing);
ref_miss_mean = mean(mean_firing(ind_miss),'omitnan');

if ref_miss_mean >=10
    upper_bound = ref_miss_mean + ref_miss_mean*0.3;
    lower_bound = ref_miss_mean - ref_miss_mean*0.3;
else
    upper_bound = ref_miss_mean + ref_miss_mean*0.5;
    lower_bound = ref_miss_mean - ref_miss_mean*0.5;
end

if ~isempty(avg_missing)
    
    spike_times_presence = zeros(size(theseSpikeTimes));

    for iChunk = 1:length(timeChunks) - 1
    
        if mean_firing(iChunk) >= lower_bound && mean_firing(iChunk) <= upper_bound
            % Find start index (first spike in the chunk)
            startIdx = find(theseSpikeTimes >= timeChunks(iChunk), 1, 'first');
            % Find stop index (last spike in the chunk)
            stopIdx = find(theseSpikeTimes < timeChunks(iChunk + 1), 1, 'last');
            spike_times_presence(startIdx:stopIdx+1) = 1;
        end
    end
else
    spike_times_presence = zeros(size(theseSpikeTimes));
end

subplot(4,1,4)  
plot(theseSpikeTimes,spike_times_presence)
prettify_plot


fig.WindowState = 'maximized';
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
        warning('spike_times is empty. Returning empty spike rate.');
        spike_rate = [];
        return;
    end

    % Ensure spike_times is sorted
    spike_times = sort(spike_times);

    % Calculate the total duration of the recording
    total_duration = spike_times(end) - spike_times(1);

    % Handle edge case where total_duration is less than the time window
    if total_duration < time_window
        warning('Total duration is less than the time window. Returning empty spike rate.');
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
