function plot_example_waveform(data_path)

load(fullfile(data_path,'population_data','fig1_waveforms.mat'),"cell_list_name","cells");

prefixes   = extractBefore(cell_list_name, '_');
categories = {'pc','mli1','mli2','mf'};

out_dir = fullfile(data_path,'population_figs','ephys');

for c = 1:numel(categories)
    idx = strcmp(prefixes, categories{c});
    if ~any(idx)
        continue; % skip if none
    end
    
    this_list = cell_list_name(idx);
    n_cells   = numel(this_list);

    ncols = ceil(sqrt(n_cells));
    nrows = ceil(n_cells / ncols);

    fig = figure;
 
    for i = 1:n_cells
        s(i) = subplot(nrows, ncols, i);
        Neural_Prop_tot = MAF_combineNeuralProp([cells.(this_list{i}).Neural_Prop]);
        JDM_plot_waveform(Neural_Prop_tot.waveform,3e4);
        title(this_list{i}, 'Interpreter','none','FontSize',9);
    end

    sgtitle(upper(categories{c}), 'FontWeight','bold');
    ESN_Beautify_Plot(fig, [ncols*3,nrows*4]);
    
    saveas(fig, fullfile(out_dir, sprintf('waveforms_%s.pdf', categories{c})));
end

end
%%
function JDM_plot_waveform(waveform,sample_rate,color_SS,color_CS)

if nargin < 3
    color_SS = 'blue';
    color_CS = 'red';
end

hold on
ss_wave = waveform.ss_wave;
cs_wave = waveform.cs_wave;

x_ = waveform.ch_map.x * 4;
y_ = waveform.ch_map.y * 100;

if not(isnan(ss_wave))
    n_ch = size(ss_wave,1);
    n_sig = length(ss_wave);
    wave_ch = waveform.ss_wave_ch;

    x = x_(wave_ch+1);
    y = y_(wave_ch+1);
    ch_num = waveform.ch_map.map(wave_ch+1);

    span_ind = (0:n_sig-1)/sample_rate;
    span_group_ = repmat([span_ind,nan],n_ch,1);
    span_group = reshape((span_group_*1e3+x)',1,n_ch*(n_sig+1));

    ss_wave_ = [ss_wave, nan(n_ch,1)];
    ss_wave_group = reshape((ss_wave_+y)',1,n_ch*(n_sig+1));
    plot(span_group, ss_wave_group, 'color', color_SS, 'linewidth', 1)
end

if not(isnan(cs_wave))
    n_ch = size(cs_wave,1);
    n_sig = length(cs_wave);
    wave_ch = waveform.cs_wave_ch;

    x = x_(wave_ch+1);
    y = y_(wave_ch+1);
    ch_num = waveform.ch_map.map(wave_ch+1);

    span_ind = (0:n_sig-1)/sample_rate;
    span_group_ = repmat([span_ind,nan],n_ch,1);
    span_group = reshape((span_group_*1e3+x)',1,n_ch*(n_sig+1));

    cs_wave_ = [cs_wave, nan(n_ch,1)];
    cs_wave_group = reshape((cs_wave_+y)',1,n_ch*(n_sig+1));
    plot(span_group, cs_wave_group, 'color', color_CS, 'linewidth', 1);
end

ch_map = arrayfun(@num2str, ch_num+1,'UniformOutput', 0);

text(x-3,y,ch_map)

axis off

plot([0,1],[0,0]+min(y)-100,'black','LineWidth',1);
plot([0,0],[0,300]+min(y)-100,'black','LineWidth',1);

end