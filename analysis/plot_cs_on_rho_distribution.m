function plot_cs_on_rho_distribution(data_path)
    JDM_params_funcs;
load(fullfile(data_path,'population_data','cs_on_rate_spe_rpe_prim_sac'));

% --------- vectors for visual & saccade ----------
cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_data, 'UniformOutput', false));
cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_data, 'UniformOutput', false));

cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, cs_on_data, 'UniformOutput', false));
cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, cs_on_data, 'UniformOutput', false));

% --------- Figure 1: scatter + histogram ----------
fig = figure('Position',[100 100 1200 500],'Color','w');

% (Left) individual vectors (angle vs rho)
subplot(1,2,1);
polarscatter(deg2rad(cs_on_ang_sac), cs_on_rho_sac, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.7);
ax1 = gca; ax1.ThetaDir='counterclockwise'; ax1.ThetaZeroLocation='right';
ax1.ThetaTick = 0:45:315; ax1.ThetaTickLabel = arrayfun(@(x) sprintf('%d°',x),0:45:315,'uni',false);
rho_max = prctile(cs_on_rho_sac,95); ax1.RLim = [0, rho_max*1.2];
ax1.FontSize=12; ax1.GridColor=[.7 .7 .7]; ax1.LineWidth=1.2;
title('Individual cs-on','FontSize',14,'FontWeight','bold');

% (Right) counts or proportions per direction bin
subplot(1,2,2);
nBins = 24; 
binEdges = linspace(0, 2*pi, nBins+1);
[~,~,binIdx] = histcounts(deg2rad(cs_on_ang_sac), binEdges);
valid = binIdx ~= 0;

% ----- choose representation -----
showProportion = false;  % set true to plot proportion instead of count

counts = accumarray(binIdx(valid), 1, [nBins 1], @sum, 0);
if showProportion
    binVals = counts / numel(cs_on_ang_sac);
    radialTitle = 'Proportion of cells';
else
    binVals = counts;
    radialTitle = 'Cell count';
end

% Polar histogram with explicit BinCounts = counts/proportions
polarhistogram('BinEdges', binEdges, 'BinCounts', binVals, ...
               'FaceColor', '#D95319', 'FaceAlpha', 0.85, 'EdgeColor', 'k');

ax2 = gca; 
ax2.ThetaDir = 'counterclockwise'; 
ax2.ThetaZeroLocation = 'right'; 
ax2.ThetaTick = 0:45:315;

% Radial scale
rmax = max(binVals);
if rmax == 0, rmax = 1; end
ax2.RLim = [0, rmax*1.05];

% Nice integer ticks for counts, or % ticks for proportions
if showProportion
    % use 0, 0.05, 0.10 ... up to rmax
    step = max(0.05, round(rmax/4, 2));
    rt = 0:step:rmax; 
    ax2.RTick = rt; 
    ax2.RTickLabel = strcat(string(round(rt*100)), "%");
else
    % choose ~4 ticks
    step = max(1, ceil(rmax/4));
    rt = 0:step:ceil(rmax);
    ax2.RTick = rt; 
    ax2.RTickLabel = string(rt);
end

ax2.FontSize = 12; ax2.GridColor = [.7 .7 .7]; ax2.LineWidth = 1;
title(sprintf('cs-on per %d° bin (%s)', 360/nBins, radialTitle), 'FontSize',12,'FontWeight','bold');

% Global annotation: N and bin width + what radius means
annotation(fig,'textbox',[0.35 0.92 0.3 0.06], ...
    'String', sprintf('N = %d cells | Bin width = %d° | Radius = %s', ...
                      numel(cs_on_ang_sac), 360/nBins, radialTitle), ...
    'FitBoxToText','on','EdgeColor','none','HorizontalAlignment','center');

ESN_Beautify_Plot(fig,[10 4]);
img_save_path = fullfile(data_path,'population_figs','ephys');
if ~exist(img_save_path,'dir'), mkdir(img_save_path); end
print(fig, fullfile(img_save_path,'cs-on_distribution.pdf'), '-dpdf', '-bestfit');


    % --------- choose 3 example cells (low, mid, high ρ) ----------
    rho_levels=[0.15 0.30];
    hi_thr=0.40; 
    used=false(size(cs_on_rho_sac));
    pick=@(t) find(~used & (abs(cs_on_rho_sac-t)==min(abs(cs_on_rho_sac(~used)-t))),1);
    idxs(1)=pick(rho_levels(1)); used(idxs(1))=true;
    idxs(2)=pick(rho_levels(2)); used(idxs(2))=true;
    cand=find(~used & cs_on_rho_sac>hi_thr,1); 
    if isempty(cand), [~,m]=max(cs_on_rho_sac(~used)); cand=find(~used,m,'first'); end
    idxs(3)=cand;

    % --------- robust event/time windows ----------
    nEvents = max(1,size(cs_rate,6));
    ev_vis  = min(1,nEvents);
    ev_sac  = min(2,nEvents);

    time_window_vis = [40,85];
    time_window_sac = [-70,30];

    T0=250;
    t_vis = (T0+time_window_vis(1)):(T0+time_window_vis(2));
    t_sac = (T0+time_window_sac(1)):(T0+time_window_sac(2));
    
    Tmax=size(cs_rate,2);
    t_vis=t_vis(t_vis>=1 & t_vis<=Tmax); if isempty(t_vis), t_vis=1:Tmax; end
    t_sac=t_sac(t_sac>=1 & t_sac<=Tmax); if isempty(t_sac), t_sac=1:Tmax; end

    % --------- per-cell figures (SEPARATE) ----------
    for k = 1:numel(idxs)
        c = idxs(k);

        % mean rate vs direction (avg across time/reward/vel)
        dir_rate_vis = mean_over_dims(squeeze(cs_rate(c, t_vis, :, :, :, ev_vis)), [1 3 4]); % -> [dir]
        dir_rate_sac = mean_over_dims(squeeze(cs_rate(c, t_sac, :, :, :, ev_sac)), [1 3 4]); % -> [dir]

        % angles (adjust if your order differs)
        thetas = deg2rad(0:45:315);
        [thetas_vis, rates_vis] = close_curve(thetas, dir_rate_vis);
        [thetas_sac, rates_sac] = close_curve(thetas, dir_rate_sac);

        % common radial limit for the two panels of this cell
        rmax = 1.1*max([rates_vis; rates_sac], [], 'omitnan'); if ~isfinite(rmax)||rmax<=0, rmax=1; end

        % make figure
        figC = figure('Position',[300 200 800 360],'Color','w');
        t = tiledlayout(figC,1,2,'TileSpacing','compact','Padding','compact');
        
        % ---- VISUAL panel (tile 1) ----
        axV = polaraxes(t);                 % parent is the layout
        axV.Layout.Tile = 1;                % put it in tile #1
        hold(axV,'on');
        style_polar(axV, rmax);
        axV.Color = 'none';
        fill_polar(axV, thetas_vis, rates_vis, [0.85 0.2 0.2], 0.25);
        %polarplot(axV, thetas_vis, rates_vis, 'k', 'LineWidth', 1.8);
        %polarscatter(axV, thetas, dir_rate_vis, 28, 'k', 'filled', 'MarkerFaceAlpha', 0.9);
        arrow_len_vis = cs_on_rho_vis(c) * rmax * 0.95;
        arrow_ang_vis = deg2rad(cs_on_ang_vis(c));
        polarplot(axV, [0 arrow_ang_vis], [0 arrow_len_vis], 'LineWidth', 2.2, 'Color', [0.95 0.7 0.1]);
        title(axV, sprintf('visual  |  \\rho=%.2f @%d°', cs_on_rho_vis(c), round(cs_on_ang_vis(c))), ...
              'FontSize',11,'FontWeight','bold');
        
        % ---- MOTOR panel (tile 2) ----
        axM = polaraxes(t);
        axM.Layout.Tile = 2;                % put it in tile #2
        hold(axM,'on');
        style_polar(axM, rmax);
        axM.Color = 'none';
        fill_polar(axM, thetas_sac, rates_sac, [0.2 0.35 0.85], 0.25);
        %polarplot(axM, thetas_sac, rates_sac, 'k', 'LineWidth', 1.8);
        %polarscatter(axM, thetas, dir_rate_sac, 28, 'k', 'filled', 'MarkerFaceAlpha', 0.9);
        arrow_len_sac = cs_on_rho_sac(c) * rmax * 0.95;
        arrow_ang_sac = deg2rad(cs_on_ang_sac(c));
        polarplot(axM, [0 arrow_ang_sac], [0 arrow_len_sac], 'LineWidth', 2.2, 'Color', [0.95 0.7 0.1]);
        title(axM, sprintf('motor  |  \\rho=%.2f @%d°', cs_on_rho_sac(c), round(cs_on_ang_sac(c))), ...
              'FontSize',11,'FontWeight','bold');

        % ---- Neuron ID under the figure ----
        neuron_id = get_cell_id(cs_cell_ids{c});
        subject   = infer_subject_from_id(cs_cell_ids{c}); 

        % ---- Save with ID in filename ----
        ESN_Beautify_Plot(figC,[10 4]);   % do this first
        
        id_str = sprintf('Subject: %s  |  Neuron ID: %s', subject, neuron_id);
        annotation(figC,'textbox', ...
            'Units','normalized', ...
            'Position',[0 0.005 1 0.04], ...
            'String', id_str, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'EdgeColor','none', ...
            'Interpreter','none', ...
            'FontSize',9);
        
        safe_id      = regexprep(neuron_id,'[^a-zA-Z0-9_.-]','_');
        safe_subject = regexprep(subject,  '[^a-zA-Z0-9_.-]','_');
        
        print(figC, fullfile(data_path,'population_figs','ephys', ...
            sprintf('%s_%s.pdf', safe_subject, safe_id)), '-dpdf','-bestfit');
 
    end
end

% ---- helpers ----
function X = mean_over_dims(X, dims)
    for d = sort(dims,'descend')
        if ndims(X) >= d && size(X,d) > 1
            X = mean(X, d, 'omitnan');
        end
    end
    X = squeeze(X);
end

function [th_closed, r_closed] = close_curve(thetas, r)
    r = r(:);
    if numel(r) ~= numel(thetas)
        error('Direction count mismatch. Expected %d, got %d.', numel(thetas), numel(r));
    end
    th_closed = [thetas(:); thetas(1)];
    r_closed  = [r; r(1)];
end

function fill_polar(axPolar, thetas_closed, r_closed, colorRGB, alphaVal)
    % Underlay a filled polygon by overlaying a cartesian axes
    rmax = axPolar.RLim(2);
    fig  = ancestor(axPolar,'figure');

    axCart = axes('Parent', fig, 'Units', axPolar.Units, ...
                  'Position', axPolar.Position, 'Color', 'none', ...
                  'Visible', 'off', 'HitTest', 'off');
    [x, y] = pol2cart(thetas_closed, r_closed);

    patch('XData', x, 'YData', y, 'FaceColor', colorRGB, ...
          'FaceAlpha', alphaVal, 'EdgeColor', 'none', 'Parent', axCart);

    axis(axCart, 'equal');
    xlim(axCart, [-rmax rmax]);
    ylim(axCart, [-rmax rmax]);

    uistack(axCart, 'bottom'); % keep under polar
end

function style_polar(ax, rmax)
    ax.ThetaZeroLocation='right';
    ax.ThetaDir='counterclockwise';
    ax.ThetaTick=0:45:315;
    ax.RLim=[0 rmax];
    ax.GridColor=[.7 .7 .7];
    ax.LineWidth=1.1;
    ax.FontSize=11;
end

function out = get_cell_id(fname)
      [~,b,~] = fileparts(fname);
    m = regexp(b,'^(.*?)(?=_?combine\b)','tokens','once','ignorecase');
    if isempty(m), m = {b}; end
    out = regexprep(m{1},'[^0-9_]','');
end
function subj = infer_subject_from_id(cell_id_str)
    dt = extract_date_from_string(cell_id_str);
    if ~isnat(dt) && dt >= datetime(2024,3,11) && dt <= datetime(2024,9,16)
        subj = '132F';
    elseif ~isnat(dt)
        subj = '65F';
    else
        subj = 'Unknown';
    end
end

function dt = extract_date_from_string(s)
    % Return datetime from strings like:
    % YYYY[-_.]?MM[-_.]?DD, YYYYMMDD, or YYMMDD (e.g., '240712...')
    dt = NaT;

    % 1) YYYY[-_.]?MM[-_.]?DD
    tok = regexp(s,'(?<!\d)(20\d{2})[-_\.]?(\d{1,2})[-_\.]?(\d{1,2})(?!\d)','tokens','once');
    if ~isempty(tok)
        y = str2double(tok{1}); m = str2double(tok{2}); d = str2double(tok{3});
        if isvalid_ymd(y,m,d), dt = datetime(y,m,d); return; end
    end

    % 2) YYYYMMDD
    tok = regexp(s,'(?<!\d)(20\d{2})(\d{2})(\d{2})(?!\d)','tokens','once');
    if ~isempty(tok)
        y = str2double(tok{1}); m = str2double(tok{2}); d = str2double(tok{3});
        if isvalid_ymd(y,m,d), dt = datetime(y,m,d); return; end
    end

    % 3) YYMMDD (prefer at start, then anywhere)
    tok = regexp(s,'^(?<!\d)(\d{2})(\d{2})(\d{2})(?!\d)','tokens','once');
    if isempty(tok)
        tok = regexp(s,'(?<!\d)(\d{2})(\d{2})(\d{2})(?!\d)','tokens','once');
    end
    if ~isempty(tok)
        yy = str2double(tok{1}); mm = str2double(tok{2}); dd = str2double(tok{3});
        % pivot: map 00–69 -> 2000–2069, 70–99 -> 1970–1999 (adjust if you prefer)
        if mm>=1 && mm<=12 && dd>=1 && dd<=31
            if yy <= 69, y = 2000 + yy; else, y = 1900 + yy; end
            if isvalid_ymd(y,mm,dd), dt = datetime(y,mm,dd); end
        end
    end
end

function tf = isvalid_ymd(y,m,d)
    tf = ~isnan(y) && ~isnan(m) && ~isnan(d) && m>=1 && m<=12 && d>=1 && d<=31;
    if tf
        try, datetime(y,m,d); catch, tf = false; end
    end
end

