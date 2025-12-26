function plot_saccade_amp_vel_countour(data_path)
% Amplitude–Velocity only, using 2D KDE (ksdensity) + SURF.
% Each subject shows two tiles: Task-relevant (red) and Task-irrelevant (blue).
% Colorbar shows probability mass per grid cell (P sums ~1 over the plotted window).

    JDM_params_funcs;
    S = load(fullfile(data_path, 'population_data', 'behave_data.mat'));
    if isfield(S,'data_behav'), data_behav = S.data_behav; else, error('behave_data missing'); end

    subject_names = params.animal_list;
    amp_edges = params.pop.sac.amp_edges;  % linear units
    vel_edges = params.pop.sac.vel_edges;  % linear units

    % Groups
    tag_groups = {[1,4,6,7],10};
    group_titles = {'Targeted saccades (Tags 1,4,6,7)', 'Task-irrelevant (Tags 10)'};
    group_colors = {[1 0 0], [0 0 1]};      % red, blue

    % Log10 axis limits & ticks
    xL  = [-1, 1.5];                 % log10 amp (deg)
    yL  = [1, 3];                    % log10 vel (deg/s)
    xts = xL(1):xL(2);
    yts = yL(1):yL(2);

    % KDE grid (in log space)
    nx = 160; ny = 160;
    gx = linspace(xL(1), xL(2), nx);
    gy = linspace(yL(1), yL(2), ny);
    [GX, GY] = meshgrid(gx, gy);
    XI = [GX(:), GY(:)];
    dx = mean(diff(gx)); dy = mean(diff(gy));

    maxN = 2e5;  % optional cap for speed

    % Figure & layout: rows = subjects, columns = 2 groups
    fig = figure('Name','Saccade Amplitude–Velocity (KDE SURF w/ probability)');
    nSubjects = numel(subject_names);
    tl = tiledlayout(fig, nSubjects, 2, 'TileSpacing','compact', 'Padding','compact');

    for s = 1:nSubjects
        amp = data_behav.sac_amp{s};
        vel = data_behav.sac_vel{s};
        tag = data_behav.sac_tag{s};

        valid = isfinite(amp) & isfinite(vel) & amp>0 & vel>0;
        if ~any(valid), continue; end

        % Global GM regression (for overlay)
        b_ms = gmregress(log10(amp(valid)), log10(vel(valid)));
        xx = linspace(xL(1), xL(2), 200);
        yy = b_ms(1) + b_ms(2)*xx;

        % Precompute both group probability grids to unify color scaling
        P_all = cell(1,2);
        for g = 1:2
            mask = valid & ismember(tag, tag_groups{g});
            if ~any(mask)
                P_all{g} = zeros(ny,nx);
                continue;
            end
            idx = find(mask);
            if numel(idx) > maxN
                idx = idx(randperm(numel(idx), maxN));
            end
            lx = log10(amp(idx));
            ly = log10(vel(idx));

            % 2D KDE pdf over (log amp, log vel), then mass per cell
            f = ksdensity([lx ly], XI, 'Function','pdf');
            F = reshape(f, size(GX));
            P = F * dx * dy;
            % Re-normalize to sum to 1 over the plotted window (for a clean colorbar)
            P = P / max(1e-12, sum(P(:)));
            P_all{g} = P;
        end
       

        % --- user-tunable contrast knobs ---
        gamma_colormap     = 0.55;   % <1 boosts midtones; try 0.4–0.7
        clim_percentile    = 99.5;   % clip top % across both groups for better contrast

       
        % Shared CLim using high-percentile (not max) for better contrast
        hi1 = prctile(P_all{1}(:), clim_percentile);
        hi2 = prctile(P_all{2}(:), clim_percentile);
        clim_max = max([hi1, hi2, 1e-9]);   % avoid zero upper bound
        for g = 1:2
            ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
            title(ax, sprintf('%s — %s', subject_names{s}, group_titles{g}));
        
            % SURF with boosted colormap
            hs = surf(ax, GX, GY, zeros(size(P_all{g})), P_all{g}, ...
                'EdgeColor','none', 'FaceAlpha', 1);
            view(ax, 2); shading(ax, 'interp');
        
            colormap(ax, cmap_to_color(group_colors{g}, 256, gamma_colormap));  % << gamma
            caxis(ax, [0, clim_max]);                                           % << percentile clip
            cb = colorbar(ax, 'Location','eastoutside');
            cb.Label.String = 'Probability mass';

            % Regression line & bin guides
            plot(ax, xx, yy, 'k', 'LineWidth', 1.1, 'HandleVisibility','off');
            if numel(amp_edges) > 2, xline(ax, log10(amp_edges(2:end-1)), '--k', 'HandleVisibility','off'); end
            if numel(vel_edges) > 2, yline(ax, log10(vel_edges(2:end-1)), '--k', 'HandleVisibility','off'); end

            P = P_all{g};
            hpd = 0.50;                               % half the mass
            thr = hpd_threshold_from_grid(P, hpd);    % mass-per-cell threshold
            contour(ax, GX, GY, P, [thr thr], 'k-', 'LineWidth', 0.8, ...
                    'HandleVisibility','off');

            % Axis formatting: show linear units on ticks
            xlim(ax, xL); ylim(ax, yL);
            xticks(ax, xts); yticks(ax, yts);
            xticklabels(ax, arrayfun(@(v) sprintf('%g',10.^v), xts, 'UniformOutput', false));
            yticklabels(ax, arrayfun(@(v) sprintf('%g',10.^v), yts, 'UniformOutput', false));
            xlabel(ax, 'amp (deg)'); ylabel(ax, 'vel (deg/s)');
        end
    end

    title(tl, 'Saccade Amplitude–Velocity — KDE Probability Mass');
    ESN_Beautify_Plot(fig, [12, 4 + 2*nSubjects]);  
    print(fig, fullfile(data_path, 'population_figs\behavior','saccade_amp_vel_countour.pdf'), '-dpdf', '-bestfit')
end

% --------- helpers ---------
function cmap = cmap_to_color(rgb, n, gamma)
% Monochrome colormap from white -> rgb with optional gamma boost
    if nargin < 2, n = 256; end
    if nargin < 3 || isempty(gamma), gamma = 1; end
    t = linspace(0,1,n)'.^gamma;   % gamma<1 boosts midrange
    white = ones(n,3);
    cmap  = white.*(1 - t) + t.*(ones(n,1)*rgb(:)');  % linear blend with gamma
end

function thr = hpd_threshold_from_grid(P, mass)
    p = P(:); [ps,~] = sort(p,'descend'); cs = cumsum(ps);
    j = find(cs >= mass - eps, 1, 'first');
    thr = ps(j);
end