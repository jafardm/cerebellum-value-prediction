function plot_cs_on_rate_velocityBins_prim_sac(data_path)
    % Plots CS rate for visual, saccade, and offset,
    % split by High/Low reward and velocity bins, aligned to CS-on.
    %
    % Assumes:
    %   cs_rate: [cell x time x dir x reward(2) x vel x event(3)]
    %   cs_vm  : same size (velocity magnitude)
    % events: 1=visual, 2=saccade, 3=offset
    % reward: 1=High, 2=Low  (edit rew_names if reversed)

    JDM_params_funcs;
    S = load(fullfile(data_path,'population_data','cs_on_rate_velocityBins_prim_sac.mat'));  
    cs_rate    = S.cs_rate;
    cs_vm      = S.cs_vm;
    cs_on_data = S.cs_on_data;

    % --- Extract rho/ang to choose CS-on vector per cell
    cs_on_rho_vis = cell2mat(cellfun(@(x) x.vis.rho_avg, cs_on_data,'UniformOutput',false));
    cs_on_ang_vis = cell2mat(cellfun(@(x) x.vis.ang_avg, cs_on_data,'UniformOutput',false));
    cs_on_rho_sac = cell2mat(cellfun(@(x) x.sac.rho_avg, cs_on_data,'UniformOutput',false));
    cs_on_ang_sac = cell2mat(cellfun(@(x) x.sac.ang_avg, cs_on_data,'UniformOutput',false));

    num_cs = numel(cs_on_data);

    % --- Choose larger vector (visual vs saccade) for CS-on
    ind_m = cs_on_rho_sac > cs_on_rho_vis;
    cs_on_ang = cs_on_ang_vis/180*pi;
    cs_on_ang(ind_m) = cs_on_ang_sac(ind_m)/180*pi;

    % --- Determine CS-on bin for each cell
    ang_bins = params.sac.ang_values/180*pi;        
    [~, cs_on_bin] = min(abs(angdiff(repmat(ang_bins,num_cs,1), cs_on_ang*ones(size(ang_bins)))), [], 2);

    order_ang = nan(num_cs,8);
    for i = 1:num_cs
        order_ang(i,:) = circshift(1:8, 1 - cs_on_bin(i)); % put CS-on at index 1
    end

    % ===== Rotate to CS-on–aligned space, preserving reward/vel/event =====
    sz = size(cs_rate);
    assert(numel(sz) == 6, 'cs_rate must be 6-D [cell,time,dir,reward,vel,event]');
    n_time   = sz(2);
    n_dir    = sz(3); %#ok<NASGU>
    n_rew    = sz(4);
    n_vel    = sz(5);
    n_event  = sz(6); %#ok<NASGU>

    cs_rate_rot = nan(num_cs, n_time, 8, n_rew, n_vel, 3);
    cs_vm_rot   = nan(num_cs, n_time, 8, n_rew, n_vel, 3);

    for i = 1:num_cs
        cs_rate_rot(i,:,:,:,:,:) = cs_rate(i,:, order_ang(i,1:8),:,:,:);
        cs_vm_rot(  i,:,:,:,:,:) = cs_vm(  i,:, order_ang(i,1:8),:,:,:);
    end

    % ===== Plotting settings =====
    labels_vel = {'low','medium','high'};
    n_bins     = numel(labels_vel);
    cmap_bins  = flipud(cool(n_bins));    % velocity colors slow→fast
    leg_dirs   = {'cs-on','cs+180'};
    dir_map    = {1, 5};
    rew_names  = {'High','Low'};
    rew_styles = {'-','--'};

    % Time indices
    t_ind_vis = (-50:200)  + 250;
    t_ind_sac = (-100:150) + 250;
    t_ind_off = (-100:200) + 250;
    T = n_time;
    t_ind_vis = t_ind_vis(t_ind_vis>=1 & t_ind_vis<=T);
    t_ind_sac = t_ind_sac(t_ind_sac>=1 & t_ind_sac<=T);
    t_ind_off = t_ind_off(t_ind_off>=1 & t_ind_off<=T);

    vm_color = [204, 174, 98]/256;
    vm_alpha = 0.25;

    % === Figures ===
    make_event_figure(1, 'Visual',  t_ind_vis, cs_rate_rot, cs_vm_rot, dir_map, leg_dirs, ...
                      rew_names, rew_styles, labels_vel, cmap_bins, vm_color, vm_alpha, num_cs);
    make_event_figure(2, 'Saccade', t_ind_sac, cs_rate_rot, cs_vm_rot, dir_map, leg_dirs, ...
                      rew_names, rew_styles, labels_vel, cmap_bins, vm_color, vm_alpha, num_cs);
%     make_event_figure(3, 'Offset',  t_ind_off, cs_rate_rot, cs_vm_rot, dir_map, leg_dirs, ...
%                       rew_names, rew_styles, labels_vel, cmap_bins, vm_color, vm_alpha, num_cs);
end

% ---------- Local helpers ----------
function make_event_figure(event_idx, event_name, t_ind, cs_rate_rot, cs_vm_rot, ...
                           dir_map, leg_dirs, rew_names, rew_styles, labels_vel, cmap_bins, ...
                           vm_color, vm_alpha, num_cs)
    n_rew  = numel(rew_names);
    n_bins = numel(labels_vel);
    n_tiles = numel(leg_dirs) + 1; % +1 for velocity summary panel

    fig = figure('Name', event_name);
    tiledlayout(fig, 1, n_tiles, 'TileSpacing','compact','Padding','compact');

    % ----- First tiles: CS rate panels -----
    for d = 1:numel(leg_dirs)
        dir_idx    = dir_map{d};
        block_rate = avg_over_dirs(cs_rate_rot(:, t_ind, :, :, :, event_idx), dir_idx);
        block_vm   = avg_over_dirs(cs_vm_rot(:,   t_ind, :, :, :, event_idx), dir_idx);

        nexttile; hold on; box on; colororder({'k','k'});

        % VM overlay
        vm_mean = squeeze(mean(block_vm, [1 4], 'omitnan')); % [time x reward]
        vm_mean_collapsed = mean(vm_mean, 2, 'omitnan');

        yyaxis right
        ylim([0 450])
        area((t_ind - 250), vm_mean_collapsed, ...
             'FaceColor', vm_color, 'FaceAlpha', vm_alpha, ...
             'EdgeColor','none','HandleVisibility','off');
        ylabel('Velocity (deg/s)');

        % CS rate traces
        yyaxis left
        for r = 1:n_rew
            for vb = 1:n_bins
                X  = squeeze(block_rate(:,:,r,vb));           
                mu = mean(X, 1, 'omitnan');
                n_eff = sum(~any(isnan(X),2));
                se = std(X, 0, 1, 'omitnan') ./ max(1, sqrt(n_eff));

                [hl1, hl2] = boundedline((t_ind - 250), mu, se, 'alpha');
                set(hl1, 'Color',   cmap_bins(vb,:), ...
                         'LineStyle', rew_styles{r}, ...
                         'DisplayName', sprintf('%s | %s', rew_names{r}, labels_vel{vb}));
                set(hl2, 'HandleVisibility','off', 'FaceColor', cmap_bins(vb,:));
            end
        end

        title(leg_dirs{d});
        if d == 1, ylabel('CS rate (Hz)'); end
        xlabel(sprintf('Time from %s onset (ms)', lower(event_name)));
        xline(0,'--','HandleVisibility','off');
    end

    % ----- Last tile: velocity summary (low/med/high) -----
        all_dirs = [dir_map{:}];
    % vm_all_dirs: [cell x time x reward x vel]
    vm_all_dirs = avg_over_dirs(cs_vm_rot(:, t_ind, :, :, :, event_idx), all_dirs);

    nexttile; hold on; box on;
    for r = 1:n_rew
        for vb = 1:n_bins
            % Xv: [cell x time]
            Xv  = squeeze(vm_all_dirs(:,:,r,vb));
            muv = mean(Xv, 1, 'omitnan');
            n_eff_v = sum(~any(isnan(Xv),2));
            sev = std(Xv, 0, 1, 'omitnan') ./ max(1, sqrt(n_eff_v));

            [hl1, hl2] = boundedline((t_ind - 250), muv, sev, 'alpha');
            set(hl1, 'Color', cmap_bins(vb,:), ...
                     'LineStyle', rew_styles{r}, ...
                     'DisplayName', sprintf('%s | %s', rew_names{r}, labels_vel{vb}));
            set(hl2, 'HandleVisibility','off', 'FaceColor', cmap_bins(vb,:));
        end
    end
    xlabel(sprintf('Time from %s onset (ms)', lower(event_name)));
    ylabel('Velocity (deg/s)');
    title('Velocity (mean \pm SEM) by reward');
    xline(0,'--','HandleVisibility','off');
    legend('Box','off','Location','northeast');

    sgtitle(sprintf('%s: reward × velocity bins   |   CS Cells (n=%d)', event_name, num_cs));
    ESN_Beautify_Plot(fig, [16, 6]); 

    % Save
    outdir = 'C:\Users\Jafar\Documents\reward\population_figs\ephys';
    print(fig, fullfile(outdir, [event_name '-velocity_bins.pdf']), '-dpdf', '-bestfit');
end

function Y = avg_over_dirs(X, dir_idx)
    if numel(dir_idx) == 1
        Y = squeeze(X(:,:,dir_idx,:,:));
    else
        Y = squeeze(mean(X(:,:,dir_idx,:,:), 3, 'omitnan'));
    end
end
