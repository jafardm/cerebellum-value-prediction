function plot_projected_space_prim_sac_vel(data_path)
% Plot SS rate by explicit velocity bins saved in dataset meta.vel_edges.
% Makes two PDFs:
%   - ss-rate-Projection_prim-sac_byVel_pY.pdf  (proj onto sac+90)
%   - ss-rate-Projection_prim-sac_byVel_pX.pdf  (proj onto sac)
%
% Works with any number of velocity bins (uses size of 5th dim or meta.vel_edges).

img_save_path = fullfile(data_path,'population_figs','ephys');
save_path     = fullfile(data_path,'population_data');
if ~exist(img_save_path,'dir'); mkdir(img_save_path); end

JDM_params_funcs

S = load(fullfile(save_path,'SS_population_clique_prim_sac.mat'),'data','meta');
data_ss = S.data;
meta    = [];
if isfield(S,'meta'); meta = S.meta; end

% ---- essentials & weights ----
load(fullfile(save_path,'purkinje_cell_ids.mat'),'ind_p','ind_b');
sess_SS    = cell2mat(cellfun(@(x) str2double(x(1:6)), data_ss.cell_ids_tot,'UniformOutput',false));
cliques_SS = sess_SS.*10 + data_ss.cliques;
num_SS     = numel(cliques_SS);
num_b      = sum(ind_b);
num_p      = sum(ind_p);
ss_cs_rho  = data_ss.cs_on_rho_tot;

% ---- arrays (velocity-binned) ----
if ~isfield(data_ss,'rate_tot_sac_vel')
    error('SS file missing rate_tot_sac_vel. Rebuild the dataset with velocity bins.');
end
rate_vel = data_ss.rate_tot_sac_vel;                 % [cell x time x dir x cond x nVel]

have_vm_vel = isfield(data_ss,'vm_tot_sac_vel');
if have_vm_vel
    vm_vel = data_ss.vm_tot_sac_vel;                 % [cell x time x dir x cond x nVel]
else
    vm_vel = [];
    warning('vm_tot_sac_vel not found. Will fallback to vm_tot_sac (no per-vel shading).');
end

% ---- how many velocity bins? ----
n_vel_bins = size(rate_vel,5);
VE = [];
if isfield(meta,'vel_edges') && ~isempty(meta.vel_edges)
    VE = meta.vel_edges(:)';                          % edges length = n_vel_bins+1
    % Safety: if edges and data disagree, prefer data count and clip edges
    if numel(VE) ~= (n_vel_bins+1)
        VE = [];  % avoid misleading labels
    end
end

% ---- labels for velocity bins ----
if ~isempty(VE)
    vel_bin_labels = cell(1,n_vel_bins);
    for v = 1:n_vel_bins
        lo = VE(v); hi = VE(v+1);
        if isfinite(hi)
            vel_bin_labels{v} = sprintf('%.0f–%.0f°/s', lo, hi);
        else
            vel_bin_labels{v} = sprintf('≥%.0f°/s', lo);
        end
    end
else
    % generic labels if edges unavailable
    vel_bin_labels = arrayfun(@(k) sprintf('Vel bin %d',k), 1:n_vel_bins, 'uni', 0);
end

% ---- grouping & plotting params ----
cond_labels = {'High','Low'};
high_idx = [1 2];
low_idx  = [3 4];

vm_color = [204,174,98]/256;
vm_lim   = [0,600];
time_ind = (-100:150)+250;           % -100..+150 around decel
proj_lim = [-35,55];

% reducer for VM: [cell x time x dir] -> [1 x numel(time_ind)]
vm_reduce = @(A, idx) squeeze(mean(mean(A(:,idx,:),3,'omitnan'),1));

%% =========================
%  Figure A: p_y(t) (proj onto sac+90) for all velocity bins
%  =========================
fig = figure;
t = tiledlayout(2, n_vel_bins, 'TileSpacing','compact','Padding','compact');
title(t, 'Purkinje Projection — SS by Velocity Bin (p_y(t) = cs-on+90)')

for i_group = 1:2
    if i_group == 1
        cond_idx = high_idx; cond_name = cond_labels{1};
    else
        cond_idx = low_idx;  cond_name = cond_labels{2};
    end

    for v = 1:n_vel_bins
        % Average selected conditions first
        current_data = mean(rate_vel(:,:,:,cond_idx,v), 4, 'omitnan');   % [cell x time x dir]
        if ~isempty(vm_vel)
            current_vm = mean(vm_vel(:,:,:,cond_idx,v), 4, 'omitnan');   % [cell x time x dir]
        elseif isfield(data_ss,'vm_tot_sac')
            current_vm = mean(data_ss.vm_tot_sac(:,:,:,cond_idx), 4, 'omitnan'); % coarse fallback
        else
            current_vm = [];
        end

        ax = nexttile((i_group-1)*n_vel_bins + v); hold(ax,'on'); colororder({'k','k'})
        proj_w = -cosd((0:45:360-45)+90);  % p_y

        % --- Bursters ---
        rate_pr_b = zeros(num_b,numel(time_ind));
        for d = 1:8
            rate_pr_b = rate_pr_b + current_data(ind_b, time_ind, d) .* proj_w(d);
        end
        rate_pr_b = rate_pr_b * num_b .* ss_cs_rho(ind_b) / sum(ss_cs_rho(ind_b));

        yyaxis left
        [hl,hp] = boundedline(time_ind-250, ...
            mean(rate_pr_b,'omitnan'), ...
            std(rate_pr_b,0,'omitnan')/sqrt(num_b), 'r','alpha'); %#ok<*NBRAK>
        hl.DisplayName = sprintf('SS burst (n=%d)', num_b); hp.HandleVisibility='off';
        ylabel('Firing rate (Hz)')

        % --- Pausers ---
        rate_pr_p = zeros(num_p,numel(time_ind));
        for d = 1:8
            rate_pr_p = rate_pr_p + current_data(ind_p, time_ind, d) .* proj_w(d);
        end
        rate_pr_p = rate_pr_p * num_p .* ss_cs_rho(ind_p) / sum(ss_cs_rho(ind_p));

        [hl,hp] = boundedline(time_ind-250, ...
            mean(rate_pr_p,'omitnan'), ...
            std(rate_pr_p,0,'omitnan')/sqrt(num_p), 'b','alpha');
        hl.DisplayName = sprintf('SS pause (n=%d)', num_p); hp.HandleVisibility='off';

        % --- All SS ---
        rate_pr_ss = zeros(num_SS,numel(time_ind));
        for d = 1:8
            rate_pr_ss = rate_pr_ss + current_data(:, time_ind, d) .* proj_w(d);
        end
        rate_pr_ss = rate_pr_ss * num_SS .* ss_cs_rho / sum(ss_cs_rho);

        [hl,hp] = boundedline(time_ind-250, ...
            mean(rate_pr_ss,'omitnan'), ...
            std(rate_pr_ss,0,'omitnan')/sqrt(num_SS), 'k','alpha');
        hl.DisplayName = sprintf('SS all (n=%d)', num_SS); hp.HandleVisibility='off';

        xline(0,'--','HandleVisibility','off'); yline(0,'--','HandleVisibility','off')
        xlim([-100,100]); ylim(proj_lim)
        xlabel('Time from deceleration')

        yyaxis right
        if ~isempty(current_vm)
            vel_trace = vm_reduce(current_vm, time_ind);
            if any(isfinite(vel_trace))
                area(time_ind-250, vel_trace, 'FaceColor', vm_color, ...
                    'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                ylim(vm_lim); ylabel('Velocity (°/s)')
            end
        end

        if i_group == 1 && v == 1
            legend('Box','off','Location','northwest')
        end

        if ~isempty(vel_bin_labels)
            title(sprintf('%s — %s', cond_name, vel_bin_labels{v}))
        else
            title(sprintf('%s — Vel bin %d', cond_name, v))
        end
    end
end

% scale figure width with n_vel_bins to avoid clipped titles
ESN_Beautify_Plot(fig, [3.8*n_vel_bins, 7.5]);
print(fig, fullfile(img_save_path,'ss-rate-Projection_prim-sac_byVel_pY.pdf'), '-dpdf','-bestfit');

%% =========================
%  Figure B: p_x(t) (proj onto sac) for all velocity bins
%  =========================
fig = figure;
t = tiledlayout(2, n_vel_bins, 'TileSpacing','compact','Padding','compact');
title(t, 'Purkinje Projection — SS by Velocity Bin (p_x(t) = cs-on)')

for i_group = 1:2
    if i_group == 1
        cond_idx = high_idx; cond_name = cond_labels{1};
    else
        cond_idx = low_idx;  cond_name = cond_labels{2};
    end

    for v = 1:n_vel_bins
        current_data = mean(rate_vel(:,:,:,cond_idx,v), 4, 'omitnan');   % [cell x time x dir]
        if ~isempty(vm_vel)
            current_vm = mean(vm_vel(:,:,:,cond_idx,v), 4, 'omitnan');
        elseif isfield(data_ss,'vm_tot_sac')
            current_vm = mean(data_ss.vm_tot_sac(:,:,:,cond_idx), 4, 'omitnan');
        else
            current_vm = [];
        end

        ax = nexttile((i_group-1)*n_vel_bins + v); hold(ax,'on'); colororder({'k','k'})
        proj_w = -cosd((0:45:360-45));            % p_x

        % --- Bursters ---
        rate_pr_b = zeros(num_b,numel(time_ind));
        for d = 1:8
            rate_pr_b = rate_pr_b + current_data(ind_b, time_ind, d) .* proj_w(d);
        end
        rate_pr_b = rate_pr_b * num_b .* ss_cs_rho(ind_b) / sum(ss_cs_rho(ind_b));

        yyaxis left
        [hl,hp] = boundedline(time_ind-250, ...
            mean(rate_pr_b,'omitnan'), ...
            std(rate_pr_b,0,'omitnan')/sqrt(num_b), 'r','alpha');
        hl.DisplayName = 'SS burst'; hp.HandleVisibility='off';

        % --- Pausers ---
        rate_pr_p = zeros(num_p,numel(time_ind));
        for d = 1:8
            rate_pr_p = rate_pr_p + current_data(ind_p, time_ind, d) .* proj_w(d);
        end
        rate_pr_p = rate_pr_p * num_p .* ss_cs_rho(ind_p) / sum(ss_cs_rho(ind_p));

        [hl,hp] = boundedline(time_ind-250, ...
            mean(rate_pr_p,'omitnan'), ...
            std(rate_pr_p,0,'omitnan')/sqrt(num_p), 'b','alpha');
        hl.DisplayName = 'SS pause'; hp.HandleVisibility='off';

        % --- All SS ---
        rate_pr_ss = zeros(num_SS,numel(time_ind));
        for d = 1:8
            rate_pr_ss = rate_pr_ss + current_data(:, time_ind, d) .* proj_w(d);
        end
        rate_pr_ss = rate_pr_ss * num_SS .* ss_cs_rho / sum(ss_cs_rho);

        [hl,hp] = boundedline(time_ind-250, ...
            mean(rate_pr_ss,'omitnan'), ...
            std(rate_pr_ss,0,'omitnan')/sqrt(num_SS), 'k','alpha');
        hl.DisplayName = 'SS all'; hp.HandleVisibility='off';

        xline(0,'--','HandleVisibility','off'); yline(0,'--','HandleVisibility','off')
        xlim([-100,100]); ylim(proj_lim)
        xlabel('Time from deceleration')

        yyaxis right
        if ~isempty(current_vm)
            vel_trace = vm_reduce(current_vm, time_ind);
            if any(isfinite(vel_trace))
                area(time_ind-250, vel_trace, 'FaceColor', vm_color, ...
                    'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                ylim(vm_lim); ylabel('Velocity (°/s)')
            end
        end

        if i_group == 1 && v == 1
            legend('Box','off','Location','northwest')
        end

        if ~isempty(vel_bin_labels)
            title(sprintf('%s — %s', cond_name, vel_bin_labels{v}))
        else
            title(sprintf('%s — Vel bin %d', cond_name, v))
        end
    end
end

ESN_Beautify_Plot(fig, [3.8*n_vel_bins, 7.5]);
print(fig, fullfile(img_save_path,'ss-rate-Projection_prim-sac_byVel_pX.pdf'), '-dpdf','-bestfit');

end
