function res = JDM_load_sess_adj_matrix(path_data_monkey_sorted,sess,selected_types,keep_simul_cell,is_sort,is_jitt)

JDM_params_funcs;

if nargin  == 2
    selected_types = {'PC','SS','MLI'};
    keep_simul_cell = 1;
    is_sort = 1;
    is_jitt = 1;
end

units_info = MAF_extract_cell_metadata(path_data_monkey_sorted,sess);

% jitter for the position of the cells
if ~is_jitt
    jitt_x = 0;
else
    if numel(unique(units_info.elec.x)) == 1
        jitt_x = 10*randn(size(units_info.cell_x));
    else
        jitt_x = 5*randn(size(units_info.cell_x));
    end
end


load([units_info.units_path...
    units_info.sess_name '_interactions.mat'],...
    'cross_prob','cross_prob_jitt');

num_cells = size(cross_prob,1);

Adj = nan(num_cells);

xprob_corrected = nan(num_cells,num_cells,101);
xprob_corrected_cs = nan(num_cells,num_cells,101);
xprob_corrected_cs_ss = nan(num_cells,num_cells,101);

for counter_cell1 = 1:num_cells
    parfor counter_cell2 = 1:num_cells
        current_cross_prob      = cross_prob{counter_cell1,counter_cell2};
        current_cross_prob_jitt = cross_prob_jitt{counter_cell1,counter_cell2};

        if isempty(current_cross_prob)
            continue;
        end

        span = current_cross_prob.span;
        ind_span = span<=2 & span>=0;

        xprob_corrected_ = ...
            current_cross_prob.ss1xss2 -...
            current_cross_prob_jitt.ss1xss2;
        ss1xss2_corrected = xprob_corrected_(ind_span);
        [~,ind_inter] = max(abs(ss1xss2_corrected));

        xprob_corrected(counter_cell1,counter_cell2,:) = xprob_corrected_;
        xprob_corrected_cs(counter_cell1,counter_cell2,:) = ...
            current_cross_prob.cs1xcs2 -...
            current_cross_prob_jitt.cs1xcs2;

        xprob_corrected_cs_ss(counter_cell1,counter_cell2,:) = ...
            current_cross_prob.cs1xss2 -...
            current_cross_prob_jitt.cs1xss2;

        Adj(counter_cell1,counter_cell2) = ...
            ss1xss2_corrected(ind_inter);
    end
end
Adj = max(abs(Adj),abs(Adj)').*sign((Adj+Adj')/2);

cell_types = units_info.cell_type;

selected_cells = find(ismember(cell_types,selected_types));

if keep_simul_cell
    % only keep simultaneous cells
    criteria = 1;
    while(criteria == 1)
        Adj_ = Adj(selected_cells,selected_cells);
        num_nan_units = sum(isnan(Adj_));
        criteria = any(num_nan_units > 1);
        if criteria == 0
            break;
        end
        [~,ind_rm] = max(sum(isnan(Adj_)));
        selected_cells(ind_rm) = [];
    end
end


Adj = Adj(selected_cells,selected_cells);

cell_types = cellfun(@(x) x,units_info.cell_type,'UniformOutput',false);
cell_names = cellfun(@(x,y) [x, ' ' regexprep(y(15:end-14),'_','-')],cell_types,...
    units_info.cell_list,'UniformOutput',false);

cell_ids = units_info.cell_list(selected_cells);

cell_x = (units_info.cell_x(selected_cells) + jitt_x(selected_cells) + 30)*10;
cell_y = units_info.cell_y(selected_cells);

cell_names = cell_names(selected_cells);
cell_types = cell_types(selected_cells);

xprob_corrected = xprob_corrected(selected_cells,selected_cells,:);
xprob_corrected_cs = xprob_corrected_cs(selected_cells,selected_cells,:);
xprob_corrected_cs_ss = xprob_corrected_cs_ss(selected_cells,selected_cells,:);

Adj_ = abs(Adj);
Adj_(isnan(Adj_)) = 0;
if (length(Adj_) > 3) && (is_sort == 1)
    D = diag(sum(Adj_));
    % Dinvs = D^(-1/2);
    % L = eye(size(Adj_)) - Dinvs*Adj_*Dinvs;
    L = D - Adj_;
    [u,D] = eig(L);
    [~,ind_sorted] = sort(diag(D));
    dim_Sp = u(:,ind_sorted(2:3));

    cell_cliques= units_info.cell_clique(selected_cells);
    [AdjMat_xl,node_order] = sort_adj_matrix(Adj,cell_cliques,is_sort);
    dim_Sp              = dim_Sp(node_order,:);
else
    dim_Sp = 0;
    cell_cliques= units_info.cell_clique(selected_cells);
    AdjMat_xl = 0;
    node_order = 1:length(Adj_);
end

Adj                 = Adj(node_order,node_order);
cell_names          = cell_names(node_order);
cell_cliques        = cell_cliques(node_order);
cell_types          = cell_types(node_order);
xprob_corrected     = xprob_corrected(node_order,node_order,:);
xprob_corrected_cs  = xprob_corrected_cs(node_order,node_order,:);
xprob_corrected_cs_ss  = xprob_corrected_cs_ss(node_order,node_order,:);

res.Adj                = Adj;
res.Adj_xline          = AdjMat_xl;
res.cell_names         = cell_names;
res.cell_types         = cell_types;
res.dim_Sp             = dim_Sp;
res.cell_cliques       = cell_cliques;
res.xprob_corrected    = xprob_corrected;
res.xprob_corrected_cs = xprob_corrected_cs;
res.xprob_corrected_cs_ss = xprob_corrected_cs_ss;
res.cell_x             = cell_x(node_order);
res.cell_y             = cell_y(node_order);
res.cell_id            = cellfun(@(x) x(end-5:end),cell_names,'UniformOutput',false);
res.full_cell_id       = cell_ids(node_order);
res.selected_cells     = selected_cells;
res.node_order         = node_order;

end

function [AdjMat_xl,node_order] = sort_adj_matrix(current_Adj,current_labels,is_sort)

clusters = unique(current_labels);
num_clust = numel(clusters);

order = zeros(size(current_labels));
xl_ = zeros(num_clust+1,1);
for counter_clust = 1:num_clust
    current_clust = clusters(counter_clust);

    ind_ = find(current_labels == current_clust)';
    Adj = current_Adj(ind_,ind_);
    Adj_ = Adj;
    Adj_(isnan(Adj)) = 0;

    if numel(ind_)>2 && is_sort
        fig = figure;
        Z = linkage(Adj_,'single','cosine');
        [~,~,order_] = dendrogram(Z,200);
        close(fig);
    else
        order_ = 1:numel(ind_);
    end
    xl_(counter_clust+1) = numel(ind_) + xl_(counter_clust);
    order(xl_(counter_clust)+(1:numel(ind_))) = ind_(order_);
end
node_order = order;
xl_ = xl_ + .5;
AdjMat_xl = xl_;
end
