function JDM_build_clique_interaction_summary
JDM_params_funcs

user = 'JDM';
loc  = 'ctx';

path_data_monkey_sorted = params.path_data_monkey_sorted;
animal_list = params.animal_list;
session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);

num_animal = numel(session_list);

% 1:SS,2:ML,3:ML2,4:CS,5:SS|ML,6:SS|ML2,7:ML|ML2,8:ML|CS,9:ML2|CS,10:SS|CS
bet_inter = cell(10,1);
wit_inter = cell(10,1);

bet_inter_ids = cell(10,1);
wit_inter_ids = cell(10,1);

for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess = numel(current_sess_list);

    bet_inter_sess_ss = cell(num_sess,1);
    bet_inter_sess_ml = cell(num_sess,1);
    bet_inter_sess_ml2 = cell(num_sess,1);
    bet_inter_sess_cs = cell(num_sess,1);
    bet_inter_sess_ss_ml = cell(num_sess,1);
    bet_inter_sess_ss_ml2 = cell(num_sess,1);
    bet_inter_sess_ml_ml2 = cell(num_sess,1);
    bet_inter_sess_ml_cs = cell(num_sess,1);
    bet_inter_sess_ml2_cs = cell(num_sess,1);
    bet_inter_sess_ss_cs = cell(num_sess,1);

    wit_inter_sess_ss = cell(num_sess,1);
    wit_inter_sess_ml = cell(num_sess,1);
    wit_inter_sess_ml2 = cell(num_sess,1);
    wit_inter_sess_cs = cell(num_sess,1);
    wit_inter_sess_ss_ml = cell(num_sess,1);
    wit_inter_sess_ss_ml2 = cell(num_sess,1);
    wit_inter_sess_ml_ml2 = cell(num_sess,1);
    wit_inter_sess_ml_cs = cell(num_sess,1);
    wit_inter_sess_ml2_cs = cell(num_sess,1);
    wit_inter_sess_ss_cs = cell(num_sess,1);

    bet_inter_sess_ss_ids = cell(num_sess,1);
    bet_inter_sess_ml_ids = cell(num_sess,1);
    bet_inter_sess_ml2_ids = cell(num_sess,1);
    bet_inter_sess_cs_ids = cell(num_sess,1);
    bet_inter_sess_ss_ml_ids = cell(num_sess,1);
    bet_inter_sess_ss_ml2_ids = cell(num_sess,1);
    bet_inter_sess_ml_ml2_ids = cell(num_sess,1);
    bet_inter_sess_ml_cs_ids = cell(num_sess,1);
    bet_inter_sess_ml2_cs_ids = cell(num_sess,1);
    bet_inter_sess_ss_cs_ids = cell(num_sess,1);

    wit_inter_sess_ss_ids = cell(num_sess,1);
    wit_inter_sess_ml_ids = cell(num_sess,1);
    wit_inter_sess_ml2_ids = cell(num_sess,1);
    wit_inter_sess_cs_ids = cell(num_sess,1);
    wit_inter_sess_ss_ml_ids = cell(num_sess,1);
    wit_inter_sess_ss_ml2_ids = cell(num_sess,1);
    wit_inter_sess_ml_ml2_ids = cell(num_sess,1);
    wit_inter_sess_ml_cs_ids = cell(num_sess,1);
    wit_inter_sess_ml2_cs_ids = cell(num_sess,1);
    wit_inter_sess_ss_cs_ids = cell(num_sess,1);
   parfor counter_sess = 1:num_sess
        current_sess = current_sess_list{counter_sess};
        disp(['     ' current_sess]);
        res = JDM_load_sess_adj_matrix(current_path,current_sess,{'SS','PC','MLI','MLI2'},0,0,1);
        cell_types_   = res.cell_types;
        cell_cliques_ = res.cell_cliques;
        cell_ids_     = res.full_cell_id;

        cliques_ = unique(cell_cliques_);
        n_cell   = numel(cell_cliques_);

        [ids_x,ids_y] = ndgrid(1:n_cell,1:n_cell);
        ids_mat = cat(3,cell_ids_(ids_x),cell_ids_(ids_y));

        xprob_corrected_ = reshape(res.xprob_corrected,n_cell^2,101);
        xprob_corrected_cs_ = reshape(res.xprob_corrected_cs,n_cell^2,101);
        xprob_corrected_cs_ss_ = reshape(res.xprob_corrected_cs_ss,n_cell^2,101);
        ids_mat_ = reshape(ids_mat,n_cell^2,2);

        wit_ind_ = zeros(n_cell,n_cell);

        for counter_clique = 1:numel(cliques_)
            current_clique = cliques_(counter_clique);
            clique_ind = cell_cliques_ == current_clique;
            if current_clique == 0
                continue;
            end
            wit_ind_(clique_ind,clique_ind) = 1;
        end
        bet_ind_ = not(wit_ind_);
        wit_ind_ = wit_ind_ - diag(diag(wit_ind_));
        bet_ind_ = bet_ind_ - diag(diag(bet_ind_));
        bet_ind_(cell_cliques_ == 0,:) = 0;
        bet_ind_(:,cell_cliques_ == 0) = 0;

        ss_ind = ismember(cell_types_,{'SS','PC'});
        ml_ind = ismember(cell_types_,{'MLI'});
        cs_ind = ismember(cell_types_,{'PC'});
        ml2_ind = ismember(cell_types_,{'MLI2'});

        ind_ss_wit = wit_ind_;
        ind_ss_wit(~ss_ind,:) = 0;
        ind_ss_wit(:,~ss_ind) = 0;
        wit_inter_sess_ss{counter_sess} = xprob_corrected_(ind_ss_wit(:) == 1,:);
        wit_inter_sess_ss_ids{counter_sess} = ids_mat_(ind_ss_wit(:) == 1,:);

        ind_ss_bet = bet_ind_;
        ind_ss_bet(~ss_ind,:) = 0;
        ind_ss_bet(:,~ss_ind) = 0;
        bet_inter_sess_ss{counter_sess} = xprob_corrected_(ind_ss_bet(:) == 1,:);
        bet_inter_sess_ss_ids{counter_sess} = ids_mat_(ind_ss_bet(:) == 1,:);

        ind_ml_wit = wit_ind_;
        ind_ml_wit(~ml_ind,:) = 0;
        ind_ml_wit(:,~ml_ind) = 0;
        wit_inter_sess_ml{counter_sess} = xprob_corrected_(ind_ml_wit(:) == 1,:);
        wit_inter_sess_ml_ids{counter_sess} = ids_mat_(ind_ml_wit(:) == 1,:);

        ind_ml_bet = bet_ind_;
        ind_ml_bet(~ml_ind,:) = 0;
        ind_ml_bet(:,~ml_ind) = 0;
        bet_inter_sess_ml{counter_sess} = xprob_corrected_(ind_ml_bet(:) == 1,:);
        bet_inter_sess_ml_ids{counter_sess} = ids_mat_(ind_ml_bet(:) == 1,:);

        ind_ml2_wit = wit_ind_;
        ind_ml2_wit(~ml2_ind,:) = 0;
        ind_ml2_wit(:,~ml2_ind) = 0;
        wit_inter_sess_ml2{counter_sess} = xprob_corrected_(ind_ml2_wit(:) == 1,:);
        wit_inter_sess_ml2_ids{counter_sess} = ids_mat_(ind_ml2_wit(:) == 1,:);

        ind_ml2_bet = bet_ind_;
        ind_ml2_bet(~ml2_ind,:) = 0;
        ind_ml2_bet(:,~ml2_ind) = 0;
        bet_inter_sess_ml2{counter_sess} = xprob_corrected_(ind_ml2_bet(:) == 1,:);
        bet_inter_sess_ml2_ids{counter_sess} = ids_mat_(ind_ml2_bet(:) == 1,:);

        ind_cs_wit = wit_ind_;
        ind_cs_wit(~cs_ind,:) = 0;
        ind_cs_wit(:,~cs_ind) = 0;
        wit_inter_sess_cs{counter_sess} = xprob_corrected_cs_(ind_cs_wit(:) == 1,:);
        wit_inter_sess_cs_ids{counter_sess} = ids_mat_(ind_cs_wit(:) == 1,:);

        ind_cs_bet = bet_ind_;
        ind_cs_bet(~cs_ind,:) = 0;
        ind_cs_bet(:,~cs_ind) = 0;
        bet_inter_sess_cs{counter_sess} = xprob_corrected_cs_(ind_cs_bet(:) == 1,:);
        bet_inter_sess_cs_ids{counter_sess} = ids_mat_(ind_cs_bet(:) == 1,:);

        ind_ss_ml_wit = wit_ind_;
        ind_ss_ml_wit(~ml_ind,:) = 0;
        ind_ss_ml_wit(:,~ss_ind) = 0;
        wit_inter_sess_ss_ml{counter_sess} = xprob_corrected_(ind_ss_ml_wit(:) == 1,:);
        wit_inter_sess_ss_ml_ids{counter_sess} = ids_mat_(ind_ss_ml_wit(:) == 1,:);

        ind_ss_ml_bet = bet_ind_;
        ind_ss_ml_bet(~ml_ind,:) = 0;
        ind_ss_ml_bet(:,~ss_ind) = 0;
        bet_inter_sess_ss_ml{counter_sess} = xprob_corrected_(ind_ss_ml_bet(:) == 1,:);
        bet_inter_sess_ss_ml_ids{counter_sess} = ids_mat_(ind_ss_ml_bet(:) == 1,:);

        ind_ss_ml2_wit = wit_ind_;
        ind_ss_ml2_wit(~ml2_ind,:) = 0;
        ind_ss_ml2_wit(:,~ss_ind) = 0;
        wit_inter_sess_ss_ml2{counter_sess} = xprob_corrected_(ind_ss_ml2_wit(:) == 1,:);
        wit_inter_sess_ss_ml2_ids{counter_sess} = ids_mat_(ind_ss_ml2_wit(:) == 1,:);

        ind_ss_ml2_bet = bet_ind_;
        ind_ss_ml2_bet(~ml2_ind,:) = 0;
        ind_ss_ml2_bet(:,~ss_ind) = 0;
        bet_inter_sess_ss_ml2{counter_sess} = xprob_corrected_(ind_ss_ml2_bet(:) == 1,:);
        bet_inter_sess_ss_ml2_ids{counter_sess} = ids_mat_(ind_ss_ml2_bet(:) == 1,:);

        ind_ml_ml2_wit = wit_ind_;
        ind_ml_ml2_wit(~ml2_ind,:) = 0;
        ind_ml_ml2_wit(:,~ml_ind) = 0;
        wit_inter_sess_ml_ml2{counter_sess} = xprob_corrected_(ind_ml_ml2_wit(:) == 1,:);
        wit_inter_sess_ml_ml2_ids{counter_sess} = ids_mat_(ind_ml_ml2_wit(:) == 1,:);

        ind_ml_ml2_bet = bet_ind_;
        ind_ml_ml2_bet(~ml2_ind,:) = 0;
        ind_ml_ml2_bet(:,~ml_ind) = 0;
        bet_inter_sess_ml_ml2{counter_sess} = xprob_corrected_(ind_ml_ml2_bet(:) == 1,:);
        bet_inter_sess_ml_ml2_ids{counter_sess} = ids_mat_(ind_ml_ml2_bet(:) == 1,:);

        ind_ml_cs_wit = wit_ind_;
        ind_ml_cs_wit(~cs_ind,:) = 0;
        ind_ml_cs_wit(:,~ml_ind) = 0;
        wit_inter_sess_ml_cs{counter_sess} = xprob_corrected_cs_ss_(ind_ml_cs_wit(:) == 1,:);
        wit_inter_sess_ml_cs_ids{counter_sess} = ids_mat_(ind_ml_cs_wit(:) == 1,:);

        ind_ml_cs_bet = bet_ind_;
        ind_ml_cs_bet(~cs_ind,:) = 0;
        ind_ml_cs_bet(:,~ml_ind) = 0;
        bet_inter_sess_ml_cs{counter_sess} = xprob_corrected_cs_ss_(ind_ml_cs_bet(:) == 1,:);
        bet_inter_sess_ml_cs_ids{counter_sess} = ids_mat_(ind_ml_cs_bet(:) == 1,:);

        ind_ml2_cs_wit = wit_ind_;
        ind_ml2_cs_wit(~cs_ind,:) = 0;
        ind_ml2_cs_wit(:,~ml2_ind) = 0;
        wit_inter_sess_ml2_cs{counter_sess} = xprob_corrected_cs_ss_(ind_ml2_cs_wit(:) == 1,:);
        wit_inter_sess_ml2_cs_ids{counter_sess} = ids_mat_(ind_ml2_cs_wit(:) == 1,:);

        ind_ml2_cs_bet = bet_ind_;
        ind_ml2_cs_bet(~cs_ind,:) = 0;
        ind_ml2_cs_bet(:,~ml2_ind) = 0;
        bet_inter_sess_ml2_cs{counter_sess} = xprob_corrected_cs_ss_(ind_ml2_cs_bet(:) == 1,:);
        bet_inter_sess_ml2_cs_ids{counter_sess} = ids_mat_(ind_ml2_cs_bet(:) == 1,:);

        ind_ss_cs_wit = wit_ind_;
        ind_ss_cs_wit(~cs_ind,:) = 0;
        ind_ss_cs_wit(:,~ss_ind) = 0;
        wit_inter_sess_ss_cs{counter_sess} = xprob_corrected_cs_ss_(ind_ss_cs_wit(:) == 1,:);
        wit_inter_sess_ss_cs_ids{counter_sess} = ids_mat_(ind_ss_cs_wit(:) == 1,:);

        ind_ss_cs_bet = bet_ind_;
        ind_ss_cs_bet(~cs_ind,:) = 0;
        ind_ss_cs_bet(:,~ss_ind) = 0;
        bet_inter_sess_ss_cs{counter_sess} = xprob_corrected_cs_ss_(ind_ss_cs_bet(:) == 1,:);
        bet_inter_sess_ss_cs_ids{counter_sess} = ids_mat_(ind_ss_cs_bet(:) == 1,:);
    end
    wit_inter{1}  = [wit_inter{1}; cell2mat(wit_inter_sess_ss)];
    wit_inter{2}  = [wit_inter{2}; cell2mat(wit_inter_sess_ml)];
    wit_inter{3}  = [wit_inter{3}; cell2mat(wit_inter_sess_ml2)];
    wit_inter{4}  = [wit_inter{4}; cell2mat(wit_inter_sess_cs)];
    wit_inter{5}  = [wit_inter{5}; cell2mat(wit_inter_sess_ss_ml)];
    wit_inter{6}  = [wit_inter{6}; cell2mat(wit_inter_sess_ss_ml2)];
    wit_inter{7}  = [wit_inter{7}; cell2mat(wit_inter_sess_ml_ml2)];
    wit_inter{8}  = [wit_inter{8}; cell2mat(wit_inter_sess_ml_cs)];
    wit_inter{9}  = [wit_inter{9}; cell2mat(wit_inter_sess_ml2_cs)];
    wit_inter{10} = [wit_inter{10}; cell2mat(wit_inter_sess_ss_cs)];

    bet_inter{1}  = [bet_inter{1}; cell2mat(bet_inter_sess_ss)];
    bet_inter{2}  = [bet_inter{2}; cell2mat(bet_inter_sess_ml)];
    bet_inter{3}  = [bet_inter{3}; cell2mat(bet_inter_sess_ml2)];
    bet_inter{4}  = [bet_inter{4}; cell2mat(bet_inter_sess_cs)];
    bet_inter{5}  = [bet_inter{5}; cell2mat(bet_inter_sess_ss_ml)];
    bet_inter{6}  = [bet_inter{6}; cell2mat(bet_inter_sess_ss_ml2)];
    bet_inter{7}  = [bet_inter{7}; cell2mat(bet_inter_sess_ml_ml2)];
    bet_inter{8}  = [bet_inter{8}; cell2mat(bet_inter_sess_ml_cs)];
    bet_inter{9}  = [bet_inter{9}; cell2mat(bet_inter_sess_ml2_cs)];
    bet_inter{10} = [bet_inter{10}; cell2mat(bet_inter_sess_ss_cs)];

    wit_inter_ids{1}  = [wit_inter_ids{1}; vertcat(wit_inter_sess_ss_ids{:})];
    wit_inter_ids{2}  = [wit_inter_ids{2}; vertcat(wit_inter_sess_ml_ids{:})];
    wit_inter_ids{3}  = [wit_inter_ids{3}; vertcat(wit_inter_sess_ml2_ids{:})];
    wit_inter_ids{4}  = [wit_inter_ids{4}; vertcat(wit_inter_sess_cs_ids{:})];
    wit_inter_ids{5}  = [wit_inter_ids{5}; vertcat(wit_inter_sess_ss_ml_ids{:})];
    wit_inter_ids{6}  = [wit_inter_ids{6}; vertcat(wit_inter_sess_ss_ml2_ids{:})];
    wit_inter_ids{7}  = [wit_inter_ids{7}; vertcat(wit_inter_sess_ml_ml2_ids{:})];
    wit_inter_ids{8}  = [wit_inter_ids{8}; vertcat(wit_inter_sess_ml_cs_ids{:})];
    wit_inter_ids{9}  = [wit_inter_ids{9}; vertcat(wit_inter_sess_ml2_cs_ids{:})];
    wit_inter_ids{10} = [wit_inter_ids{10}; vertcat(wit_inter_sess_ss_cs_ids{:})];

    bet_inter_ids{1}  = [bet_inter_ids{1}; vertcat(bet_inter_sess_ss_ids{:})];
    bet_inter_ids{2}  = [bet_inter_ids{2}; vertcat(bet_inter_sess_ml_ids{:})];
    bet_inter_ids{3}  = [bet_inter_ids{3}; vertcat(bet_inter_sess_ml2_ids{:})];
    bet_inter_ids{4}  = [bet_inter_ids{4}; vertcat(bet_inter_sess_cs_ids{:})];
    bet_inter_ids{5}  = [bet_inter_ids{5}; vertcat(bet_inter_sess_ss_ml_ids{:})];
    bet_inter_ids{6}  = [bet_inter_ids{6}; vertcat(bet_inter_sess_ss_ml2_ids{:})];
    bet_inter_ids{7}  = [bet_inter_ids{7}; vertcat(bet_inter_sess_ml_ml2_ids{:})];
    bet_inter_ids{8}  = [bet_inter_ids{8}; vertcat(bet_inter_sess_ml_cs_ids{:})];
    bet_inter_ids{9}  = [bet_inter_ids{9}; vertcat(bet_inter_sess_ml2_cs_ids{:})];
    bet_inter_ids{10} = [bet_inter_ids{10}; vertcat(bet_inter_sess_ss_cs_ids{:})];
end

save_path = 'C:\Users\Jafar\Documents\reward\';
file_name = 'clique_interactions';
save(fullfile(save_path,'population_data',sprintf('%s.mat',file_name)), 'bet_inter','wit_inter','bet_inter_ids','wit_inter_ids', '-v7.3');

end