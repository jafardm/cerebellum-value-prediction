function plot_saccade_counts(data_path)

    % Dependencies & data
    JDM_params_funcs;
    S = load(fullfile(data_path, 'population_data', 'behave_data.mat'));  % loads behave_data, params, etc.
    if isfield(S,'data_behav'), behave_data = S.data_behav; else, error('behave_data missing'); end
 
    % Output folder
    img_save_path = fullfile('C:\Users\Jafar\Documents\reward\population_figs\behavior');
    if ~exist(img_save_path,'dir'), mkdir(img_save_path); end

    % Subjects
    subj_names   = params.animal_list;
    num_subjects = numel(subj_names);

    % -------------------------------
    % Condition definitions & labels
    % -------------------------------
    % [idx_tsk ,idx_tag, idx_tgt, idx_jmp, idx_rew]
    forced_conds = [ 1 1 1 0 1;  % HH
                     1 1 1 0 0;  % HL
                     1 1 0 0 0;  % LL
                     1 1 0 0 1;  % LH
                     1 1 1 1 1;  % HHJ
                     1 1 1 1 0;  % HLJ
                     1 1 0 1 0;  % LLJ
                     1 1 0 1 1]; % LHJ
    labels_forced = {'HH','HL','LL','LH','HHJ','HLJ','LLJ','LHJ'};

    % [idx_tsk ,idx_tag, idx_choice]  (choice task has idx_tsk = 0)
    choice_conds  = [0 1 1;  % High
                     0 1 0]; % Low
    labels_choice = {'High','Low'};


    % Colors
    col_forced = lines(numel(labels_forced));
    col_choice = lines(numel(labels_choice));

    % ------------------------
    % Iterate through subjects
    % ------------------------
    for jj = 1:num_subjects
        subj = subj_names{jj};
    
        % Pull required index vectors (must be 0/1 vectors per trial)
        idx_tsk   = behave_data.task_cond{jj};
        idx_tag   = behave_data.sac_tag{jj};
        idx_tgt   = behave_data.tgt_cond{jj};
        idx_jmp   = behave_data.jump_cond{jj};
        idx_rew   = behave_data.rew_cond{jj};
        % For choice:
        idx_choice = behave_data.choice{jj};

        istuned_HH = ismember([idx_tsk idx_tag idx_tgt idx_jmp idx_rew],...
            forced_conds(1,:),'rows');
        istuned_HL = ismember([idx_tsk idx_tag idx_tgt idx_jmp idx_rew],...
            forced_conds(2,:),'rows');
        istuned_LL = ismember([idx_tsk idx_tag idx_tgt idx_jmp idx_rew],...
            forced_conds(3,:),'rows');
        istuned_LH = ismember([idx_tsk idx_tag idx_tgt idx_jmp idx_rew],...
            forced_conds(4,:),'rows');
    
        istuned_HHJ = ismember([idx_tsk idx_tag idx_tgt idx_jmp idx_rew],...
            forced_conds(5,:),'rows');
        istuned_HLJ = ismember([idx_tsk idx_tag idx_tgt idx_jmp idx_rew],...
            forced_conds(6,:),'rows');
        istuned_LLJ = ismember([idx_tsk idx_tag idx_tgt idx_jmp idx_rew],...
            forced_conds(7,:),'rows');
        istuned_LHJ = ismember([idx_tsk idx_tag idx_tgt idx_jmp idx_rew],...
            forced_conds(8,:),'rows');
        
       istuned_choice_high = ismember([idx_tsk idx_tag idx_choice ],...
            choice_conds(1,:),'rows');
       istuned_choice_low = ismember([idx_tsk idx_tag idx_choice ],...
            choice_conds(2,:),'rows');

           % -------------------------
    % Count & Plot per subject
    % -------------------------
    N_forced = [ ...
        sum(istuned_HH);  sum(istuned_HL);  sum(istuned_LL);  sum(istuned_LH); ...
        sum(istuned_HHJ); sum(istuned_HLJ); sum(istuned_LLJ); sum(istuned_LHJ) ...
    ];

    N_choice = [ ...
        sum(istuned_choice_high); ...
        sum(istuned_choice_low)  ...
    ];

    % Pretty labels with n
  

    fig = figure('Name', sprintf('Saccade counts — %s', subj), ...
                 'Color','w', 'Position',[100 100 1150 420]);
    tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

    % ---- Forced ----
    ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');
    bh1 = bar(ax1, N_forced, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1);
    for i = 1:numel(N_forced), bh1.CData(i,:) = col_forced(i,:); end
    text(1:numel(N_forced), N_forced, compose('%d',N_forced), ...
         'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
    set(ax1,'XTick',1:numel(labels_forced), 'XTickLabel',labels_forced, ...
            'TickDir','out','LineWidth',1);
    ylabel(ax1,'# Saccades');
    title(ax1, sprintf('%s — Forced', subj), 'FontWeight','bold');

    % ---- Choice ----
    ax2 = nexttile; hold(ax2,'on'); box(ax2,'on');
    bh2 = bar(ax2, N_choice, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1);
    for i = 1:numel(N_choice), bh2.CData(i,:) = col_choice(i,:); end
    text(1:numel(N_choice), N_choice, compose('%d',N_choice), ...
         'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
    set(ax2,'XTick',1:numel(labels_choice), 'XTickLabel',labels_choice, ...
            'TickDir','out','LineWidth',1);
    ylabel(ax2,'# Saccades');
    title(ax2, sprintf('%s — Choice', subj), 'FontWeight','bold');

    title(tl, 'Saccade Counts by Condition', 'FontWeight','bold');

    % Save
    base = sprintf('saccade_counts_%s', subj);
    ESN_Beautify_Plot(fig,[10 4])
    print(fig, fullfile(img_save_path, [base '.pdf']), '-dpdf', '-bestfit');
    close(fig);

    end

end

