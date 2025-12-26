function plot_clique_interaction(data_path)

load(fullfile(data_path,'population_data','clique_interactions.mat'),'bet_inter','wit_inter','wit_inter_ids');

title_list = {'SS|SS','MLI|MLI','MLI2|MLI2','CS|CS','SS|MLI',...
    'SS|MLI2','MLI|MLI2','MLI|CS','MLI2|CS','SS|CS'};
delay_list = {0,0,0,0,1,3,2,0:10,0:10,1};

t_span = -50:50;

fig = figure;
for counter_ = 1:numel(title_list)-2
    subplot(2,5,counter_)
    sig_wit = wit_inter{counter_};
    cell_ids_wit = wit_inter_ids{counter_};
    rm_ind = any(isnan(sig_wit),2);
    sig_wit(rm_ind,:) = [];
    cell_ids_wit(rm_ind,:) = [];
    n_wit = size(sig_wit,1);
    [hl1,~] = boundedline(t_span,mean(sig_wit),2*std(sig_wit)/sqrt(n_wit),'r','alpha');
    hold on
    sig_bet = bet_inter{counter_};
    sig_bet(any(isnan(sig_bet),2),:) = [];
    n_bet = size(sig_bet,1);
    [hl2,~] = boundedline(t_span,mean(sig_bet),2*std(sig_bet)/sqrt(n_bet),'b','alpha');
    if counter_<= 4
        legend([hl1,hl2],{['wit: ', num2str(n_wit/2)], ['bet:', num2str(n_bet/2)]},'Box','off');
    else
        legend([hl1,hl2],{['wit: ', num2str(n_wit)], ['bet:', num2str(n_bet)]},'Box','off');
    end
    xlim([-25,25])
    xline(0,'--','HandleVisibility','off')
    yline(0,'--','HandleVisibility','off')

    sig_bet_ = mean(sig_bet(:,delay_list{counter_} + 51),2);
    sig_wit_ = mean(sig_wit(:,delay_list{counter_} + 51),2);

    [~,p] = ttest2(sig_bet_,sig_wit_);
    title({title_list{counter_},[num2str(mean(sig_bet_),'%.2f'),'±',num2str(std(sig_bet_)/sqrt(n_bet),'%.2f'),', '...
        num2str(mean(sig_wit_),'%.2f'),'±',num2str(std(sig_wit_,1)/sqrt(n_wit),'%.2f'),', '...
        '(p = ', num2str(p,'%.2e'),')']});

end

ESN_Beautify_Plot(fig,[16,8])
save_path1 = fullfile(data_path,'population_figs','ephys','clique-interaction.pdf');
print(fig,save_path1,'-dpdf','-bestfit')
end