function JDM_build_neural_properties(data_path)
JDM_params_funcs

user = 'JDM';
loc  = 'ctx';

path_data_monkey_sorted = params.path_data_monkey_sorted(1:2);
animal_list = params.animal_list(1:2);
session_list = MAF_extract_sess_list(path_data_monkey_sorted,user,loc);

num_animal = numel(session_list);

% SS, CS, PC, MLI, MLI2
wave_tot_ss  = [];
wave_tot_cs  = [];
acorr_tot_ss = [];
acorr_tot_cs = [];
xcorr_tot    = [];
rate_bl_ss   = [];
rate_bl_cs   = [];
cell_types   = cell(0);
cell_ids     = cell(0);


for counter_animal = 1:num_animal
    disp(animal_list{counter_animal});
    current_path = path_data_monkey_sorted{counter_animal};
    current_sess_list = session_list{counter_animal};
    num_sess = numel(current_sess_list);

    waveform_ss  = cell(num_sess,1);
    waveform_cs  = cell(num_sess,1);
    acorr_ss     = cell(num_sess,1);
    acorr_cs     = cell(num_sess,1);
    xcorr        = cell(num_sess,1);
    cell_type    = cell(num_sess,1);
    cell_id      = cell(num_sess,1);
    bl_ss        = cell(num_sess,1);
    bl_cs        = cell(num_sess,1);

    parfor counter_sess = 1:num_sess
        current_sess = current_sess_list{counter_sess};
        disp(['     ' current_sess]);

        units_info = MAF_extract_cell_metadata(current_path,current_sess);
        
        cell_ids_     = units_info.cell_list;
        cell_types_   = units_info.cell_type;

        cell_ind = ismember(cell_types_,{'MLI','MLI2','PC','CS','SS'});

        res = par_loop_neural_prop(cell_ids_(cell_ind),cell_types_(cell_ind),current_path,current_sess);
        
        waveform_ss{counter_sess}  = res.wave_tot_ss;
        waveform_cs{counter_sess}  = res.wave_tot_cs;
        acorr_ss{counter_sess}     = res.acorr_tot_ss;
        acorr_cs{counter_sess}     = res.acorr_tot_cs;
        xcorr{counter_sess}        = res.xcorr_tot;
        cell_type{counter_sess}    = res.cell_types;
        cell_id{counter_sess}      = res.cell_ids;
        bl_ss{counter_sess}        = res.bl_ss;
        bl_cs{counter_sess}        = res.bl_cs;
    end
    wave_tot_ss  = [wave_tot_ss;cell2mat(waveform_ss)];
    wave_tot_cs  = [wave_tot_cs;cell2mat(waveform_cs)];
    acorr_tot_ss = [acorr_tot_ss;cell2mat(acorr_ss)];
    acorr_tot_cs = [acorr_tot_cs;cell2mat(acorr_cs)];
    xcorr_tot    = [xcorr_tot;cell2mat(xcorr)];
    rate_bl_ss   = [rate_bl_ss;cell2mat(bl_ss)];
    rate_bl_cs   = [rate_bl_cs;cell2mat(bl_cs)];
    cell_types   = [cell_types;vertcat(cell_type{:})];
    cell_ids     = [cell_ids;vertcat(cell_id{:})];
end

save(fullfile(data_path, 'population_data','neural_prop_clique.mat'),...
    'wave_tot_ss','wave_tot_cs',...
    'acorr_tot_ss','acorr_tot_cs','xcorr_tot','cell_types','cell_ids',...
    'rate_bl_cs',"rate_bl_ss");

end

function res = par_loop_neural_prop(cell_ids,cell_types,current_path,current_sess)

num_cell = length(cell_types);
len_x = 50;
len_y = 150;

wave_tot_ss  = nan(num_cell,len_x*len_y);
wave_tot_cs  = nan(num_cell,len_x*len_y);
acorr_tot_ss = nan(num_cell,101);
acorr_tot_cs = nan(num_cell,201);
xcorr_tot    = nan(num_cell,101);
bl_ss        = nan(num_cell,1);
bl_cs        = nan(num_cell,1);

parfor counter_cell = 1:num_cell
    current_cell = cell_ids{counter_cell};
    current_type = cell_types{counter_cell};
    data = MAF_load_cell(current_path,current_sess,current_cell);
    data_recordings = data.data_recordings;
    Neural_Prop = MAF_combineNeuralProp([data_recordings.Neural_Prop]);
    waveform = Neural_Prop.waveform;
    if strcmp(current_type,'CS')
        wave_ = MAF_waveform_to_image(waveform, len_x, len_y, 'CS');
        wave_tot_cs(counter_cell,:) = wave_(:);
        acorr_tot_cs(counter_cell,:) = Neural_Prop.cs_aprob.xprob;

        sp_time = Neural_Prop.CS_time;
        sp_ISI = diff(sp_time);
        sp_ISI(abs(sp_ISI)>50)=0;
        bl_cs(counter_cell,:) = (length(sp_ISI)+1) ./ sum(sp_ISI); % 1ms probability

    elseif strcmp(current_type,'PC')
        wave_ = MAF_waveform_to_image(waveform, len_x, len_y, 'SS');
        wave_tot_ss(counter_cell,:) = wave_(:);
        acorr_tot_ss(counter_cell,:) = Neural_Prop.ss_xprob.xprob;
        sp_time = Neural_Prop.SS_time;
        sp_ISI = diff(sp_time);
        sp_ISI(abs(sp_ISI)>5)=0;
        bl_ss(counter_cell,:) = (length(sp_ISI)+1) ./ sum(sp_ISI); % 1ms probability

        wave_ = MAF_waveform_to_image(waveform, len_x, len_y, 'CS');
        wave_tot_cs(counter_cell,:) = wave_(:);
        acorr_tot_cs(counter_cell,:) = Neural_Prop.cs_aprob.xprob;
        sp_time = Neural_Prop.CS_time;
        sp_ISI = diff(sp_time);
        sp_ISI(abs(sp_ISI)>50)=0;
        bl_cs(counter_cell,:) = (length(sp_ISI)+1) ./ sum(sp_ISI); % 1ms probability

        xcorr_tot(counter_cell,:) = Neural_Prop.cs_xprob.xprob;
    else
        wave_ = MAF_waveform_to_image(waveform, len_x, len_y, 'SS');
        wave_tot_ss(counter_cell,:) = wave_(:);
        acorr_tot_ss(counter_cell,:) = Neural_Prop.ss_xprob.xprob;
        sp_time = Neural_Prop.SS_time;
        sp_ISI = diff(sp_time);
        sp_ISI(abs(sp_ISI)>5)=0;
        bl_ss(counter_cell,:) = (length(sp_ISI)+1) ./ sum(sp_ISI); % 1ms probability
    end
end

res.wave_tot_ss  = wave_tot_ss;
res.wave_tot_cs  = wave_tot_cs;
res.acorr_tot_ss = acorr_tot_ss;
res.acorr_tot_cs = acorr_tot_cs;
res.xcorr_tot    = xcorr_tot;
res.cell_types   = cell_types;
res.cell_ids     = cell_ids;
res.bl_ss        = bl_ss;
res.bl_cs        = bl_cs;

end