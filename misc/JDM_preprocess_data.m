%% LOAD CS POPULATION METADATA

data_path      = 'C:\Users\Jafar\Documents\new';
save_file_path = 'C:\Users\Jafar\Documents\new\data_population\';
metadata_dir   = dir2(fullfile(data_path,'*id_new*'));
pcell_path     = 'C:\Users\Jafar\Documents\new\ALL_PCELL\cell_data';
%% variables
JSP_global_variables(1)
global event_type_list amp_edges vel_edges length_trace tag_name_list min_rxn_time max_rxn_time ang_values ang_edges inds_span

num_tag_bin   = length(tag_name_list);
num_ang_bin   = length(ang_values);
num_amp_bin   = length(amp_edges) - 1;
num_vel_bin   = length(vel_edges) - 1;
num_event_bin = length(event_type_list);

time_win              = [-90 10];
tag_                  = 1;
min_sac            = 300;
tag_local_list        = 1; % for building pair data, which tags to look at
event_type_local_list = {'onset'};
num_local_event_bin   = length(event_type_local_list);
aligh_event_type      = 'vmax'; % vamx or onsst
local_amp_edges       = [0,9];
num_local_amp_bin     = length(local_amp_edges) - 1;
data_type_local_list  = {'CS'}; 

%% METADATA & CELL LIST

% Get the cells list
metadata  = readtable(fullfile(metadata_dir.folder, metadata_dir.name)); 
idx_eye   = logical(metadata.eye_mod_iter_2);
trail_num = metadata.num_trial >= min_sac;
idx_selected_cells = logical(trail_num .* idx_eye);
cell_list = metadata.id(idx_selected_cells); % take the eye-modulated neurons

num_pCells = size(cell_list,1);
%% Extract 125d data

monkeyname = '125d';
idx_monkey = cell2mat(cellfun(@(x) strcmp(x,monkeyname),table2cell(metadata(:,1)),'UniformOutput',false));
cell_list  = table2cell(metadata(idx_monkey,:));


for icell = 1:size(cell_list,1)

    fprintf('***LOADING PCELL NUMBER %d ***\n' ,icell)

   load(fullfile(cell_list{icell,3},cell_list{icell,2}));
   CS_on_data = JSP_CS_on_analysis_fr(SACS_ALL_DATA,0);
   file_name = fullfile(save_file_path,cell_list{icell,2});
   save(file_name,'EXPERIMENT_PARAMS','SACS_ALL_DATA','Neural_Properties','CS_on_data',"-v7.3")
   clearvars EXPERIMENT_PARAMS SACS_ALL_DATA Neural_Properties CS_on_data

   fprintf('***PCELL NUMBER %d DONE***\n' ,icell)
end
%% Create metadata (replicate Salomon's paper)

taglist = 1;
data_all_tags_SS = struct([]);

for tag_idx = 1:length(taglist)

    for icell = 1:num_pCells

        load(fullfile(data_path,'ALL_PCELL\cell_data\',cell_list{icell}));

        tag_ = SACS_ALL_DATA.tag ==tag_idx;

        CS_ang = CS_on_data.CS_ang_avg;
        CS_ang_round = round(CS_on_data.CS_ang_avg/45)*45;

        if ~isfield(SACS_ALL_DATA,'eye_r_px_onset') 
            SACS_ALL_DATA.eye_r_px_onset  = SACS_ALL_DATA.eye_l_px_onset;
            SACS_ALL_DATA.eye_r_px_offset = SACS_ALL_DATA.eye_l_px_offset;
            SACS_ALL_DATA.eye_r_py_onset  = SACS_ALL_DATA.eye_l_py_onset;
            SACS_ALL_DATA.eye_r_py_offset = SACS_ALL_DATA.eye_l_py_offset;
        end

        data_all_tags_SS(tag_idx).data_all(icell).CS_ang           = CS_ang;
        data_all_tags_SS(tag_idx).data_all(icell).CS_ang_round     = CS_ang_round;
        data_all_tags_SS(tag_idx).data_all(icell).neuro_CS_onset   = SACS_ALL_DATA.neuro_CS_onset(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).neuro_CS_visual  = SACS_ALL_DATA.neuro_CS_visual(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).neuro_CS_visual  = SACS_ALL_DATA.neuro_CS_visual(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).neuro_SS_onset   = SACS_ALL_DATA.neuro_SS_onset(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).neuro_SS_visual  = SACS_ALL_DATA.neuro_SS_visual(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).eye_r_px_onset   = SACS_ALL_DATA.eye_r_px_onset(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).eye_r_px_offset  = SACS_ALL_DATA.eye_r_px_offset(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).visual_px_onset  = SACS_ALL_DATA.visual_px_onset(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).visual_px_offset = SACS_ALL_DATA.visual_px_offset(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).eye_r_py_onset   = SACS_ALL_DATA.eye_r_py_onset(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).eye_r_py_offset  = SACS_ALL_DATA.eye_r_py_offset(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).visual_py_onset  = SACS_ALL_DATA.visual_py_onset(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).visual_py_offset = SACS_ALL_DATA.visual_py_offset(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).eye_vx_vmax      = SACS_ALL_DATA.eye_vx_vmax(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).eye_vy_vmax      = SACS_ALL_DATA.eye_vy_vmax(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).reaction_time    = SACS_ALL_DATA.reaction(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).sac_time_offset  = SACS_ALL_DATA.time_offset(:,tag_);
        data_all_tags_SS(tag_idx).data_all(icell).nr_success_sac   = sum([SACS_ALL_DATA.tag ==1 SACS_ALL_DATA.tag ==4 SACS_ALL_DATA.tag ==6]);

         clearvars EXPERIMENT_PARAMS SACS_ALL_DATA Neural_Properties CS_on_data

        fprintf('***PCELL NUMBER %d FOR TAG %d DONE***\n' ,icell,tag_idx)

    end
    
end

save_name = fullfile(save_file_path,'data_all_tags_SS');
save(save_name,'data_all_tags_SS',"-v7.3")


%% MANUAL CATEGORIZATION OF CS WAVEFORM TYPE
% Need to have built CS population data first
% Manually initialize 'cell_type' group in metadata as '4'
% Manually update metadata at the end
% 3 Classification:
%   1: early
%   2: early & late (not used)
%   3: late
%   4: unclassified (e.g., not eye-mod.)

% 4 Classification:
%   1: early
%   2: weak early
%   3: weak late
%   4: late
%   5: unclassified (e.g., not eye-mod.)
visual_resp_period = -24:250; % period of CS resp., wrt. visual alignment, to be used to categorize cell types
early_range = 1:100; 
late_range  = 101:200;
unclassifed_num = 5;
% Data smoothing parameters
smooth_params.num_points           = 31;
smooth_params.method               = 'sgolay';

% Read single and paired metadata
metadata_dir = dir2(fullfile(parent_path,'*xls*'));
metadata_single_id   = readtable(fullfile(metadata_dir.folder, metadata_dir.name),'Sheet','id');
metadata_single_list = readtable(fullfile(metadata_dir.folder, metadata_dir.name),'Sheet','list');
metadata_single_list = table2cell(metadata_single_list);
% Processing the list
for i = 1 : size(metadata_single_list,1)
    for j = 1 : size(metadata_single_list,2)
        if contains(metadata_single_list{i,j},'[]')
            metadata_single_list{i,j} = [];
        end
    end
end
metadata_pair_id     = readtable(fullfile(metadata_dir.folder, metadata_dir.name),'Sheet','pair_id');
metadata_pair_list   = readtable(fullfile(metadata_dir.folder, metadata_dir.name),'Sheet','pair_list');
metadata_pair_list   = table2cell(metadata_pair_list);
pCell_ids            = metadata_single_id.id;
pCell_list_isstr     = arrayfun(@iscellstr,metadata_single_list);

% Specify cells to analyze
cell_to_analyze = find(metadata_single_id.num_trial > 300 & metadata_single_id.eye_mod_iter_2 == 1);

% Load CS population data
num_pCells      = length(cell_to_analyze);
file_name = sprintf('CS_population_tuned_%d_cell',num_pCells);
loaded_data      = load(fullfile(parent_path,'data_population',sprintf('%s.mat',file_name)),'population_tuned', 'num_sac_tuned', 'population_neural_properties');
population_tuned = loaded_data.population_tuned;
CS_resp_orig     = population_tuned.amp(1).visual{1,1}.*1000;
CS_resp_smoothed = JSP_smooth(CS_resp_orig, 2, smooth_params);
CS_resp          = CS_resp_smoothed(:,visual_resp_period+length_trace./2);

% Init. var
single_cell_type = ones(size(metadata_single_id.cell_type))*unclassifed_num;
pair_cell_type   = ones(size(metadata_pair_id.cell_type))*unclassifed_num;

% Semi-manually classify cell type
% Determine based on early and late phase response
figure;
resp_plot = plot([0],[0],'k','LineWidth',2); % init. plot
xline(0);
xlabel('Time from visual stimulus (ms)');
ylim([0,10]);
xlim([min(visual_resp_period), max(visual_resp_period)]);   
ylabel('CS firing (change, Hz');
for counter_pCell = 1 : num_pCells
    cell_idx = cell_to_analyze(counter_pCell);
    
    resp_plot.XData = visual_resp_period; resp_plot.YData = CS_resp(counter_pCell,:);
    
    early_resp = CS_resp_smoothed(counter_pCell, early_range+length_trace./2);
    early_peak = max(early_resp);
    early_area = sum(early_resp);
    late_resp = CS_resp_smoothed(counter_pCell, late_range+length_trace./2);
    late_peak = max(late_resp);
    late_area = sum(late_resp);

    early_vs_late_area = early_area > late_area;
    early_vs_late_peak = early_peak > late_peak;

    if early_vs_late_area == early_vs_late_peak
        resp_plot.Color = 'k';
        if early_vs_late_area == 1
            cell_type = 1;
        else
            cell_type = 3;
        end
    else
        resp_plot.Color = 'r';
        cell_type = 4;
    end
    title(sprintf('early: %.0f, late: %.0f, type: %d', early_area, late_area, cell_type));
    fprintf('\nCell # %d; %s,',counter_pCell, pCell_ids{cell_idx});
    cell_type = input('type: ');
    single_cell_type(cell_idx,1) = cell_type;

    % Find shared recordings in the pair list
    for rec_idx = find(pCell_list_isstr(cell_idx,:))
        rec_name = metadata_single_list{cell_idx, rec_idx};
        pair_cell_type(logical(sum(arrayfun(@(str_) contains(str_, rec_name), metadata_pair_list),2))) = cell_type;
    end  
end
%% 
for ii = 1:size(data_all_tags_SS.data_all,2)
   temp(ii) =  sum(sum(data_all_tags_SS.data_all(ii).neuro_SS_onset)) ; 
end
sum(temp>0)