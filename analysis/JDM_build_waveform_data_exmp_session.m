function JDM_build_waveform_data_exmp_session
JDM_params_funcs;

current_path = params.path_data_monkey_sorted{2};
current_sess = '2024-12-30';

PC_list  =  {'46_087','197_274','220_313',...
             '221_309','228_318','293_379','294_367'};

ML_list  = {'18_016','46_089','47_090' ,'48_091',... 
          '180_234','182_228','178_231','180_233',...
          '182_245', '195_271', '217_305', '224_314',...
          '224_315','287_374','287_383','286_372',...
          '284_363','287_366', '289_519','286_462',...
          '284_518','282_362'};

ML2_list = {'196_272','88_113','197_273'};
MF_list  = {'27_050','31_055','36_062',...
    '163_201','173_219','171_216',...
    '232_320','234_325','279_359',...
    '275_508'};



cells = struct();

for counter_cell = 1:numel(PC_list)
    current_cell = PC_list{counter_cell};
    cells.(['pc_' num2str(counter_cell)]) = [MAF_load_cell(current_path,current_sess,current_cell).data_recordings];
end

for counter_cell = 1:numel(ML_list)
    current_cell = ML_list{counter_cell};
    cells.(['mli1_' num2str(counter_cell)]) = [MAF_load_cell(current_path,current_sess,current_cell).data_recordings];
end

for counter_cell = 1:numel(ML2_list)
    current_cell = ML2_list{counter_cell};
    cells.(['mli2_' num2str(counter_cell)]) = [MAF_load_cell(current_path,current_sess,current_cell).data_recordings];
end

for counter_cell = 1:numel(MF_list)
    current_cell = MF_list{counter_cell};
    cells.(['mf_' num2str(counter_cell)]) = [MAF_load_cell(current_path,current_sess,current_cell).data_recordings];
end

cell_lists = {PC_list; ML_list; ML2_list; MF_list};  % vertical concat

% Generate corresponding names automatically
cell_list_name = [ ...
    arrayfun(@(x) sprintf('pc_%d', x), 1:numel(PC_list), 'UniformOutput', false), ...
    arrayfun(@(x) sprintf('mli1_%d', x), 1:numel(ML_list), 'UniformOutput', false), ...
    arrayfun(@(x) sprintf('mli2_%d', x), 1:numel(ML2_list), 'UniformOutput', false), ...
    arrayfun(@(x) sprintf('mf_%d', x), 1:numel(MF_list), 'UniformOutput', false) ...
];

path_save = 'C:\Users\Jafar\Documents\reward\population_data';
save(fullfile(path_save,'fig1_waveforms.mat'),"cell_list_name","cells","cell_lists");
end