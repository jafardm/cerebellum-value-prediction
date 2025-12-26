function [chanMapFile,probe_flag,Cat_answer] = JDM_probe_channel_map

electList = {'M1', 'M2', 'npx', 'RHS'};
answer_ = listdlg('PromptString', 'Select Electrode:', ...
    'SelectionMode', 'single', ...
    'ListString', {'M1', 'M2', 'npx', 'RHS'});

answer_elec = electList{answer_};

switch(answer_elec)
    case('M1')
        answer_cmap = questdlg('Which Channel Map?','Channel Map',...
            '0: Checker pre june 2022','2: Checker post june 2022',...
            '2: Checker post june_2022');
        Cat_answer = questdlg('Would you like concatenate data?','Concatenate','1:Yes','0:No','1:Yes');
        Cat_answer = str2double(Cat_answer(1));
    case('M2')
        answer_cmap = questdlg('Which Channel Map','Channel Map',...
            '1: Linear pre june 2022','3: Linear post june_2022',...
            '1: Linear pre june 2022');
        Cat_answer = questdlg('Would you like concatenate the data?','Concatenate','1:Yes','0:No','1:Yes');
        Cat_answer = str2double(Cat_answer(1));
    case('npx')
        answer_cmap = questdlg('Which Channel Map','Channel Map',...
            '4: npx 1.0 staggered OE','5: SGLX ','5: SGLX ');
        Cat_answer = questdlg('What version data?','Recorded Format:','1:supercat','0:single','1:supercat');
        Cat_answer = str2double(Cat_answer(1));
    case('RHS')
        answer_cmap = questdlg('Which Channel Map?','Channel Map',...
            '6: Checker post june 2022','7: linear post june_2022',...
            '6: Checker post june 2022');
         Cat_answer = questdlg('What version data?','Recorded Format:','1:supercat','0:single','1:supercat');
         Cat_answer = str2double(Cat_answer(1));
        
end

 geom = str2double(answer_cmap(1));

switch(geom)
    % 0: checkerboard M1
    case(0)
        chanMapFile = 'Camb_Check_kilosortChanMap.mat';
        probe_flag = 0;
    % 1: linear M2
    case(1)
        chanMapFile = 'Camb_lin_kilosortChanMap.mat';
        probe_flag = 1;
    % 2: checkerboard M1 post 0622
    case(2)
        chanMapFile = 'Camb_Check_kilosortChanMap.mat';
        probe_flag = 2;
    % npx checker
    case(4)
        probe_flag = 4;
        chanMapFile = 'neuropixPhase3B1_kilosortChanMap.mat';
    % npx SGLX 
    case(5)
        probe_flag = 5;
        chanMapFile = 'neuropixPhase3B2_kilosortChanMap.mat';
    case(6)
        probe_flag = 6;
        chanMapFile = 'Camb_Check_kilosortChanMap.mat';
    case (7)
        probe_flag = 7;
        chanMapFile = 'Camb_lin_kilosortChanMap.mat';
end

