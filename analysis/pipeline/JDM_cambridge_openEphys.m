function JDM_cambridge_openEphys (path_data_monkey_sorted,sess,probe_flag,catFlag)

 % 1: linear M2
id_linear = [28,26,24,19,21,15,32,30,18,29,31,22,20,23,25,27,1,3,6,8,10,12,14,16,17,13,11,...
            9,7,5,4,2,64,62,59,57,55,53,51,47,50,52,54,56,58,60,61,63,37,39,41,46,44,33,35,48,...
            36,34,49,43,45,42,40,38];

 % 0: checkerboard M1
id_checker = [28,26,24,19,21,15,32,30,18,29,31,22,20,23,25,27,1,3,6,8,10,12,14,16,17,13,11,...
            9,7,5,4,2,64,62,59,57,55,53,51,47,50,52,63,49,54,33,37,43,56,35,39,45,58,48,41,42,...
            60,36,46,40,61,34,44,38];

 % 2: checkerboard M1 post 0622
id_checker_post_0622 = [27,25,23,20,22,31,29,18,30,32,15,21,19,24,26,28,...
    2,4,5,7,9,11,13,17,16,14,12,10,8,6,3,1,64,62,59,57,55,53,51,47,50,...
    52,63,49,54,33,37,43,56,35,39,45,58,48,41,42,60,36,46,40,61,34,...
    44,38];


if probe_flag == 0
    geom = id_checker;
elseif probe_flag == 1
    geom = id_linear;
elseif probe_flag == 2
    geom = id_checker_post_0622;
elseif probe_flag == 4
    geom = nan;
else
    error('The Cambridge Channel map is nor correct!')
end

 if catFlag && (probe_flag == 0 || probe_flag == 1 || probe_flag == 2)
    % data path
    sess_path = [path_data_monkey_sorted, sess(1:7), filesep, sess, filesep];
    
    % number of recorded sessions
    rec_list = dir2(sess_path);
    rec_list([rec_list.isdir] == 0) = [];
    rec_name = {rec_list.name};
    rec_path = {rec_list.folder};
    
    % checking if there is binary data exist
    containsLetter = @(str) any(isletter(str));
    id_let = cellfun(containsLetter, rec_name);
    rec_name = rec_name(~id_let);
    rec_path = rec_path(~id_let);
    
    rec_num = length(rec_name);

    if isempty(rec_num)
        error("There is no Electropysiology data!")
    end

    % Initialize arrays for data and timestamps
    allData = [];

    smp_imap  = nan(rec_num,1);
    sec_imap  = nan(rec_num,1);

    for counter_rec = 1:rec_num
        pat_to_raw = fullfile(rec_path{counter_rec}, rec_name{counter_rec}, 'raw_data');
        if ~isfolder(pat_to_raw)
            error("There is no raw_data folder in %s", fullfile(rec_path{counter_rec}, rec_name{counter_rec}));
        end
           fprintf('')
        % get the ephys data    
        ch_data = extract_openEphys(pat_to_raw,geom);
    
            if sum(ch_data(:)) == 0
                fprintf('No Ephys \n');
                return;
            end
        
        allData = [allData, ch_data];
        smp_imap(counter_rec,1) = size(ch_data,2);
        sec_imap(counter_rec,1) = size(ch_data,2)/3e4; 
    end
    
    fprintf(' \n Concatenation of the data completed. Converting to Binary \n');
    
    cat_recs_path = fullfile(rec_path{counter_rec},'cat_recs\analyzed_data\sorted_data');
    if ~isfolder(cat_recs_path)
      mkdir(cat_recs_path);
    end
    
    data_name   = [sess '_tcat.camb0.ap'];
    offset_name = [sess '_camb_offsets'];
    
    fid = fopen([cat_recs_path filesep data_name '.bin'],'w');
    dat = int16(allData);
    fwrite(fid,dat,'int16');
    fclose(fid);
 
% creates metadata file

 % Initialize the first sample of the concatenated file to 0
smp_imap0 = zeros(length(smp_imap),1); 
smp_imap0(1) = 0;  % First sample is 0

% Calculate cumulative sums for sample indices
for ii = 2:length(smp_imap)
    smp_imap0(ii) = smp_imap0(ii-1) + smp_imap(ii-1);
end

% Corrected time offset initialization
sec_imap0 = zeros(length(sec_imap), 1); 
sec_imap0(1) = 0;

for ii = 2:length(sec_imap)
    sec_imap0(ii) = sec_imap0(ii-1) + sec_imap(ii-1);
end
  
% Sampling rate (assuming constant rate of 30kHz)
srate_imap0 = repmat(3e4, length(smp_imap0), 1);  
% Total recording time
fileTimeSecs = sum(sec_imap);
fileName = [cat_recs_path filesep offset_name '.txt'];

createMetadataTextfile(fileName,'smp_imap0', smp_imap0,...
    'sec_imap0',sec_imap0,'srate_imap0',srate_imap0,...
    'fileTimeSecs',fileTimeSecs);


fprintf('-> Created a Metadata of the recording.  \n');

 elseif not(catFlag) && (probe_flag == 0 || probe_flag == 1 || probe_flag == 2)
    % no concatenate 
    pattern_sess = '(.*\\\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})'; 
    sess_name = regexp(path_data_monkey_sorted, pattern_sess, 'match', 'once');
    
    pattern_rec = '\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}'; 
    rec_name = regexp(path_data_monkey_sorted, pattern_rec, 'match', 'once');

    ch_data = extract_openEphys(path_data_monkey_sorted,geom);

    if sum(ch_data(:)) == 0
        fprintf('No Ephys \n');
        return;
    end
       
    fprintf(' \n done. Converting to Binary \n');
    
    save_path = fullfile(sess_name,'\analyzed_data\sorted_data');
    if ~isfolder(save_path)
      mkdir(save_path);
    end
     
    fid = fopen([save_path filesep rec_name '.bin'],'w');
    dat = int16(ch_data);
    fwrite(fid,dat,'int16');
    fclose(fid);

 elseif probe_flag == 4 && catFlag == 0
     
   pattern_sess = '(.*\\\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})'; 
   sess_name = regexp(path_data_monkey_sorted, pattern_sess, 'match', 'once');
   save_path = fullfile(sess_name,'analyzed_data\sorted_data');
   % copy the dat file
    fprintf('copying file ')
    file_path = [path_data_monkey_sorted filesep 'experiment1' filesep...
        'recording1' filesep 'continuous' filesep ...
        'Neuropix-PXI-100.0' filesep 'continuous.dat'];
    mkdir(save_path)
    copyfile(file_path,[save_path filesep]);
 else
     error("Unkown probe type!")

 end

end
%% Get the data recorded with openEphys
function ch_data = extract_openEphys(pat_to_raw,id_)

nofCh = length(id_);

flag_new_OE = 0;
    file_name = dir([pat_to_raw filesep '*_CH*.continuous']);
    if isempty(file_name)
        file_name = dir([pat_to_raw filesep '*_*.continuous']);
        flag_new_OE = 1;
    end
    init_num = file_name(1).name(1:3);
    
    
    fprintf('Loading Continuous Files: \n');
    if flag_new_OE == 1
        [ch_data_, ~, ~] = ...
            load_open_ephys_data_faster([pat_to_raw filesep init_num...
            '_' num2str(id_(1)) '.continuous']);
    else
        [ch_data_, ~, ~] = ...
            load_open_ephys_data_faster([pat_to_raw filesep init_num...
            '_CH' num2str(id_(1)) '.continuous']);
    end
    
    ch_data = zeros(nofCh,size(ch_data_,1),'int16');
    
    for i = 1:nofCh
        fprintf('-')
        if i == 1
            ch_data(i,:) = int16(ch_data_);
            continue;
        end
        try
            if flag_new_OE
                [ch_data_, ~, ~] =  load_open_ephys_data_faster([pat_to_raw filesep init_num '_' num2str(id_(i)) '.continuous']);
            else
                [ch_data_, ~, ~] =  load_open_ephys_data_faster([pat_to_raw filesep init_num '_CH' num2str(id_(i)) '.continuous']);
            end
        catch
            ch_data_ = zeros(1,size(ch_data_,1),'int16');
        end
        ch_data(i,:) = int16(ch_data_);
        
    end
    fprintf('-> \n')
end
%% Create metadata file
function createMetadataTextfile(filename, varargin)
    % Writes variables to a text file with one variable per row
    % Usage: writeVariablesToTxt(filename, 'var1', values1, 'var2', values2, ...)
    %
    % Inputs:
    %   filename  - Output text file name (e.g., '2024-03-11_camb_offsets.txt')
    %   varargin  - Alternating variable names and their values
    % Example:
    %   writeVariablesToTxt('output.txt', ...
    %       'smp_imap0', [0, 49224704], ...
    %       'sec_imap0', [1640.823467, 2516.514133], ...
    %       'srate_imap0', [30000, 30000], ...
    %       'fileTimeSecs', 2516.514133);

    % Validate inputs
    if mod(nargin-1, 2) ~= 0
        error('Inputs must be: filename followed by variable name-value pairs');
    end
    
    % Open file for writing
    fid = fopen(filename, 'wt');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    % Process each variable
    for i = 1:2:length(varargin)
        varName = varargin{i};
        values = varargin{i+1};
        
        % Write variable name with colon
        fprintf(fid, '%s:\t', varName);
        
        % Format values
        if isnumeric(values)
            % Auto-detect integer/float formatting
            fmt = arrayfun(@(x) num2str(x, chooseFormat(x)), values, 'UniformOutput', false);
        else
            % Handle string/cell array inputs
            fmt = cellstr(string(values));
        end
        
        % Write values separated by tabs
        fprintf(fid, '%s\t', fmt{:});
        
        % Newline with Windows-style line ending
        fprintf(fid, '\r\n');
    end
    
    fclose(fid);
    fprintf('File "%s" created successfully.\n', filename);
end

% Helper function to choose numeric format
function fmt = chooseFormat(x)
    if abs(x - round(x)) < eps(x)
        fmt = '%.0f';  % Integer format
    else
        fmt = '%.6f'; % Float format (6 decimal places)
    end
end