function JDM_concate_RHS(path_data_monkey_sorted,sess,probe_flag,cat_flag)


cambridge_M1_rhs = [6, 5, 4, 31, 30, 8, 7, 32, 26, 25, 64, 3, 2, ...
    29, 28,	27,	40,	39,	59,	60,	61,	62,	63,	1, 33, 34, 35, 36,...
    37,	38,	58,	57,	41,	42,	54,	53,	52,	51,	50,	16,	48,	47,	56,...
    49,	46,	9,	11,	14,	45,	10,	12,	15,	44,	17,	13,	20,	43,	23,	18,...
    21,	55,	24,	19,	22];

% 1: linear M2
id_linear = [28,26,24,19,21,15,32,30,18,29,31,22,20,23,25,27,1,3,6,8,10,12,14,16,17,13,11,...
            9,7,5,4,2,64,62,59,57,55,53,51,47,50,52,54,56,58,60,61,63,37,39,41,46,44,33,35,48,...
            36,34,49,43,45,42,40,38];


if probe_flag == 6
    chan_map = cambridge_M1_rhs;
elseif probe_flag == 7
    chan_map = id_linear;
else
    error('Channel map is not corrects')
end

if cat_flag

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
    
    cat_recs_path = fullfile(rec_path{1},'cat_recs\analyzed_data\sorted_data');
    if ~isfolder(cat_recs_path)
      mkdir(cat_recs_path);
    end
    
    % Initialize arrays for data and timestamps
    cat_data = []; 
    smp_imap  = nan(rec_num,1); 
    sec_imap  = nan(rec_num,1);
    
    for counter_rec = 1:rec_num
    
           fprintf(['      ', num2str(counter_rec), ' / ' num2str(rec_num)...
                ' Analyzing recording ', rec_name{counter_rec} '\n']);
    
        pat_to_raw = fullfile(rec_path{counter_rec}, rec_name{counter_rec}, 'raw_data', filesep);
        if ~isfolder(pat_to_raw)
            fprintf("Warning: There is no raw_data folder in %s\n", fullfile(rec_path{counter_rec}, rec_name{counter_rec}));
            continue;
        end
        [amplifier_data,record_time,num_amplifier_samples] = JDM_read_Intan_RHS2000_file_v1(pat_to_raw);
        
        cat_data = [cat_data, amplifier_data(chan_map,:)]; 
        smp_imap(counter_rec)  = num_amplifier_samples;
        sec_imap(counter_rec)  = record_time;
        clear amplifier_data
    end
    
    fprintf(' \n -> Concatenation of the data completed. Converting to Binary \n');
    
    data_name   = [sess '_tcat.camb0.ap'];
    offset_name = [sess '_camb_offsets'];
    
    
    % Save as binary file
    fid = fopen([cat_recs_path filesep data_name '.bin'], 'w');
    fwrite(fid, cat_data, 'int16');
    fclose(fid);
    
    clear cat_data
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
    
    
    fprintf(1, 'Generating a Metadata for the recording.  \n');
    fprintf(1, 'Done! \n');
else
    % no concatenate 
    pattern_sess = '(.*\\\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})'; 
    sess_name = regexp(path_data_monkey_sorted, pattern_sess, 'match', 'once');
    
    pattern_rec = '\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}'; 
    rec_name = regexp(path_data_monkey_sorted, pattern_rec, 'match', 'once');

   [amplifier_data,~,~] = JDM_read_Intan_RHS2000_file_v1(path_data_monkey_sorted);

   amplifier_data = amplifier_data(chan_map,:);

    if sum(amplifier_data(:)) == 0
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
    fprintf(1,'Data conversion to binary is done. \n')
end
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
