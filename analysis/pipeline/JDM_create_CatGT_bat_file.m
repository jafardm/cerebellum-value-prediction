function JDM_create_CatGT_bat_file(path_data_monkey_sorted,sess)
%JDM_create_CatGT_bat_file Generates a .bat file with specified variables.
%
% Usage:
% create_bat_file(dir, run, g, t, prb, butter_freqs, gfix, xd, xid)
%
% Example:
% JDM_create_CatGT_bat_file('Y:\Ephys\data_132F\2024-03\2024-03-19)

sess_path = [path_data_monkey_sorted, sess(1:7), filesep, sess, filesep];

% path to the CatGT folder contains .bat file
Path_Str       = regexp(path,pathsep,'split');
idx_gt         = contains(Path_Str,[filesep,'CatGT-win']);
all_GT_paths   = Path_Str(idx_gt);
cat_GT_path    = all_GT_paths{1};

% cat_GT_path = 'C:\Users\Jafar\Desktop\CatGT';

cd (cat_GT_path)

rec_list = dir2(sess_path);
rec_list([rec_list.isdir] == 0) = [];

rec_name = {rec_list.name};
rec_path = {rec_list.folder};

rec_num = numel(rec_name);

for counter_rec = 1:rec_num
    current_rec = rec_name{counter_rec};
    path_to_sorted = fullfile(rec_path{counter_rec} ,current_rec, 'analyzed_data', 'sorted_data', filesep);
     if ~isfolder(path_to_sorted)
         mkdir(path_to_sorted);
     end

    % check if the destination folder has .tact file to delete it
    
    ap_file_chk   = dir2(fullfile(path_to_sorted,'*_tcat*.ap.bin'));
    meta_file_chk = dir2(fullfile(path_to_sorted,'*_tcat*.ap.meta'));

    if ~isempty(ap_file_chk)
        delete(fullfile(ap_file_chk.folder,ap_file_chk.name))
    end

    if ~isempty(meta_file_chk) 
        delete(fullfile(meta_file_chk.folder,meta_file_chk.name))
    end

end


disp('----->> CatGT is running, check the raw data directory to see the progress!! ')

for counter_rec = 1:rec_num

     current_rec = rec_name{counter_rec};
     
     fprintf(['   ' num2str(counter_rec) '/' num2str(rec_num) ': ' current_rec '\n'])
    
     newDir = fullfile(rec_path{counter_rec}, current_rec,'raw_data',filesep);
    
     ap_file = dir2(fullfile(rec_path{counter_rec},current_rec,'raw_data','*.ap.bin*'));

     if isempty(ap_file)
         disp("Data .bin does not exist, check the raw data")
         continue
     end
    
     newRun = ap_file(1).name(1:10);
    
     newG = regexp(ap_file(1).name,['(?<=' newRun '_g*)[0-9]*'],'match');
     newG = newG{1};
    
     path_to_sorted = fullfile(rec_path{counter_rec} ,current_rec, 'analyzed_data', 'sorted_data', filesep);

    

      % Read the content of the batch file
       filePath = fullfile(cat_GT_path, 'runit.bat');
       fileID = fopen(filePath, 'r');
       fileContent = fread(fileID, '*char')';
       fclose(fileID);

%        % Display the original content (optional)
%        disp('Original Content:');
%        disp(fileContent);
%         
        % Extract the current values using regular expressions
        dirPattern = 'set LOCALARGS=-dir=([^\s]+) \^';
        runPattern = '-run=([^\s]+) \^';
        gPattern = '-g=([^\s]+)';
        
        currentDir = regexp(fileContent, dirPattern, 'tokens', 'once');
        currentRun = regexp(fileContent, runPattern, 'tokens', 'once');
        currentG   = regexp(fileContent, gPattern, 'tokens', 'once');
        
        if isempty(currentDir)
            currentDir = {''};
        end
        if isempty(currentRun)
            currentRun = {''};
        end
        if isempty(currentG)
            currentG = {''};
        end
                
        % Escape backslashes in the new directory path
        escapedNewDir = strrep(newDir, '\', '\\');
        
        % Replace the variables in the content
        fileContent = regexprep(fileContent, ['-dir=', regexptranslate('escape', currentDir{1})], ['-dir=', escapedNewDir]);
        fileContent = regexprep(fileContent, ['-run=', regexptranslate('escape', currentRun{1})], ['-run=', newRun]);
        fileContent = regexprep(fileContent, ['-g=', regexptranslate('escape', currentG{1})], ['-g=', newG]);

        
        % Write the modified content back to the file
        fileID = fopen(filePath, 'w');
        fwrite(fileID, fileContent);
        fclose(fileID);

        % Run the batch file using the system command
        status = system(filePath);

        % Check the status
        if status == 0

        % move the .tcat and meta data to the sorted folder

        ap_file = dir2(fullfile(sess_path,current_rec,'raw_data', '*_tcat*.ap.bin'));
        md_file = dir2(fullfile(sess_path,current_rec,'raw_data', '*_tcat*.ap.meta'));
        movefile([ap_file.folder filesep ap_file.name], path_to_sorted);
        movefile([md_file.folder filesep md_file.name], path_to_sorted);
       
            fprintf('--> Batch file executed successfully.\n');
        else
            disp('Error executing batch file.');
        end
        
end
fprintf('----->> Done!')
 
end


