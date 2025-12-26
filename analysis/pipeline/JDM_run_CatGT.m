function JDM_run_CatGT(path_data_monkey_sorted,sess)

answer_ = questdlg('What CatGT option?', ...
 'Select option', ...
 'Single', 'Supercat', 'Single');  % Default option is 'Single'

% find the path in the computer where the CatGT folder contains .bat file
Path_Str       = regexp(path,pathsep,'split');
idx_gt         = contains(Path_Str,[filesep,'CatGT-win']);
all_GT_paths   = Path_Str(idx_gt);
batFilePath    = all_GT_paths{1};

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


rec_num       = numel(rec_name);
gValues       = cell(1,rec_num);
dirPaths      = cell(1,rec_num);
supercatPaths = cell(1,rec_num);
runValues     = cell(1,rec_num);


for counter_rec = 1:rec_num
    pat_to_raw = fullfile(rec_path{counter_rec}, rec_name{counter_rec}, 'raw_data');
    if ~isfolder(pat_to_raw)
        error("There is no raw_data folder in %s", fullfile(rec_path{counter_rec}, rec_name{counter_rec}));
    end
    
    % Check for .ap.bin and .ap.meta files
    ap_files = dir2(fullfile(pat_to_raw, '*.ap.bin*'));
    meta_files = dir2(fullfile(pat_to_raw, '*.ap.meta'));
    
    if isempty(ap_files) || isempty(meta_files)
        error('Binary data or Metadata missing in %s', pat_to_raw);
    end
    
    % Extract newRun and new_g using regex
    filename = ap_files(1).name;
    match = regexp(filename, '^(.*?)_g(\d+)', 'tokens', 'once');
    if isempty(match) || numel(match) < 2
        error('Filename %s does not match expected format (run_gX.ap.bin)', filename);
    end
    newRun = match{1};
    new_g = match{2};
    
    % Assign extracted values
    gValues{counter_rec} = new_g;
    dirPaths{counter_rec} = pat_to_raw;
    runValues{counter_rec} = newRun;
    supercatPaths{counter_rec} = sprintf('{%s,%s_g%s}', pat_to_raw, newRun, new_g);
end

switch answer_
    case 'Single'
    cd (batFilePath)
    tic;
        for counter_rec = 1:rec_num
           current_rec = rec_name{counter_rec};
           fprintf('   %d/%d: %s\n', counter_rec, rec_num, current_rec);
           
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
         ap_file = dir2(fullfile(rec_path{counter_rec},current_rec,'raw_data','*.ap.bin*'));
    
         if isempty(ap_file)
             disp("Data .bin does not exist, check the raw data")
             continue
         end
          cmd_ = ['CatGT -dir=' dirPaths{counter_rec} ' -run=' runValues{counter_rec} ' '...
             '-g=' gValues{counter_rec} ' -t=0 -no_run_fld -t_miss_ok -ap -prb=0 '...
             '-xd=2,0,384,6,0 ' '-apfilter=butter,12,300,9000 ' '-gblcar '...
             '-xid=2,0,384,6,0'];
    
         [status, result] = system(cmd_);
         if status ~= 0
            fprintf('Error executing the .bat file.\n');
            fprintf('System result: %s\n', result);
         else
            ap_file = dir2(fullfile(rec_path{counter_rec},current_rec,'raw_data', '*_tcat*.ap.bin'));
            md_file = dir2(fullfile(rec_path{counter_rec},current_rec,'raw_data', '*_tcat*.ap.meta'));
            movefile([ap_file.folder filesep ap_file.name], path_to_sorted);
            movefile([md_file.folder filesep md_file.name], path_to_sorted);
            fprintf('Successfully executed the .bat file.\n');
         end

        end
       fprintf('Elapsed Time for the creating CatGT and TTL files is %d\n',toc)
       
    case 'Supercat'
     tic;
     cd (batFilePath)
    for counter_bat = 1:rec_num

         fprintf(['   ' num2str(counter_bat) '/' num2str(rec_num) ': ' rec_name{counter_bat} '\n'])

        ap_file_cat   = dir2(fullfile(pat_to_raw,'*_tcat*.ap.bin'));

        if ~isempty(ap_file_cat)
             warning('The CatGT file exist, overwriting on the old data'); 
        end
    
         cmd_ = ['CatGT -dir=' dirPaths{counter_bat} ' -run=' runValues{counter_bat} ' '...
             '-g=' gValues{counter_bat} ' -t=0 -no_run_fld -t_miss_ok -ap -prb=0 '...
             '-xd=2,0,384,6,0 ' ...
             '-xid=2,0,384,6,0'];
    
         [status, result] = system(cmd_);
         if status ~= 0
            fprintf('Error executing the .bat file.\n');
            fprintf('System result: %s\n', result);
         else
            fprintf('Successfully executed the .bat file.\n');
         end
           
     end
  
     % run supercat on the .cat data
 
    fprintf('Elapsed Time for the creating CatGT and TTL files is %d\n',toc)

     tic;
         cat_recs_path = fullfile(rec_path{counter_rec},'cat_recs\analyzed_data\sorted_data');
          if ~isfolder(cat_recs_path)
              mkdir(cat_recs_path);
          end
    
          for counter_rec = 1:rec_num
                pat_to_raw    = fullfile(rec_path{counter_rec},rec_name{counter_rec},'raw_data');
                
                % Check the data contain binary NPX data and Metadata
                % get the name of .bin and . meta
                ap_file_chk   = dir2(fullfile(pat_to_raw,'*_tcat*.ap.bin'));
                meta_file_chk = dir2(fullfile(pat_to_raw,'*_tcat*.ap.meta'));
                
                % Display error if there is no data
                if isempty(ap_file_chk.name) || isempty(meta_file_chk.name)
                    error('There is no Binary data or Meta data file!') 
                end
          end
    
          temp_cats = strjoin(supercatPaths, ' ');
    
          % Split the string by spaces
          parts = strsplit(temp_cats, ' ');
            
          % Remove any empty strings
          parts = parts(~cellfun('isempty', parts));
            
          % Concatenate all the parts into a single string
          concatenatedResult = [parts{:}];
            
         cmd_cat = ['CatGT -ap -prb=0 -no_run_fld ' '-xd=2,0,384,6,0 '...
             '-apfilter=butter,12,300,9000 ' '-gblcar '...
             '-xid=2,0,384,6,0 ' '-no_catgt_fld ' '-zerofillmax=500 '...
             '-dest=' cat_recs_path ' ' '-supercat=' concatenatedResult ' '];
    
         [status, result] = system(cmd_cat);
    
         if status ~= 0
            fprintf('Error executing supercat file.\n');
            fprintf('System result: %s\n', result);
         else
                for counter_rec = 1:rec_num
                    pat_to_raw    = fullfile(rec_path{counter_rec},rec_name{counter_rec},'raw_data');
                    % Check the data contain binary NPX data and Metadata and
                    % remove them
                    % get the name of .bin and . meta
                    ap_file_chk   = dir2(fullfile(pat_to_raw,'*_tcat*.ap.bin'));
                    delete(fullfile(pat_to_raw,ap_file_chk.name))
                end
            
            fprintf('Successfully executed supercat file.\n');
         end
         cat_folder = dir2(cat_recs_path);
         fileList = dir2(fullfile(cat_folder.folder,cat_folder.name));
    
             % Loop through each file and move it to the destination folder
        for i = 1:length(fileList)
            % Skip '.' and '..' directories
            if ~fileList(i).isdir
                % Get the full path of the source file
                sourceFile = fullfile(fileList(i).folder, fileList(i).name);
                % Move the file to the destination folder
                movefile(sourceFile, cat_recs_path);
            end
        end
        rmdir(fullfile(cat_folder.folder,cat_folder.name))
     fprintf('Elapsed Time for the concatenation files is %d\n',toc)
    otherwise
    disp('No selection made.');
        
end

end
 
