function values = JDM_read_meta_var(meta_file, var_names)
   % Read the metadata file line by line
    fid = fopen(meta_file, 'rt');
    if fid == -1
        error('Cannot open the file.');
    end

    % Initialize output as a cell array
    values = cell(size(var_names));
    
    % Read the entire file into a cell array
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = lines{1};

    % Loop through each requested variable
    for i = 1:length(var_names)
        var_name = var_names{i};
        
        % Search for the variable in the metadata lines
        match = startsWith(lines, [var_name '=']);
        
        if any(match)
            line = lines{match};
            [~, val] = strtok(line, '=');
            val = strtrim(val(2:end)); % Remove '='
            
            % Extract the first number if commas are present
            if contains(val, ',')
                values{i} = str2double(extractBefore(val, ','));
            else
                values{i} = val;
            end
        else
            warning('Variable "%s" not found in the metadata file.', var_name);
            values{i} = [];
        end
    end
end

