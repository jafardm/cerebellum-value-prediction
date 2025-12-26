%% Function to extract metadata from the meta file
function meta_data = JDM_extract_metadata(filename)
fid = fopen(filename, 'r');
meta_data = struct();

while ~feof(fid)
    line = fgetl(fid);
    if contains(line, 'fileTimeSecs')
        meta_data.fileTimeSecs = str2double(extractAfter(line, '='));
    elseif contains(line, 'firstSample')
        meta_data.firstSample = str2double(extractAfter(line, '='));
    elseif contains(line, 'imSampRate')
        meta_data.sampRate = str2double(extractAfter(line, '='));
    end
end
fclose(fid);
end