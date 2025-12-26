%% function generate_session_metadata
function generate_session_metadata(path_data_monkey_sorted, sess_list, progress)

if ~strcmp(path_data_monkey_sorted(end), filesep)
    path_data_monkey_sorted = [path_data_monkey_sorted filesep];
end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : num_sess
    progress.Value = counter_sess/num_sess;
    current_sess = sess_list{counter_sess};
    metadata_app = MAF_session_metadata_app(path_data_monkey_sorted, current_sess);
    waitfor(metadata_app);
end

end