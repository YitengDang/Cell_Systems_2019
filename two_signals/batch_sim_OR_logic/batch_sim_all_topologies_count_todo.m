function [sim_todo, filecount] = batch_sim_all_topologies_count_todo(...
    max_trials, folder, pattern)
    
    % Count how many simulations have already been done
    if exist(folder, 'dir') ~= 7
        warning('Folder does not exist! ');
    end

    filecount = 0;
    %pattern = 'all_topologies_simulate-v(\d+)';
    
    % original method
    %
    listing = dir(folder);
    num_files = numel(listing)-2;
    names = {};
    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~, name, ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            match = regexp(name, pattern, 'match');
            %disp(match);
            if ~isempty(match)
                filecount = filecount + 1;
                names{end+1} = name;
            end
        end
    end
    %}
    
    % more accurate method (simulate file loading), slower
    %{
    for idx2=1:max_trials
        % load file
        fname_str = sprintf('all_topologies_simulate-v%d.mat', idx2);
        fname = fullfile(folder, fname_str);
        if exist(fname, 'file')~=2
            fname_str = sprintf('all_topologies_simulate-v%d_tmax_reached.mat', idx2);
            fname = fullfile(folder, fname_str);
        end
        % if new file name also doens't exist, continue
        if exist(fname, 'file')~=2
            continue;
        end
        filecount = filecount + 1;
    end
    %}
    
    sim_todo = max_trials-filecount;
    fprintf('Sim to do: %d \n', sim_todo);
    
end