function [sim_todo, filecount] = all_topologies_simulate_count_todo(...
    max_trials, folder, pattern)
    
    % Count how many simulations have already been done
    if exist(folder, 'dir') ~= 7
        warning('Folder does not exist! ');
    end

    filecount = 0;
    %pattern = 'all_topologies_simulate-v(\d+)';

    listing = dir(folder);
    num_files = numel(listing)-2;
    names = {};
    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            match = regexp(name, pattern, 'match');
            %disp(match);
            if ~isempty(match)
                filecount = filecount + 1;
                names{end+1} = name;
            end
        end
    end
    
    sim_todo = max_trials-filecount;
    fprintf('Sim to do: %d \n', max_trials-filecount);
    
end