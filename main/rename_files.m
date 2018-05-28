%% Renames stored files with given pattern in folder according to given rule
% Here applied to saved paths for the finite Hill case

% Lattice parameters
gridsize = 11;
N = gridsize^2;
initialID = 'uniform';

% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient

%a0 = 5.0;
all_a0 = 5.0:0.1:6.0;

for j=1:11
    a0 = all_a0(j);
    
    
    % search for data
    path = fullfile(pwd, 'data', 'dynamics', 'finiteHill'); 
    straux = '(\d+)';
    straux2 = '(\w+)';
    fpattern = strrep(sprintf('N%d_a0_%d_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_%s-v%s',...
                N, 10*a0, K, Con, hill, straux, straux2, initialID, straux), '.', 'p');

    % Get all file names in the directory
    listing = dir(path);
    num_files = numel(listing)-2; %first two entries are not useful
    count = 0;

    for i = 1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            count = count + 1;
            names{count} = name;
        end
    end

    for i = 1:numel(names)
        % first get the filename given by the student
        [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
        if numel(tokens) > 0
            disp(names{i}) % displays the file name
            old_str = sprintf('a0_%d', a0*10);
            new_str = strrep(sprintf('a0_%.2f', a0), '.', 'p');
            new_name = strrep(names{i}, old_str, new_str);
            disp(new_name);
            % load the data
            load(fullfile(path,strcat(names{i},'.mat')),...
                'cells_hist', 'I', 'mom', 'xmean', 'xstd');
            save(fullfile(path, 'temp', new_name));
        end
    end

end
