% Move files
clear all
close all
clc
%%
parent_folder = 'H:\My Documents\Multicellular automaton\data\two_signals\parameter_gradient';
subfolder = 'vertical_step_function_Ay_0p0';
%subfolder = 'vertical_step_function_Ay_0p1';
%subfolder = strrep(sprintf('vertical_step_function_Ay_%.1f', Ay), '.', 'p');
folder = fullfile(parent_folder, subfolder);
listing = dir(folder);
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
%% destination files
subfolder_dest = 'negative_control';
folder_dest = fullfile(parent_folder, subfolder_dest);
listing = dir(folder_dest);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
names_dest = {};
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names_dest{count} = name;
    end
end
%% match name with name in destination

% get original name
name_orig = names{1};
name_exists = check_fname(name_orig, names_dest);
pattern = 'Parameter_gradient_K_2_1_square_wave_t_out_(\d+)_period_(\d+|Inf)-v(\d+)';
[tokens, match] = regexp(name_orig, pattern, 'tokens', 'match');
v = str2double(tokens{1}{3});
%%
while name_exists
    v = v+1;
    disp(v);
    new_name = sprintf('Parameter_gradient_K_2_1_square_wave_t_out_%d_period_%d-v%d', ...
        str2double(tokens{1}{1}), str2double(tokens{1}{2}), v);
    name_exists = check_fname(new_name, names_dest);
end
%% Move files
fname_orig = fullfile(folder, subfolder, strcat(name_orig, '.mat'));
fname_dest = fullfile(folder, subfolder_dest, new_name);
copyfile(fname_orig, fname_dest);

%%
function name_exists = check_fname(name_orig, names_dest)
    name_exists = 0;
    for i=1:numel(names_dest)
        if strcmp(names_dest{i}, name_orig)
            name_exists = 1;
            return
        end
    end
end