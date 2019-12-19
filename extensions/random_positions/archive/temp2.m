clear all
close all
clc
%% folder to save simulations in
network = 19;
%parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\randomized lattice';
parent_folder = 'W:\staff-bulk\tnw\BN\HY\Shared\Yiteng\two_signals\randomized lattice';

subfolder = sprintf('TW_propagation_network_%d', network);
%subfolder = sprintf('TW_formation_network_%d', network);
save_folder = fullfile(parent_folder, subfolder);

files = dir(save_folder);
files = files(3:end);
for ii=1:numel(files)
    old_name = files(ii).name;
    new_name = strrep(old_name, 'TW', 'TW_params_1');
    disp(new_name);
    copyfile(fullfile(save_folder, old_name),...
        fullfile(save_folder, new_name));
end