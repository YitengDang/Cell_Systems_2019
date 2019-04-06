clear all
close all
clc
%% 
folder = 'H:\My Documents\Multicellular automaton\figures\trav_wave_stability\spider_plots_analytical';
content = dir(folder);
for i=3:size(content, 1)
    name = content(i).name;
    
    new_name = strrep(name, 'simulations', 'analytical');
    
    % copies the file name into a new file, does not delete old file
    copyfile(fullfile(folder, name), fullfile(folder, new_name));
end
