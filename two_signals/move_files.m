% Warning: overwrites existing files in target folder with same name
% Warning: check if target folder exists (otherwise, files are lost)
clear all

pattern = 'v2[0-9].mat';
file_ext = '.mat'; % extension of the files to move

%path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution';
path = 'H:\My Documents\Multicellular automaton\data\two_signals\time_evolution\III_chaotic_scan_p_ini';
target_path = fullfile(path, 'additional');

listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, file_ext)
        count = count + 1;
        names{count} = name;
    end
end

%%
for i=1:count
    fname = strcat(names{i}, file_ext);
    if ~isempty(regexp(fname, pattern, 'once'))
        disp(fname);
        [~, msg] = movefile(fullfile(path, fname), target_path);
        %disp(status);
        disp(msg);
    end
end