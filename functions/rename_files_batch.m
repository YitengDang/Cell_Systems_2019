clear all

pattern = 'propagation';
new_pattern = 'formation';

current_path = pwd;
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\trav_wave_reliability';
names = dir(folder);
cd(folder); % change to folder with files

for i=3:numel(names)
    name = names(i).name;
    if ~isempty(regexp(name, pattern, 'match', 'once'))
        disp(name)
        new_name = strrep(name, pattern, new_pattern);
        [status, msg] = copyfile(  name, new_name );
        if status
            disp(new_name);
        end
    end
end

cd(current_path); % change back to original folder
