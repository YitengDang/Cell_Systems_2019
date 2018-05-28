clear variables
close all
warning off

% Plot the transition matrix for a saved transition matrix result (computed
% from transition_matrix.m)

% Get analytical data
data_path = pwd;
[fname, path, ~] = uigetfile(fullfile(data_path,'data','transition_matrix'));
load(fullfile(path,fname));
p = (0:N)/N;
[path, name, ~] = fileparts(fname);
disp(name)

% Save Figure?
save_all = 1;
% Plot analytical data
h1 = figure(1);
imfig = imagesc(p,p,t_mat');
%set(imfig, 'AlphaData', out_map>0)
set(gca, 'ydir', 'normal')
set(gca,'FontSize', 20)
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{t}', 'FontSize', 24)
ylabel('p_{t+1}', 'FontSize', 24)

if save_all > 0
    out_file = fullfile(pwd, 'figures', 'transition_matrix', strcat(name,'_tmat'));
    save_figure_pdf(h1, 9, 6, out_file);
end
