clear variables
close all
warning off

% Plot the transition matrix for a saved branch result
% Get analytical data
data_path = pwd;
[fname, path, ~] = uigetfile(fullfile(data_path,'data','pin_pout_branch'));
load(fullfile(path,fname));
p = (0:N)/N;
[path, name, ~] = fileparts(fname);

% Save Figure?
save_all = 1;
for i = 0:N
    pmatrix(i+1, i+1) = pmatrix(i+1, i+1) + pEn(i+1);
end
% Plot analytical data
h1 = figure(1);
imfig = imagesc(p,p,pmatrix);
%set(imfig, 'AlphaData', out_map>0)
set(gca, 'ydir', 'normal')
set(gca,'FontSize', 20)
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{t}', 'FontSize', 24)
ylabel('p_{t+1}', 'FontSize', 24)

if save_all > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout_analytical', strcat(name,'_pmatrix'));
    save_figure_pdf(h1, 9, 6, out_file);
end
