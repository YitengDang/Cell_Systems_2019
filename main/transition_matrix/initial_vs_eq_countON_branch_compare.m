clear variables
close all
warning off

% Plot the saved pin_pout map calculated with the branch algorith and
% compare with the result of the exact simulation

% Save Figure?
save_all = 1;

% Get analytical data
data_path = 'C:\Users\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';
[fname, path, ~] = uigetfile(fullfile(data_path,'pin_pout_branch'));
load(fullfile(path,fname));
p = (0:N)/N;
[~, name, ~] = fileparts(fname);

% Plot analytical map
h1 = figure(1);
imfig = imagesc(p,p,out_map);
set(imfig, 'AlphaData', out_map>1e-8)
set(gca, 'ydir', 'normal')
set(gca,'FontSize', 20)
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{ini}', 'FontSize', 24)
ylabel('p_{eq}', 'FontSize', 24)
% Save the pdf
if save_all > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout_branch', strcat(name,'_ana_map'));
    save_figure_pdf(h1, 10, 6, out_file);
end

% Plot the average number of steps and compare with exact simulation
h2 = figure(2);
plot(p, t_av, 'k', 'Linewidth', 2)
hold on
% Ask for comparison simulation data
[fname, path, ~] = uigetfile(fullfile(data_path,'pin_pout'));
load(fullfile(path,fname));
p = (0:N)./N;
%prob = transpose(count./repmat(sum(count,2),1,N+1));
axes(h2.CurrentAxes)
plot(p, t_av, 'ro')
hold off
legend({'Branch'; 'Automata'},'Location','northwest')
set(gca,'FontSize', 20)
xlabel('p_{ini}', 'FontSize', 24)
ylabel('n_{steps}', 'FontSize', 24)
% Save the pdf
if save_all > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout_branch', strcat(name,'_tav'));
    save_figure_pdf(h2, 10, 6, out_file);
end

h3 = figure(3);
im_fig = imagesc(p,p,prob);
set(gca,'Ydir','normal','FontSize', 20)
set(im_fig, 'AlphaData', count' > 0);
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{ini}', 'FontSize', 24)
ylabel('p_{eq}', 'FontSize', 24)

if save_all > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout_analytical', strcat(name,'_map'));
    save_figure_pdf(h3, 10, 6, out_file);
end