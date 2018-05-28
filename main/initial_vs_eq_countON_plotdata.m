clear variables
close all
warning off

% Plot saved data of a simulated pin_pout map and save as PDF. The user is
% prompted to choose the .mat file with the saved data.

data_path = 'H:\My Documents\Multicellular automaton\data\pin_pout\noise';
[fname, path, ~] = uigetfile(data_path);
noise = 0;
load(fullfile(path,fname));

% Save Figure?
save_fig = 0;

% This makes the code compatible with new data that is saved as Con instead
% of Son
if exist('Con', 'var')
    Son = Con;
end

[path, name, ~] = fileparts(fname);
disp(name);

% Calculate the map
if ~exist('prob', 'var')
    prob = transpose(count./repmat(sum(count,2),1,N+1));
end

p = (0:N)./N;

set(0, 'defaulttextinterpreter', 'latex');
% Plot the map
h1 = figure(1);
im_fig = imagesc(p, p, prob);
% Plot lines used fot futher analysis
%{
vline_p = [0.5-0.25*(0.5*(Son+1)*(1+fN)-K)/fN];
hold on
for a = vline_p
plot([a a],[0 1], '--r', 'Linewidth', 1.5)
end
hold off
%}
set(gca,'Ydir','normal','FontSize', 24)
set(im_fig, 'AlphaData', count' > 0);
c = colorbar;
colormap(flipud(parula));
c.Label.String = 'Probability';
xlabel('$$p_{ini}$$', 'FontSize', 24)
ylabel('$$p_{eq}$$', 'FontSize', 24)
%{
if noise > 0
    title(sprintf('C_{ON}: %.1f, K: %.1f, a0: %.1f, noise: %.2f', ...
        Son, K, a0, noise), 'FontSize', 24)    
else
    title(sprintf('C_{ON}: %.1f, K: %.1f, a0: %.1f', Son, K, a0), ...
        'FontSize', 24)
end
%}
% Organize and save
if save_fig > 1
    out_file = fullfile(pwd, 'figures', 'pin_pout', 'noise', strcat(name,'_map'));
    save_figure_pdf(h1, 10, 8, out_file);
    save_figure_eps(h1, 10, 8, out_file);
end
%%
% Plot the number of average number of steps 
h2 = figure(2);
plot(p, t_av, 'r-o')
set(gca,'FontSize', 24)
xlabel('$$p_{ini}$$', 'FontSize', 24)
ylabel('$$n_{steps}$$', 'FontSize', 24)
%{
if noise > 0
    title(sprintf('C_{ON}: %.1f, K: %.1f, a0: %.1f, noise: %.1f', ...
        Son, K, a0, noise), 'FontSize', 24)    
else
    title(sprintf('C_{ON}: %.1f, K: %.1f, a0: %.1f', Son, K, a0), ...
        'FontSize', 24)
end
%}

% Organize and save
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout', 'noise', strcat(name,'_time'));
    save_figure_pdf(h2, 10, 8, out_file);
    save_figure_eps(h2, 10, 8, out_file);
end
%%
if exist('I_av', 'var')
    h3 = figure(3);
    % Plot the spatial order parameter
    plot(p, I_av, 'LineWidth', 1.5)
    set(gca,'Ydir','normal','FontSize', 24)
    xlabel('$$p_{ini}$$', 'FontSize', 24)
    ylabel('Spatial order $$I$$', 'FontSize', 24)
    box on
    legend({'Average';'Min';'Max'})
    % Organize and save
    if save_fig > 0
        out_file = fullfile(pwd, 'figures', 'pin_pout', 'noise', strcat(name,'_I'));
        save_figure_pdf(h3, 10, 8, out_file);
        save_figure_eps(h3, 10, 8, out_file);
    end
end
%%
h4 = figure(4);
plot(p, sum(count,2), 'x-');