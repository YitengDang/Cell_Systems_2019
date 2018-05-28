% Compare the analytical simulation (Monte Carlo) of the map p_in p_eq with
% the exact simulation. When running this you will be asked to choose one
% file, which is the saved mat file from the analytical simulation. The
% data from the Monte Carlo simulation is gotten automatocally.
clear variables
close all
warning off

% Get the exact simulation data
data_path = 'C:\Users\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';
[fname, path, ~] = uigetfile(fullfile(data_path,'pin_pout'));
load(fullfile(path,fname))
[~,name,~] = fileparts(fname);
% Get the Monte Carlo data
aux = load(fullfile(path,strcat(name,'_montecarlo.mat')), 't_av', 'cntmap','I_av');
t_av_analytic = aux.t_av;
prob_analytic = aux.cntmap;
I_av_analytic = aux.I_av;

% Save Figure?
save_fig = 1;

[path, name, ~] = fileparts(fname);
% Calculate the probability map
p = (0:N)./N;
prob = transpose(count./repmat(sum(count,2),1,N+1));

% Plot the map
h1 = figure(1);
im_fig = imagesc(p,p,prob);
% Plot lines that will be further analysis
vline_p = [];
hold on
for a = vline_p
plot([a a],[0 1], '--r', 'Linewidth', 1.5)
end
hold off
% Set the font, labels and the invisibility of parts that show no count
set(gca,'Ydir','normal','FontSize', 20)
set(im_fig, 'AlphaData', count' > 0);
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{ini}', 'FontSize', 24)
ylabel('p_{eq}', 'FontSize', 24)

% Organize and save pdf
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout', strcat(name,'_map_sim'));
    save_figure_pdf(h1, 10, 6, out_file);
end

% Plot and compare the average number of steps and compare both
h2 = figure(2);
plot(p, t_av, 'ko')
hold on
plot(p, t_av_analytic, 'r-')
hold off
legend({'Automata','Analytic'})
set(gca,'FontSize', 20)
xlabel('p_{ini}', 'FontSize', 24)
ylabel('n_{steps}', 'FontSize', 24)
% Organize and save pdf
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout', strcat(name,'_time'));
    save_figure_pdf(h2, 10, 6, out_file);
end

% Plot the analytical (Monte Carlo) map 
h3 = figure(3);
im_fig = imagesc(p,p,prob_analytic');
% Plot lines
vline_p = [];
hold on
for a = vline_p
plot([a a],[0 1], '--r', 'Linewidth', 1.5)
end
hold off
% Set the font, invibility and labels
set(gca,'Ydir','normal','FontSize', 20)
set(im_fig, 'AlphaData', prob_analytic' > 0);
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{ini}', 'FontSize', 24)
ylabel('p_{eq}', 'FontSize', 24)

% Organize and save
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout', strcat(name,'_map_analytic'));
    save_figure_pdf(h3, 10, 6, out_file);
end

if exist('I_av', 'var')
    h4 = figure(4);
    % Plot the spatial order parameter for simulated data
    plot(p, I_av, 'LineWidth', 1.5)
    set(gca,'Ydir','normal','FontSize', 20)
    xlabel('p_{ini}', 'FontSize', 24)
    ylabel('Spatial order I', 'FontSize', 24)
    box on
    legend({'Average';'Min';'Max'})
    % Organize and save
    if save_fig > 0
        out_file = fullfile(pwd, 'figures', 'pin_pout', strcat(name,'_I'));
        save_figure_pdf(h4, 10, 6, out_file);
    end
    yl = ylim;
    h5 = figure(5);
    % Plot the spatial order parameter for analytical data
    plot(p, I_av_analytic, 'LineWidth', 1.5)
    set(gca,'Ydir','normal','FontSize', 20)
    xlabel('p_{ini}', 'FontSize', 24)
    ylabel('Spatial order I', 'FontSize', 24)
    box on
    legend({'Average';'Min';'Max'})
    ylim(yl)
    % Organize and save
    if save_fig > 0
        out_file = fullfile(pwd, 'figures', 'pin_pout', strcat(name,'_Ianalytic'));
        save_figure_pdf(h5, 10, 6, out_file);
    end
end