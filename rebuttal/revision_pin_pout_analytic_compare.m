% Compare the analytical simulation (Monte Carlo) of the map p_in p_eq with
% the exact simulation. When running this you will be asked to choose one
% file, which is the saved mat file from the analytical simulation. The
% data from the Monte Carlo simulation is gotten automatocally.
clear variables
close all
warning off

% Parameters
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
K = 6;
Con = 15;
p = (0:N)./N;
%delta = 0.004;
%noise = 2;

a0list = [1.5 1.5 1.5 0.5 0.5 0.5];
Klist = [6 20 12 15 14 14];
Conlist = [15 15 15 8 12 5];

for idx=1:6
    a0 = a0list(idx);
    K = Klist(idx);
    Con = Conlist(idx);
    
% Get the exact simulation data
fname_str_exact = sprintf('exact_sim_pin_pout_N%d_Con_%d_K_%d_a0_%d', N, Con, K, a0*10);
fname = fullfile(pwd, 'rebuttal', 'Langevin_pin_pout', 'data', fname_str_exact);
load(fname, 'prob', 't_av');

% Load Langevin trajectory
%fname_str = strrep(sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d_delta_%.3f_noise_%d', ...
%            N, Con, K, 10*a0, delta, noise), '.','p');
%fname_str = strrep(sprintf('no_noise_params_pin_pout_N%d_Con_%d_K_%d_a0_%d', ...
%            N, Con, K, 10*a0), '.','p');
fname_str = strrep(sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d_delta_sigma_calc', ...
            N, Con, K, 10*a0), '.','p');
fname = fullfile(pwd, 'rebuttal', 'Langevin_pin_pout', 'data', strcat('Langevin_', fname_str, '.mat') );
aux = load(fname);
t_av_Langevin = aux.t_av;
prob_Langevin = aux.cntmap;
%I_av_analytic = aux.I_av;
%}

temp=sum(prob_Langevin, 2);
disp(min(temp));
disp(max(temp));
% Save Figures?
save_fig = 1;
%% Plot exact simulation map
%prob = transpose(count./repmat(sum(count,2),1,N+1));

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
set(im_fig, 'AlphaData', prob > 0);
c = colorbar;
c.Label.String = 'Probability';
caxis([0 1]);
xlabel('p_{ini}', 'FontSize', 24)
ylabel('p_{eq}', 'FontSize', 24)

% Organize and save pdf
if save_fig > 0
    out_file = fullfile(pwd, 'rebuttal', 'Langevin_pin_pout',...
        strcat(fname_str,'_map_sim'));
    %save_figure_pdf(h1, 10, 6, out_file);
    save_figure_eps(h1, 10, 6, out_file);
end
%%
% Plot and compare the average number of steps and compare both
h2 = figure(2);
plot(p, t_av, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
hold on
plot(p, t_av_Langevin, 'r-', 'LineWidth', 2)
hold off
%legend({'Automata','Langevin'})
set(gca,'FontSize', 20)
xlabel('p_{ini}', 'FontSize', 24)
ylabel('n_{steps}', 'FontSize', 24)
%ylim([0 1]);
% Organize and save pdf
if save_fig > 0
    out_file = fullfile(pwd, 'rebuttal', 'Langevin_pin_pout',...
        strcat(fname_str,'_time'));
    %save_figure_pdf(h2, 10, 6, out_file);
    save_figure_eps(h2, 10, 6, out_file);
end
%%
% Plot the analytical (Monte Carlo) map 
h3 = figure(3);
im_fig = imagesc(p,p,prob_Langevin');
% Plot lines
vline_p = [];
hold on
for a = vline_p
plot([a a],[0 1], '--r', 'Linewidth', 1.5)
end
hold off
% Set the font, invibility and labels
set(gca,'Ydir','normal','FontSize', 20)
set(im_fig, 'AlphaData', prob_Langevin' > 0);
c = colorbar;
caxis([0 1]);
c.Label.String = 'Probability';
xlabel('p_{ini}', 'FontSize', 24)
ylabel('p_{eq}', 'FontSize', 24)

% Organize and save
if save_fig > 0
    out_file = fullfile(pwd, 'rebuttal', 'Langevin_pin_pout',...
        strcat(fname_str,'_map_analytic'));
    %save_figure_pdf(h3, 10, 6, out_file);
    save_figure_eps(h3, 10, 6, out_file);
end
end