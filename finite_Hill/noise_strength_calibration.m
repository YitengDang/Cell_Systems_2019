% Generate trajectories from time_evolution_finiteHill_save_noise_v2.m
% first
close all
clear variables
warning off
%% Parameters
% Lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 5;
Rcell = 0.2*a0;

% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient
prec = 8;
%noise = 0;
noiselist = 10.^[-1:0.1:1];

% simulation 
%a0list = [5 5.4 6];
initialID = 'uniform';
tmax = 1000;
nruns = 20;

%% Load data
path = 'H:\My Documents\Multicellular automaton\temp'; 

dY_mean_all = zeros(numel(noiselist), nruns); % double averaging: (1) <dY_i> over cells, (2) <<dY_mean>> over times
dY_std_all = zeros(numel(noiselist), nruns); 

for i=1:numel(noiselist)
    noise = noiselist(i);
    for j=1:nruns
        fprintf('noise = %.2f, run = %d \n', noise, j);
        fname_str = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_noise%.2f_prec_%d_tmax%d_%s-v%d',...
            N, a0, K, Con, hill, noise, prec, tmax, initialID, j), '.','p');
        fname = fullfile(path, fname_str);
        load(fname, 'dY_mean', 'dY_std');
        dY_mean_all(i,j) = mean(dY_mean); % average over times
        dY_std_all(i,j) = mean(dY_std); % average over times
    end
end    
%% Plot results
h1=figure(1);
semilogx(noiselist, mean(dY_mean_all, 2), 'bo--', 'LineWidth', 2);
hold on
semilogx(noiselist, mean(dY_std_all, 2), 'ro--', 'LineWidth', 2);
plot([noiselist(1) noiselist(end)], [0 0], 'b-');
plot(noiselist, noiselist, 'r-');
xlabel('Noise strength');
ylabel('dY moments');
set(gca,'FontSize', 24);
%legend({'mean', 'std'});

qsave = 1;
if qsave
    dir = 'H:\My Documents\Multicellular automaton\temp';
    fname_str = strrep(sprintf('dY_moments_vs_noise_N%d_a0_%.2f_K%d_Con%d_hill%.2f_prec_%d_tmax%d_%s_blue_mean_red_std',...
            N, a0, K, Con, hill, prec, tmax, initialID), '.','p');
    out_file = fullfile(dir, fname_str); %filename
    save_figure_pdf(h1, 10, 8, out_file);
    save_figure_eps(h1, 10, 8, out_file);
end
