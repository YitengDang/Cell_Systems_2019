% This script calculates the map of p_in p_eq using exact simulation.
close all
clear all
%warning off
set(0, 'defaulttextinterpreter', 'latex');
%%
% geometric parameters
Lx = 1; % default
n = 11; % nmax = L/R (square packing)
rcell = 0.2;
ini_alg = 0; % 1: Markov MC, 0: random placement
%--------------------
N = n^2; % total number of particles
%N = round(eta*L^2/(pi*R^2));
Ly = sqrt(3)/2*Lx;
R = rcell*Lx/n; % disc radius
r_av = sqrt(Lx*Ly/N)/2; % estimate for NND of random distribution (Clark & Evans, 1954)
eta = N*pi*R^2/(Lx*Ly); %packing fraction 
if eta > pi/(2*sqrt(3)) % max. packing fraction
    disp('packing fraction too high! Abort evaluation');
    pause(10);
end
% Circuit Parameters
K_all = 16; %[6 12 20 16];
Con_all = 8; %[15 15 15 8];
noise = 0;
a0_all = 0.5; %[1.5 1.5 1.5 0.5];
rcell = 0.2;

% initial lattice options
choice = 2; %1: random, 2: Markov MC
mcsteps_all = 0; %[0 10.^[1 3] ]; %10^3;

K=K_all;
Con=Con_all;
a0=a0_all;

for idx_param = 1:numel(mcsteps_all)
    %K = K_all(idx_param);
    %Con = Con_all(idx_param);
    %a0 = a0_all(idx_param);
    mcsteps = mcsteps_all(idx_param);
    
    % save data file name
    labels = {'RandomPlacement', sprintf('MarkovMC_%d_MCsteps', mcsteps)};
    fname_str = strrep(sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%.1f_%s', ...
        N, Con, K, a0, labels{choice}), '.', 'p');
    path_out_mat = 'H:\My Documents\Multicellular automaton\data\random_positions\pin-pout';
    fname_out = fullfile(path_out_mat, strcat(fname_str,'.mat'));
    
    % save folder for figures
    path_out_fig = 'H:\My Documents\Multicellular automaton\figures\random_positions';
    
    % try to load the map
    if exist(fname_out, 'file')==2
        load(fname_out, 'count', 't_av');
    else
        % calculate the map
        disp('Does not exist!');
        %
        n_smpl = 1000;
        switch choice
            case 1
                [count, t_av] = count_eq_parallel_rand(n,...
                    Con, K, a0, rcell, noise);
            case 2
                [count, t_av] = count_eq_parallel_markovMC_temp(n,...
                    Con, K, a0, rcell, noise, mcsteps, n_smpl);
        end
        %}
    end
    
    % post-processing
    prob = transpose(count./repmat(sum(count,2),1,N+1));
    
    % save .mat file
    %save(fname_out);

    %%
    % plot the map
    h1 = figure();
    hold on
    p = (0:N)./N;
    im_fig = imagesc(p,p,prob);
    % set title and font
    %title(sprintf('N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f', ...
    %    N, K, Con, a0, Rcell),'FontSize', 18)
    set(gca,'Ydir', 'normal', 'FontSize', 20)
    % set invisible parts where count is zero
    set(im_fig, 'AlphaData', count' > 0);
    % set colorbar and labels
    c = colorbar;
    c.Label.String = 'Probability';
    xlabel('$$p_{in}$$', 'FontSize', 24)
    ylabel('$$p_{out}$$', 'FontSize', 24)
    xlim([0 1]);
    ylim([0 1]);
    
    % Plot line p = 1/2 - 4B/fN
    %{
    B = (Con+1)/2*(1+fN) - K;
    line = 1/2 - B/(4*fN);
    plot([line line], [0 1], 'r--');
    %}
    
    % Organize and save
    save_fig = 0; % save figure? 0:no, 1: yes
    if save_fig > 0
        out_file = fullfile(path_out_fig, strcat(fname_str,'_map'));
        save_figure(h1, 10, 8, out_file, '.pdf');
    end
    %%
    %{
    % Plot the average number of steps it takes to reach equilibrium
    h2 = figure();
    plot(p, t_av, 'r-o')
    set(gca,'FontSize', 20)
    %title(sprintf('N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f', ...
    %    N, K, Con, a0, Rcell),'FontSize', 18)
    xlabel('$$p_{in}$$', 'FontSize', 24)
    ylabel('$$\langle t_{eq} \rangle $$', 'FontSize', 24)

    save_fig = 1; % save figure? 0:no, 1: yes
    if save_fig > 0
        out_file = fullfile(path_out_fig, strcat(fname_str,'_teq_av'));
        save_figure(h2, 10, 8, out_file, '.pdf');
    end
    %}
    %% 
    %close all
end