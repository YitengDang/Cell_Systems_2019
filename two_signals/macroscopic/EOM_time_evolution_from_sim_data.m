%% Calculates the time evolution of the equation of motion (EOM) in terms of p^{(i,j)}(t)
% Use a Monte Carlo algorithm (stochastic)
% Uses parameters and initial conditions from simulation results
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'latex');

%% Load simulation
folder = 'H:\My Documents\Multicellular automaton\temp';
fname_str = 'plane_wave_formation_period_15';
fname = fullfile(folder, fname_str);
load(fname);

%% System parameters
s = save_consts_struct;
N = s.N;
a0 = s.a0;
K = s.K;
Con = s.Con;
Coff = s.Coff;
M_int = s.M_int;
hill = s.hill;
noise = s.noise;
rcell = s.rcell;
%cells = cells_hist{1};
lambda12 = s.lambda12;
lambda = [1 lambda12];
mcsteps = str2double(s.mcsteps);

%p0 = s.p_ini;
%tmax =  s.tmax;
gz = sqrt(N);
Rcell = rcell*a0;

% (2) new lattice
%ccsteps = 0;
[pos, dist] = initial_cells_random_markov_periodic(gz, mcsteps, rcell);

% Calculate interaction strength
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = zeros(1, 2);
gN = zeros(1, 2);
for i=1:2
    fN(i) = sum(sinh(Rcell) * sum(exp((Rcell-r)./lambda(i)).*(lambda(i)./r)) ); % calculate signaling strength
    gN(i) = sum(sinh(Rcell)^2 * sum(exp(2*(Rcell-r)./lambda(i)).*(lambda(i)./r).^2 ) ); % calculate noise variance strength
end

% initial state from sim
cells_ini = cells_hist{1};
cells_ini_idx = cells_ini*[1; 2];
iniON = reshape(histcounts(cells_ini_idx, -0.5:3.5), 2, 2);
I_in = zeros(1, 2);
I_in(1) = moranI(cells_ini(:,1), a0*dist);
I_in(2) = moranI(cells_ini(:,2), a0*dist);

%% Simulate trajectory Monte Carlo
nsim = 5;
for i=1:nsim
    % variables
    maxsteps = 1000;
    n_all = zeros(2, 2, maxsteps+1);

    % Simulate
    n_all(:,:,1) = iniON;
    n_in = iniON;
    t_out = maxsteps-1; %default value
    for t=1:maxsteps
        [n_out, Peq] = EOM_update_MC(N, M_int, Con, Coff, K, fN, gN, n_in, I_in);
        n_in = n_out;

        % system in equilibrium?
        if rand < Peq
            n_all(:,:,t+1:end) = repmat(n_out, 1, 1, maxsteps-t+1);
            t_out = t;
            break
        else
            n_all(:,:,t+1) = n_out;
        end
    end
    %}
    %% Plot trajectory
    %{
    tmin = 0;
    tmax = 100; % t_out;
    fig_pos = [1 1 7 5];

    h1=figure;
    % plot style
    plot_clrs = [1 1 1; 
                1 1 0
                0 0 1
                0 0 0];
    if tmax < 100
        ps = 'o-'; lw = 1.5;
    elseif tmax < 500
        ps = '.-'; lw = 1;
    else
        ps = '.-'; lw = 0.5;
    end

    hold on
    for idx=1:4
        clr = plot_clrs(idx, :);
        [i,j] = ind2sub(2, idx);
        plot(0:tmax, squeeze(n_all(i,j,1:tmax+1))/N, ps, 'Color', clr, 'LineWidth', lw);
        %plot(tmin:tmax+1, squeeze(pij_all(i,j,tmin:tmax+1)), ps, 'Color', clr, 'LineWidth', lw);
    end
    xlabel('$$t$$', 'Interpreter', 'latex');
    ylabel('$$p$$', 'Interpreter', 'latex');
    set(gca, 'Color', [0.8 0.8 0.8]);   
    set(gca, 'FontSize', 24);
    xlim([tmin tmax]);
    ylim([0 1]);
    legend_text = sprintfc("(%d, %d)", [0 0; 1 0; 0 1; 1 1]);
    legend(legend_text, 'Location', 'eastoutside', 'FontSize', 12);
    set(h1, 'Units', 'inches', 'Position', fig_pos);
    %}
    %% Save 
    %
    % Save data
    % filename
    M_int_s = sprintf('%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
    a0_s = sprintf('%.2f', a0);
    R_s = sprintf('%.2f', rcell);
    K_s = sprintf('%d_%d_%d_%d', K(1,1), K(1,2), K(2,1), K(2,2));
    Con_s = sprintf('%d_%d', Con(1), Con(2));
    lambda_s = sprintf('%.2f', lambda(2));
    iniON_s = sprintf('%d_%d_%d_%d', iniON(1,1), iniON(1,2), iniON(2,1), iniON(2,2));
    simID = 'EOM_stoch';
    
    %{
    if tmax~=t_out+1
        t_s = num2str(tmax);
    else
        t_s = 'tmax';
    end
    fname_str = strrep(sprintf('N%d_a0_%s_rcell_%s_lambda12_%s_M_int_%s_K_%s_Con_%s_iniON_%s_t_out_%s_%s',...
        N, a0_s, R_s, lambda_s, M_int_s, K_s, Con_s, iniON_s, t_s, simID), '.', 'p');
        
    % Save plot
    qsave = 0;
    folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\macroscopic';
    %folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\two_signals_macroscopic';
    save_figure(h1, 10, 8, fullfile(folder, fname_str), '.pdf', qsave);
    %}
    
    
    % save data
    qsave = 1;
    if qsave
        close all
        %folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\two_signals_macroscopic';
        folder = 'H:\My Documents\Multicellular automaton\temp';
        i=1;
        fname = fullfile(folder, strcat(fname_str, '_EOM-v', num2str(i), '.mat'));
        while exist(fname, 'file')==2
            i=i+1;
            fname = fullfile(folder, strcat(fname_str, '_EOM-v', num2str(i), '.mat'));
        end
        save(fname, 'n_all', 'iniON', 'I_in', 'maxsteps', 't_out', 'save_consts_struct');
        fprintf('Saved trajectory as %s \n', fname);

    end
    %}
end
