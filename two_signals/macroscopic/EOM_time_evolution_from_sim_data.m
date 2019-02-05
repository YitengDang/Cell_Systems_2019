%% Calculates the time evolution of the equation of motion (EOM) in terms of p^{(i,j)}(t)
% Use a Monte Carlo algorithm (stochastic)
% Uses parameters and initial conditions from simulation results
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'tex');

%% Load simulation
%{
folder = 'H:\My Documents\Multicellular automaton\temp';
fname_str = 'plane_wave_formation_period_15';
fname = fullfile(folder, fname_str);
%}
%folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\macroscopic\sweep_K12_N225\simulations';
%folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\filtered_TW';
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\patterns';
[file, folder] = uigetfile(fullfile(folder, '*.mat'));
fname = fullfile(folder, file);

load(fname, 'save_consts_struct', 'cells_hist')
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
if sum(strcmp('mcsteps', fields(s)))
    mcsteps = str2double(s.mcsteps);
else
    mcsteps = 0;
end

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
cells_ini_idx = cells_ini*[1; 2]; % 0=(0,0), 1=(1,0), 2=(0,1), 3=(1,1)
iniON = reshape(histcounts(cells_ini_idx, -0.5:3.5), 2, 2); % [1,1]=(0,0), [1,2] = (0,1), [2,1]=(1,0), [2,2]=(1,1)
I_in = zeros(1, 2);
%I_in(1) = moranI(cells_ini(:,1), a0*dist);
%I_in(2) = moranI(cells_ini(:,2), a0*dist);
fprintf('I_ini = [%.3f %.3f] \n', I_in(1), I_in(2) );
%% For saving
% filename for saving
M_int_s = sprintf('%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
a0_s = sprintf('%.2f', a0);
R_s = sprintf('%.2f', rcell);
K_s = sprintf('%d_%d_%d_%d', K(1,1), K(1,2), K(2,1), K(2,2));
Con_s = sprintf('%d_%d', Con(1), Con(2));
lambda_s = sprintf('%.2f', lambda(2));
iniON_s = sprintf('%d_%d_%d_%d', iniON(1,1), iniON(1,2), iniON(2,1), iniON(2,2));

% save folder
%save_data_folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\two_signals_macroscopic';
%save_data_folder = 'H:\My Documents\Multicellular automaton\temp';
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\macroscopic';
subfolder = 'sweep_K12_N225';
save_folder = fullfile(folder, subfolder);
%% Plot simulation trajectory
tmax_sim = numel(cells_hist)-1;
%tmax_sim = 1000;
pij_all_sim = zeros(tmax_sim+1, 4);
for i=1:numel(cells_hist)
    cells = cells_hist{i};
    cells_idx = cells*[1; 2];
    for ii=1:4
        pij_all_sim(i, ii) = sum(cells_idx==ii-1)/N;
    end
end

% special extension for period 4 oscillation (ad hoc)
%{
conv_table = [4 5 2 3]; %[0 1 2 3] -> [4 5 2 3]
for i=numel(cells_hist)+1:tmax_sim+1
    earlier_state = conv_table( mod(i, 4)+1 );
    pij_all_sim(i, :) = pij_all_sim(earlier_state, :) ;
end
%}

h=figure;
% plot style
plot_clrs = [1 1 1; 
            1 1 0
            0 0 1
            0 0 0];
if tmax_sim < 100
    ps = 'o-'; lw = 2;
elseif tmax_sim < 500
    ps = '.-'; lw = 1.5;
else
    ps = '.-'; lw = 0.75;
end

hold on
for idx=1:4
    clr = plot_clrs(idx, :);
    %[i,j] = ind2sub(2, idx);
    plot(0:tmax_sim, pij_all_sim(:, idx), ps, 'Color', clr, 'LineWidth', lw);
    %plot(tmin:tmax+1, squeeze(pij_all(i,j,tmin:tmax+1)), ps, 'Color', clr, 'LineWidth', lw);
end
xlabel('t', 'Interpreter', 'tex');
ylabel('p^{(i,j)}', 'Interpreter', 'tex');
set(gca, 'Color', [0.8 0.8 0.8]);   
set(gca, 'FontSize', 24);
xlim([0 tmax_sim]);
ylim([0 1]);
%legend_text = sprintfc("(%d, %d)", [0 0; 1 0; 0 1; 1 1]);
%legend(legend_text, 'Location', 'eastoutside', 'FontSize', 12);
fig_pos = [1 1 7 5];
set(h, 'Units', 'inches', 'Position', fig_pos);

% save figure
simID = 'simulation';
fname_str = strrep(sprintf('N%d_a0_%s_rcell_%s_lambda12_%s_M_int_%s_K_%s_Con_%s_iniON_%s_%s',...
    N, a0_s, R_s, lambda_s, M_int_s, K_s, Con_s, iniON_s, simID), '.', 'p');
qsave = 0;
%folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\two_signals_macroscopic';
colored_background = 1;
save_figure(h, 10, 8, fullfile(save_folder, fname_str), '.pdf', qsave, colored_background);    

%% Simulate trajectory Monte Carlo
nsim = 5;
t_out_all = zeros(nsim, 1);
pij_out_all = zeros(nsim, 4);

for i=1:nsim
    % variables
    maxsteps = 1000;
    n_all = zeros(2, 2, maxsteps+1);

    % Simulate
    n_all(:,:,1) = iniON;
    n_in = iniON;
    t_out = maxsteps-1; %default value
    for t=1:maxsteps
        % take I values from trajectory
        %{
        cells = cells_hist{ min(t, tmax_sim) };
        I_in(1) = moranI(cells(:,1), a0*dist);
        I_in(2) = moranI(cells(:,2), a0*dist);
        %}
        
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
    
    % Store data
    t_out_all(i, 1) = t_out;
    pij_out_all(i, :) = n_out(:)/N; % 1 = (0,0), 2 = (0,1), 3 = (1,0), 4 = (1,1)
    %}
    %% Plot trajectory
    %
    tmin = 0;
    tmax = t_out;
    fig_pos = [1 1 7 5];

    h1=figure;
    % plot style
    plot_clrs = [1 1 1; 
                1 1 0
                0 0 1
                0 0 0];
    if tmax < 100
        ps = 'o-'; lw = 2;
    elseif tmax < 500
        ps = '.-'; lw = 1.5;
    else
        ps = '.-'; lw = 0.75;
    end

    hold on
    for idx=1:4
        clr = plot_clrs(idx, :);
        [i,j] = ind2sub(2, idx);
        plot(0:tmax, squeeze(n_all(i,j,1:tmax+1))/N, ps, 'Color', clr, 'LineWidth', lw);
        %plot(tmin:tmax+1, squeeze(pij_all(i,j,tmin:tmax+1)), ps, 'Color', clr, 'LineWidth', lw);
    end
    xlabel('t', 'Interpreter', 'tex');
    ylabel('p', 'Interpreter', 'tex');
    set(gca, 'Color', [0.8 0.8 0.8]);   
    set(gca, 'FontSize', 24);
    xlim([tmin tmax]);
    ylim([0 1]);
    %legend_text = sprintfc("(%d, %d)", [0 0; 1 0; 0 1; 1 1]);
    %legend(legend_text, 'Location', 'eastoutside', 'FontSize', 12);
    set(h1, 'Units', 'inches', 'Position', fig_pos);
    %}
    %% Save 
    %{
    % save data and plot
    %
    if tmax~=t_out+1
        t_s = num2str(tmax);
    else
        t_s = 'tmax';
    end
    %}
    simID = 'EOM_stoch';
    fname_str = strrep(sprintf('N%d_a0_%s_rcell_%s_lambda12_%s_M_int_%s_K_%s_Con_%s_iniON_%s_%s',...
        N, a0_s, R_s, lambda_s, M_int_s, K_s, Con_s, iniON_s, simID), '.', 'p');
    
    qsave = 0;
    if qsave
        % ---- Save data ----
        i=1;
        subfolder = 'EOM_data';
        fname_str_version = strcat(fname_str, '_EOM-v', num2str(i));
        fname = fullfile(save_folder, subfolder, strcat(fname_str_version, '.mat'));
        while exist(fname, 'file')==2
            i=i+1;
            fname_str_version = strcat(fname_str, '_EOM-v', num2str(i));
            fname = fullfile(save_folder, subfolder, strcat(fname_str_version, '.mat'));
        end
        save(fname, 'n_all', 'iniON', 'I_in', 'maxsteps', 't_out', 'save_consts_struct');
        fprintf('Saved trajectory as %s \n', fname);
        
        % ---- Save plot ----
        qsave = 0;
        %folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\two_signals_macroscopic';
        colored_background = 1;
        save_figure(h1, 10, 8, fullfile(save_folder, fname_str_version), '.pdf', qsave, colored_background);
        %}
    end
    %}
end
%% Compare t_out
figure;
plot([tmax_sim tmax_sim], [0 1], 'r--');
hold on
edges = 0:50:maxsteps;
histogram(t_out_all, edges, 'Normalization', 'probability');
xlim([0 maxsteps]);
ylim([0 1]);
xlabel('t_{out}');

%%
figure;
hold on
scatter( 1:4, pij_all_sim(end, :), 'r', 'filled'  );
plot( 1:4, pij_all_sim(end, :), 'r--' );
for i=1:nsim
    scatter( 1:4, pij_out_all(i,:), 75, 'b' );
    plot(1:4, pij_out_all(i,:), 'b--', 'LineWidth', 0.5);
end
xlim([0.5 4.5]);
ylim([0 1]);
set(gca, 'XTick', 1:4);
