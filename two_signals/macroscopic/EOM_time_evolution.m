%% Calculates the time evolution of the equation of motion (EOM) in terms of p^{(i,j)}(t)
% (1) Use a Monte Carlo algorithm (stochastic)
% (2) Deterministic evolution (Markov process) using W as transition matrix
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'latex');
%warning off

%% System parameters
% lattice parameters
gz = 15;
N = gz^2;
a0 = 0.5;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters 
M_int = [0 1; -1 1];
Con = [18 16];
Coff = [1 1];
K = [0 35; 30 10];% K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
hill = Inf;

% Load from saved data

% (1) original hexagonal lattice
%[pos,ex,ey] = init_cellpos_hex(gz,gz);
%dist = dist_mat(pos,gz,gz,ex,ey);
% (2) new lattice
mcsteps = 0;
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

%% Updating step MC (tester)
%{
%n_in = [0 7; 10 13];
%p0 = n_in/N;

p_ini = [0.25 0.25; 0.25 0.25];
%p_ini = [0.8 0.2; 0 0];
I_in = [0 0];
n_in = round(p_ini*N);

% get transition matrix
W = transition_prob_states_two_signals_analytical_calc(M_int, Con, Coff, K, fN, gN, p_ini, I_in);
% (0,0), (0,1), (1,0), (1,1)

disp(W);

n_out = zeros(2);
for k=1:4
    %disp(n_in(k));
    p = [0 cumsum(W(k, :), 2)];
    r = rand(n_in(k), 1);
    for i=1:n_in(k) % can we vectorize this?
        out_idx = find(r(i) < p, 1)-1;
        n_out(out_idx) = n_out(out_idx) + 1;
    end
end

disp( n_out )
%disp( EOM_update_MC(N, M_int, Con, Coff, K, fN, gN, n_in) )
%}

%% Simulate trajectory Monte Carlo
%
% variables
maxsteps = 1000;
n_all = zeros(2, 2, maxsteps+1);

% initial state
p_ini = [0.25 0.25; 0.25 0.25];
%p_ini = [0.4 0.1; 0.1 0.4];
I_in = [0 0];

% convert input p^ij to input iniON
%
iniON = round(p_ini*N);
d = sum(p_ini(:)*N)- sum(iniON(:));
idx = randperm(4, 1);
iniON(idx) = iniON(idx) + d;
if abs(d)>2
    warning('Wrong input n!');
end
%
%iniON = [11 34; 34 146];

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

%% Updating step Markov
%{
for i=1:10
    p_in = [0.25 0.25; 0.25 0.25];
    I_in = [0 0];
    [p_out, Peq] = EOM_update_markov(N, M_int, Con, Coff, K, fN, gN, p_in, I_in);
    disp(p_out)
    %disp(Peq)
end
%}
%% Simulate trajectory Markov chain
%{
% variables
maxsteps = 100;
p_all = zeros(2, 2, maxsteps+1);

% initial state
p_ini = [0.25 0.25; 0.25 0.25];
%p_ini = [0.4 0.1; 0.1 0.4];
I_in = [0 0];

% Simulate
p_all(:,:,1) = p_ini;
p_in = p_ini;
t_out = maxsteps-1; %default value
for t=1:maxsteps
    [p_out, Peq] = EOM_update_markov(N, M_int, Con, Coff, K, fN, gN, p_in, I_in);
    p_in = p_out;
    
    p_all(:,:,t+1) = p_out;
    
    % Peq not needed in principle
    %{
    if rand < Peq
        p_all(:,:,t+1:end) = repmat(p_out, 1, 1, maxsteps-t+1);
        t_out = t;
        break
    else
        p_all(:,:,t+1) = p_out;
    end
    %}
end
%}
%% Plot trajectory
tmax = 100; % t_out+1;
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
    %plot(0:tmax, squeeze(p_all(i,j,1:tmax+1)), ps, 'Color', clr, 'LineWidth', lw);
end
xlabel('$$t$$', 'Interpreter', 'latex');
ylabel('$$p$$', 'Interpreter', 'latex');
set(gca, 'Color', [0.8 0.8 0.8]);   
set(gca, 'FontSize', 24);
xlim([0 tmax]);
ylim([0 1]);
legend_text = sprintfc("(%d, %d)", [0 0; 1 0; 0 1; 1 1]);
legend(legend_text, 'Location', 'eastoutside', 'FontSize', 12);
set(h1, 'Units', 'inches', 'Position', fig_pos);

%% Save 
%{
% Save data
% filename
M_int_s = sprintf('%d_%d_%d_%d', M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2));
a0_s = sprintf('%.2f', a0);
R_s = sprintf('%.2f', rcell);
K_s = sprintf('%d_%d_%d_%d', K(1,1), K(1,2), K(2,1), K(2,2));
Con_s = sprintf('%d_%d', Con(1), Con(2));
lambda_s = sprintf('%.2f', lambda(2));
iniON_s = sprintf('%d_%d_%d_%d', iniON(1,1), iniON(1,2), iniON(2,1), iniON(2,2));
pini_s = sprintf('%.2f_%.2f_%.2f_%.2f', p_ini(1,1), p_ini(1,2), p_ini(2,1), p_ini(2,2));

if tmax~=t_out+1
    t_s = num2str(tmax);
else
    t_s = 'tmax';
end
fname_str = strrep(sprintf('N%d_a0_%s_rcell_%s_lambda12_%s_M_int_%s_K_%s_Con_%s_iniON_%s_t_out_%s_EOM',...
    N, a0_s, R_s, lambda_s, M_int_s, K_s, Con_s, iniON_s, t_s), '.', 'p');

% Save plot
qsave = 0;
%folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\macroscopic';
folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\two_signals_macroscopic';
save_figure(h1, 10, 8, fullfile(folder, fname_str), '.pdf', qsave);

% save data
qsave = 0;
if qsave
    close all
    folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\two_signals_macroscopic';
    save(fullfile(folder, fname_str));
end
%}
%% Loop over trajectories
%--------Not updated---------
%{
p_stoch = zeros(4, nruns, maxsteps+1);
I_stoch = zeros(2, nruns, maxsteps+1);
pnoise = zeros(nruns, maxsteps+1);
Inoise = zeros(nruns, maxsteps+1);
eq_times2 = zeros(nruns, 1);

% If plotting multiple trajectories from same p_ini
%test = reshape(repmat(round(N*p_ini)/N, [10 1]), [30 1]);
%I_ini = 0;

pnoise(:,1) = 0;
Inoise(:,1) = 0; %normrnd(0, 0.02, [nruns, 1]); % take larger noise for initial spread
p_stoch(:,1) = p_ini' + pnoise(:,1);
%p_stoch(:,1) = test;
I_stoch(:,1) = I_ini + Inoise(:,1);

Peq_all = zeros(nruns, maxsteps);
for i=1:nruns
    for j=1:maxsteps 
        pt = p_stoch(i, j);
        It = I_stoch(i, j);
        [p_new, I_new, Peq, p_lim] = EOM_update_Langevin(pt, It, N, Con, K, fN, gN, delta, sigma, 0);
        Peq_all(i,j) = Peq;
        
        % system out of range?
        if p_lim == 1 || p_lim == -1
            p_stoch(i, j+1:end) = (p_lim+1)/2;
            I_stoch(i, j+1:end) = 0;
        end
        
        % system in equilibrium?
        if rand < Peq
            eq_times2(i) = j-1;
            p_stoch(i, j+1:end) = pt;
            I_stoch(i, j+1:end) = It;
            break
        else
            p_stoch(i, j+1) = p_new;
            I_stoch(i, j+1) = I_new;
        end
    end
end

% Plot Langevin trajectories
figure(h);
hold on
%h1a = plot(p_stoch', I_stoch', 'g-', 'LineWidth', 1.5);
handles = cell(numel(p_ini), 1);
for pl=1:nruns
    handles{pl} = plot(p_stoch(pl,:)', I_stoch(pl,:)', 'g-', 'LineWidth', 2);
    handles{pl}.Color(4) = 0.6;
end   

% plot start and end points
plot(p_stoch(:,1), I_stoch(:,1), 'go', 'LineWidth', 2);
plot(p_stoch(:,end), I_stoch(:,end), 'gx', 'LineWidth', 2);

% set figure properties
set(gca, 'FontSize', 40);
xticks([0 1]);
yticks([0 1]);
set(c, 'XTick', [0 1]);
%set(gcf, 'CLim', [0 1])
%}
%% Save figure
%{
qsave = 0;
if qsave
    save_path = 'H:\My Documents\Multicellular automaton\figures\random_positions';
    fname_str = strrep(sprintf('Peq_trajectories_N%d_a0_%.1f_K_%d_Con_%d_noise%.1f_%s', ...
        N, a0, K, Con, noise, labels{choice}), '.', 'p');
    fname = fullfile(save_path, fname_str);
    save_figure(h, 10, 8, fname, '.pdf');
end

%}