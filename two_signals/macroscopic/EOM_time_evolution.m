% Calculates the time evolution of the equation of motion (EOM) in terms of p^{(i,j)}(t)
% Use a Monte Carlo algorithm (stochastic)
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'latex');
%warning off

%% System parameters
% lattice parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters 
M_int = [0 1; -1 1];
Con = [18 16];
Coff = [1 1];
K = [0 15; 11 4];% K(i,j): sensitivity of type i to type j molecules
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
%% Simulate trajectory deterministic
%{
% variables
maxsteps = 1000;
pij_all = zeros(2, 2, maxsteps+1);
nij_all = zeros(2, 2, maxsteps+1);
% initial state
p0 = [0.2 0.3; 0.25 0.25];
I0 = [0 0];

% Simulate
pij_all(:,:,1) = p0;
nij_all(:,:,1) = round(p0*N); 

p_in = p0;
I_in = I0;
for t=1:maxsteps
    %disp(t);
    %disp(p_in);
    %disp(round(p_in*N));
    [p_out, Peq] = EOM_update_markov(N, M_int, Con, Coff, K, fN, gN, p_in, I_in);
    p_in = p_out;

    % system in equilibrium?
    if rand < Peq
        pij_all(:,:,t+1:end) = repmat(p_out, 1, 1, maxsteps-t+1);
        nij_all(:,:,t+1:end) = repmat(round(N*p_out), 1, 1, maxsteps-t+1);
        break
    else
        pij_all(:,:,t+1) = p_out;
        nij_all(:,:,t+1) = round(N*p_out);
    end
end
t_out = t; %default value
%}

%% Plot trajectory
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
    plot(0:tmax, squeeze(nij_all(i,j,1:tmax+1))/N, ps, 'Color', clr, 'LineWidth', lw);
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
pini_s = sprintf('%.2f_%.2f_%.2f_%.2f', p_ini(1,1), p_ini(1,2), p_ini(2,1), p_ini(2,2));
simID = 'EOM_stoch';

if tmax~=t_out+1
    t_s = num2str(tmax);
else
    t_s = 'tmax';
end
fname_str = strrep(sprintf('N%d_a0_%s_rcell_%s_lambda12_%s_M_int_%s_K_%s_Con_%s_iniON_%s_t_out_%s_%s',...
    N, a0_s, R_s, lambda_s, M_int_s, K_s, Con_s, iniON_s, t_s, simID), '.', 'p');

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
%% Updating step Markov (tester)
%{
for i=1:10
    p_in = [0.25 0.25; 0.25 0.25];
    I_in = [0 0];
    [p_out, Peq] = EOM_update_markov(N, M_int, Con, Coff, K, fN, gN, p_in, I_in);
    disp(p_out)
    %disp(Peq)
end
%}
%% Updating step MC (tester)
%{
%n_in = [0 7; 10 13];
%p0 = n_in/N;

p_ini = [0.25 0.25; 0.25 0.25];
%p_ini = [0.8 0.2; 0 0];
I_in = [0 0];
n_in = round(p_ini*N);

% get transition matrix
W = transition_prob_two_signals_pij(M_int, Con, Coff, K, fN, gN, p_ini, I_in);
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
for t=1:maxsteps
    %disp(t);
    [p_out, Peq] = EOM_update_markov(N, M_int, Con, Coff, K, fN, gN, p_in, I_in);
    p_in = p_out;
    
    p_all(:,:,t+1) = p_out;
    
    % Peq not needed in principle
    %
    if rand < Peq
        p_all(:,:,t+1:end) = repmat(p_out, 1, 1, maxsteps-t+1);
        break
    else
        p_all(:,:,t+1) = p_out;
    end
    %
end
t_out = t;
%}