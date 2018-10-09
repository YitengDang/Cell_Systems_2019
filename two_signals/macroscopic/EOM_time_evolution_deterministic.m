%% Calculates the time evolution of the equation of motion (EOM) in terms of p^{(i,j)}(t)
% Deterministic evolution (Markov process) using W as transition matrix
clear variables
close all
clc
set(0, 'defaulttextinterpreter', 'latex');
%warning off

%% System parameters
%{
% Parameters of the system
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

M_int = [0 1; -1 -1];
Con = [18 16];
Coff = [1 1];
K = [0 12; 13 8];
lambda = [1 1.2];
hill = Inf;
noise = 0;

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
%}

%% Load parameters from saved data
%
%folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\macroscopic';
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\macroscopic\period3';
%folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\macroscopic\chaotic';
fname_str = 'N225_a0_1p80_rcell_0p20_lambda12_1p20_M_int_0_1_-1_-1_K_0_12_13_8_Con_18_16_p0_0p25_0p25_0p25_0p25_t_out_20_EOM_determ';

fname = fullfile(folder, fname_str);
load(fname);

%p0 = p_ini;
p0 = iniON/N;
%tmax = t_out + 1;
%tmax = 100;
%}
%% Simulate trajectory deterministic
%
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
tmax = 50; %t_out;
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

% Test periodicity
%{
round(pij_all(:,:, end-4)*N)
round(pij_all(:,:, end-3)*N)
round(pij_all(:,:, end-2)*N)
round(pij_all(:,:, end-1)*N)
round(pij_all(:,:, end)*N)
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
%iniON_s = sprintf('%d_%d_%d_%d', iniON(1,1), iniON(1,2), iniON(2,1), iniON(2,2));
p0_s = sprintf('%.2f_%.2f_%.2f_%.2f', p0(1,1), p0(1,2), p0(2,1), p0(2,2));
SimID = 'EOM_determ';
if tmax~=t_out+1
    t_s = num2str(tmax);
else
    t_s = sprintf('tmax%d', tmax);
end

fname_str = strrep(sprintf('N%d_a0_%s_rcell_%s_lambda12_%s_M_int_%s_K_%s_Con_%s_p0_%s_t_out_%s_%s',...
    N, a0_s, R_s, lambda_s, M_int_s, K_s, Con_s, p0_s, t_s, SimID), '.', 'p');

% Save plot
qsave = 1;
folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\macroscopic';
%folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\two_signals_macroscopic';
save_figure(h1, 10, 8, fullfile(folder, fname_str), '.pdf', qsave);

% save data
qsave = 1;
if qsave
    close all
    folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\macroscopic';
    %folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\two_signals_macroscopic';
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