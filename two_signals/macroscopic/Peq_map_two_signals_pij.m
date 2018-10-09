%% Plots a map of Peq in p, I space
% Using the new approach with p_ij instead of taking genes separately
clear all
close all
clc
warning off
set(0, 'defaulttextinterpreter', 'latex');
%% Parameters
%
% lattice parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters 
M_int = [0 -1; 1 0];
Con = [18 16];
Coff = [1 1];
K = [0 15; 15 0];% K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
hill = Inf;
%}

% Load from file
%{
%folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\sample_trajectories\typeII';
folder = 'H:\My Documents\Multicellular automaton\data\two_signals\transition_matrix\sim';
fname_str = 'W_sim_M_int_0_1_-1_-1_N225_a0_1p50_rcell_0p20_K_0_12_13_8_Con_18_16_iniON_50_61_58_56_delta_t_1_trials_1000';
load(fullfile(folder,fname_str), 'gz', 'N', 'a0', 'rcell', 'Rcell', 'M_int',...
    'Con', 'Coff', 'K', 'lambda', 'hill', 'noise', 'iniON');
%}

% load dist, pos
% [dist, pos] = init_dist_hex(gz, gz);
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
%% Check custom result
pij = [0 0.3; 0.7 0];
W = transition_prob_two_signals_pij(M_int, Con, Coff, K, fN, gN, pij);
%disp(diag(W))
iniON = reshape(round(pij.*N),4,1); % iniON: (0,0), (1,0), (0,1), (1,1)
Peq = prod(diag(W).^iniON); % diag(W): W(0,0), W(1,0), W(0,1), W(1,1)
disp(Peq)

% Plot W
%
h = figure(1);
imagesc(1:4, 1:4, W')
hold on
for i=1:3
    line([0.5, 4.5], [i+0.5, i+0.5], 'Color', 'r', 'LineWidth', 1);
    line([i+0.5, i+0.5], [0.5, 4.5], 'Color', 'r', 'LineWidth', 1);
end
c = colorbar;
set(gca, 'Ydir', 'normal', 'FontSize', 20)
xlabel('$$X(t)$$', 'FontSize', 24)
ylabel('$$X(t+1)$$', 'FontSize', 24)
c.Label.String = 'Probability';
set(gca, 'XTick', 1:4, 'YTick', 1:4, 'XTickLabels', ...
    {'(0,0)','(1,0)','(0,1)','(1,1)'},...
    'YTickLabels', {'(0,0)','(1,0)','(0,1)','(1,1)'});
set(h, 'Position', [100 100 640 520]);
caxis([0 1]);
%}
%% Calculate P_eq
% Version where we consider the two genes separately
p = (0:N)/N;
%I = 0;
%I = linspace(0, 1, 100);
%[pm, Im] = meshgrid(p,I);

%Peq = zeros(numel(p), numel(p), numel(p), numel(p));
non_trivial_Peq_idx = {};
non_trivial_Peq_vals = [];
ip1 = 1; 
iI1 = 1;
for i00=N-25:N
    for i01=0:N-i00
        for i10=0:N-i00-i01
            fprintf('%d %d %d %d \n', i00, i01, i10, N-i00-i01-i10);
            pij = [i00 i01; i10 N-i00-i01-i10]./N;
            W = transition_prob_two_signals_pij(M_int, Con, Coff, K, fN, gN, pij);
            
            iniON = reshape(round(pij.*N)',4,1); % iniON: (0,0), (0,1), (1,0), (1,1)
            this_Peq =  prod(diag(W).^iniON);
            disp(this_Peq);
            if this_Peq>0
                non_trivial_Peq_idx{end+1} = [i00 i01 i10 N-i00-i01-i10];
                non_trivial_Peq_vals(end+1) = this_Peq;
            end
            %Peq(ip1, ip2, iI1, iI2) = pEn;
        end
    end
end

%% Show non-trivial Peq values/indices
celldisp(non_trivial_Peq_idx)
%celldisp(non_trivial_Peq_vals)