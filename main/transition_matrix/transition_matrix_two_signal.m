% Script to calculate the transition matrix in different conditions
% Without I: only takes into account initial p
clear variables
close all
warning off
set(0, 'defaulttextinterpreter', 'latex');
maxNumCompThreads(5);

% Parameters of the system
gz = 10;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

M_int = [0 1; -1 1];
Con = [18 16];
Coff = [1 1];
K = [0 15; 11 4];
lambda = [1 1.2];

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gz,gz);
dist = dist_mat(pos,gz,gz,ex,ey);

% Calculate interaction strength
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = zeros(1, 2);
gN = zeros(1, 2);
for i=1:2
    fN(i) = sum(sinh(Rcell) * sum(exp((Rcell-r)./lambda(i)).*(lambda(i)./r)) ); % calculate signaling strength
    gN(i) = sum(sinh(Rcell)^2 * sum(exp(2*(Rcell-r)./lambda(i)).*(lambda(i)./r).^2 ) ); % calculate noise variance strength
end

%% Load/save folder
fname_str = strrep(sprintf('t_mat_gz_%d_M_int_%d_%d_%d_%d_Con_%d_%d_K_%d_%d_%d_%d_a0_%.2f_rcell_%.2f', ...
    gz, M_int(1,1), M_int(1,2), M_int(2,1), M_int(2,2), Con,...
    K(1,1), K(1,2), K(2,1), K(2,2), a0, rcell), '.', 'p');

folder = 'H:\My Documents\Multicellular automaton\data\two_signals\transition_matrix\';
fname = fullfile(folder, strcat(fname_str, '.mat'));

if exist(fname_str, 'file') == 2
    load(fname)
else
    t_mat = zeros(N+1, N+1, N+1, N+1); % 4D transition matrix
    for n1 = 0:N
        for n2 = 0:N
            n = [n1 n2];
            disp(n)
            [ptsum, ~, ~] = transition_prob_two_signals(n, N, M_int, a0,...
                fN, gN, rcell, K, Con, Coff);
            t_mat(n1+1, n2+1, :, :) = ptsum;
        end
    end
    if max(t_mat(:))>1
       warning('Error, there is a probability > 1!');
    end
    save(fname)
end

%%
%{
n = [5 5];
t_mat(t_mat<10^(-50)) = 0; % getting rid of small values

h=figure(1);
p = (0:N)/N;
imagesc(p, p, squeeze(t_mat(n(1)+1, n(2)+1, :,:)) )
c = colorbar;
set(gca, 'Ydir', 'normal', 'FontSize', 20)
xlabel('$$p^{(1)}_{t+1}$$', 'FontSize', 24)
ylabel('$$p^{(2)}_{t+1}$$', 'FontSize', 24)
c.Label.String = 'Probability';

qsave=0;
if qsave
    folder_fig = 'H:\My Documents\Multicellular automaton\figures\main\transition_matrix';
    fname_fig = fullfile(folder_fig, strcat(fname_str, '_map'));
    save_figure(h, 10, 8, fname_fig, '.pdf')
end
%}