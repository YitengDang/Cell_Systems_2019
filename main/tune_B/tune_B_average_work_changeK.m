% Get average work over trajectories starting with random configuration by
% using transition matrix for p to calculate densities at later times.
% Works for protocols that change K only, because only in this case does the
% work depend on p only (and not on I).

clear variables
close all
warning off
clc
%%
% 1. Get transition matrix
% Refer to transition_matrix_func <-- transition_prob_func
% Parameters
gridsize = 11;
N = gridsize^2;

a0 = 1.5;
Rcell = 0.2*a0;

% Save?
qsave = [1 1]; %(1) figures, (2) data

% calculate fN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
N = numel(dist_vec);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% Circuit parameters
K0 = 12;
Con = 2*K0/(1+fN) - 1;
K = 12:-1:2;

%---Plot transition matrix---
% Calls code that imports or calculates transition matrix
[t_mat] = transition_matrix_func(gridsize, Con, K0, a0, Rcell);
% NB t_mat needs to be transposed
figure();
p = (0:N)/N;
imagesc(p, p, t_mat')
c = colorbar;
set(gca, 'Ydir', 'normal')
xlabel('p_{t}', 'FontSize', 24)
ylabel('p_{t+1}', 'FontSize', 24)

%%
% 2. Calculate (p,I) distributions as functions of time
tmax = length(K)-1;
wpt = zeros(N+1, tmax); % w(p; t) prob. distr. for p as a function of t

% Initial prob distribution for p = binomial distribution
fileid = sprintf('nchoosek_n%d', N);
fname = fullfile(pwd,'data','nchoosek', strcat(fileid, '.mat'));
if exist(fname)==2
    disp('exists!');
    load(fname);
else
    disp('Doesnt exist!');
    binom = zeros(1,N+1);
    for j=1:N+1  
        disp(j);
        binom(j) = nchoosek(N,j-1);        
    end
    save(fname, 'binom');
end

% calculate evolution of prob distribution
wp = binom'/2^N; % binomial distribution for p
wpt(:, 1) = wp; % initial distribution
for t=1:tmax
    disp('t='); disp(t);
    %calculate/load transition matrix for given value of K
    [t_mat] = transition_matrix_func(gridsize, Con, K(t+1), a0, Rcell); 
    wpt(:, t+1) = (t_mat')*wpt(:, t); % calculate new distribution
    disp('Normalized?'); disp(all(sum(wpt(:,t+1)))); % check that the distribution is still normalized
end

%%
% plot w(p; t)
h=figure();
hold on
imagesc(0:tmax, p, wpt)
c = colorbar;
set(gca, 'Ydir', 'normal')
set(gca, 'FontSize', 24);
xlabel('t', 'FontSize', 24)
ylabel('p', 'FontSize', 24)
ylabel(c, 'w(p; t)', 'FontSize', 24)
%
% plot mean p in the same figure
%figure(2);
wpt_mean = sum(p'.*wpt); % <p(t)> = sum_p p w(p; t)
plot(0:tmax, wpt_mean, 'ro-', 'LineWidth', 1.5);
xlim([-0.5 tmax+0.5]);
ylim([0-0.5/N 1+0.5/N]);

if qsave(1)
    fileid = strrep(sprintf('wpt_N%d_Con_%.1f_K%.1fto%.1f_a0%.1f_Rcell%.1f',...
        N, Con, K(1), K(end), a0, Rcell), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'distribution_pI_from_tmatrix', fileid);
    save_figure_pdf(h, 10, 8, out_file);
    savefig(h, out_file);
end
%%
% 3. Compute <W>
% NB cannot compute <Delta H> because this also depends on I
dK = K(2:end) - K(1:end-1); % changes in K from procedure
W = sum(sum((2*p'-1).*dK.*wpt(:, 1:end-1))); 
disp('<W>='); disp(W)

% save data
if qsave(2)
    close all;
    fileid = strrep(sprintf('wpt_N%d_Con_%.1f_K%.1fto%.1f_a0_%.1f_Rcell%.1f',...
        N, Con, K(1), K(end), a0, Rcell), '.', 'p');
    out_file = fullfile(pwd, 'data', 'distribution_pI_from_tmatrix', fileid);
    save(out_file, 'wpt', 'wpt_mean', 'W', 'K', 'Con', 'a0', 'fN')
end