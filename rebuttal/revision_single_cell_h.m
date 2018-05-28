% Obtain distribution of hi = -Xi(Yi-K)
clear all
close all

% Parameters
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 16;
Con = 8;

% initial conditions
p0 = 0.9;
iniON = round(p0*N);

% calculate fN
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
% get matrix of interaction strengths
%r2 = dist(dist>0);
%frij = reshape(sinh(Rcell)*exp(Rcell-r2)./r2, [N-1, N]);

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

%%
% Run simulation
% initialize ON cells
cells = zeros(N,1);
cells(randperm(N,iniON)) = 1;

hin = figure(1);
t = 0;
I = [];
Non = [];
h = [];
cells_hist = {};
hlist_all = {};
update_cell_figure(hin, pos, a0, cells, cell_type, t);
% save vars and update cells
cells_hist{end+1} = cells;
Non(end+1) = sum(cells);
I(end+1) = moranI(cells, a0*dist);
[cells_out, changed, h(end+1), hlist] = update_cells(cells, dist, Con, K, a0, Rcell);
hlist_all{end+1} = hlist; 
while changed
    pause(0.5);
    t = t+1;
    %k = waitforbuttonpress;
    update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
    cells_hist{end+1} = cells_out;
    I(end+1) = moranI(cells_out, a0*dist);
    Non(end+1) = sum(cells_out);
    cells = cells_out;
    [cells_out, changed, h(end+1), hlist] = update_cells(cells, dist, Con, K, a0, Rcell);
    hlist_all{end+1} = hlist;
end

%% Save trajectory
qsave = 0;
if qsave
    fname_str = sprintf('dynamics_N%d_n%d_neq_%d_a0_%d_K_%d_Son_%d_t_%d', ...
        N, iniON, Non(end), 10*a0, K, Con, t);
    fname = fullfile(pwd, 'rebuttal', 'single_cell_h', fname_str);
    save(fname);
end
%% Plot histogram
for i=1:t+1
    hlist = hlist_all{i};
    h1=figure();
    hbins = -5:0.5:2;
    histogram(hlist, hbins, 'Normalization', 'pdf');
    xlabel('H_i');
    ylabel('P(H_i)');
    set(gca,'FontSize', 24);

    qsave = 0;
    if qsave
        fname_str = sprintf('h_hist_N%d_n%d_neq_%d_a0_%d_K_%d_Son_%d_time_%d', ...
            N, iniON, Non(end), 10*a0, K, Con, i);
        fname = fullfile(pwd, 'rebuttal', 'single_cell_h', fname_str);
        save_figure_pdf(h1, 10, 8, fname);
    end
end
%% Plot trajectories
h_traj = zeros(N,t+1);
for i=1:t+1
    h_traj(:, i) = hlist_all{i};
end

h2=figure();
hold on
plot(0:t, h_traj);
plot(0:t, mean(h_traj, 1), 'k', 'LineWidth', 4);
xlabel('t');
ylabel('H_i');
set(gca,'FontSize', 24);

qsave = 0;
if qsave
    fname_str = sprintf('h_traj_N%d_n%d_neq_%d_a0_%d_K_%d_Son_%d_t_%d', ...
        N, iniON, Non(end), 10*a0, K, Con, t);
    fname = fullfile(pwd, 'rebuttal', 'single_cell_h', fname_str);
    save_figure_pdf(h2, 10, 8, fname);
end

%{
% Hamiltonian
E = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
piv = (0:N)/N;
Iv = -0.05:0.05:1;
[p_i,Imesh] = meshgrid(piv, Iv);
h1 = figure();
contourf(piv, Iv, E(p_i, Imesh),'LineStyle', 'none')
colormap('summer')
c = colorbar;
ylim([-0.05 1]);

% Imax 
figure(h1);
hold on
fa0 = sinh(Rcell)*exp(Rcell-a0)/a0;
Imax = @(p) (6 - 4./sqrt(p*N) - 6./(p*N) - 6*p)*fa0./(1-p)/fN;
Imax2 = @(p) (6*p - 4./sqrt((1-p)*N) - 6./((1-p)*N))*fa0./p/fN;
%plot(piv, Imax(piv), 'b')
plot(piv, max(Imax(piv), Imax2(piv)), '--b', 'Linewidth', 1.5)
%}