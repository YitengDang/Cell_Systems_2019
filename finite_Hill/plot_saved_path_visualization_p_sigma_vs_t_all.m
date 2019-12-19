% Time evolution of a system without noise and with visualization
% Finite Hill coefficient
% all: also plot mean-field and uniform lattice solutions
close all
clear variables
warning off
set(0, 'defaulttextinterpreter', 'tex');
%% Alternative: load trajectory, run mean-field and uniform lattice with initial state of loaded trajectory
% Load previous simulation
fname = fullfile('N:\tnw\BN\HY\Shared\Yiteng\one_signal_finite_Hill\sample_trajectories',...
    'N121_a0_5p60_K8_Con16_hill2p00_t172_xmeanf_0p69_prec_8_uniform-v1_partial_synchronization');
initialID = 'uniform';
load(fname);%% Run simulation
%%
% Lattice parameters
N = save_consts_struct.N;
gridsize = sqrt(N);
a0 = save_consts_struct.a0;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters
K = save_consts_struct.K;
Con = save_consts_struct.Con;
hill = save_consts_struct.hill; % Hill coefficient 
noise = save_consts_struct.noise; %10^(-2); %10^(0.25);
prec = save_consts_struct.prec;

% simulation parameters
tmax = 200;

% use hexagonal lattice
dist = round(distances, 5);
pos = round(positions, 5);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% Replay trajectory, get mean and std
xmean = zeros( numel(cells_hist), 1 );
xstd = zeros( numel(cells_hist), 1 );
t_out = numel(cells_hist) - 1;

hin = figure;
disp_mol = 1;
showI = 0;
h = reset_cell_figure(hin, pos, rcell);
cells = cells_hist{1};
xmean(1) = mean(cells);
xstd(1) = std(cells);
t=0;
update_figure_periodic_scatter(h, cells, t, disp_mol, showI, a0, dist);
for t=1:numel(cells_hist)-1
    pause(0.01);
    cells = cells_hist{t+1};
    xmean(t+1) = mean(cells);
    xstd(t+1) = std(cells);
    update_figure_periodic_scatter(h, cells, t, disp_mol, showI, a0, dist);
end

% continue simulation until tmax
cells_out = cells_hist{end};
while t<tmax
    pause(0.01);
    disp(t);
    t = t+1;
    
    cells = cells_out;
    [cells_out, ~, ~] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);

    cells_hist{end+1} = cells_out;
    xmean(end+1) = mean(cells_out);
    xstd(end+1) = std(cells_out);
    %I_mf(end+1) = moranI(cells, a0*dist);
    
    update_figure_periodic_scatter(h, cells, t, disp_mol, showI, a0, dist);
end
%% Run mean-field simulation
cells_hist_mf = {};
xmean_mf = [];
xstd_mf = [];
%I_mf = [];

t = 0;
cells = cells_hist{1};
cells_hist_mf{end+1} = cells;
xmean_mf(end+1) = mean(cells);
xstd_mf(end+1) = std(cells);
%I_mf(end+1) = moranI(cells, a0*dist);    
[cells_out, changed, H] = update_cells_continuum_mf(cells, Con, K, hill, fN, prec);
while t<tmax
    disp(t);
    t = t+1;
    %k = waitforbuttonpress;
    cells_hist_mf{end+1} = cells_out;
    xmean_mf(end+1) = mean(cells_out);
    xstd_mf(end+1) = std(cells_out);
    %I_mf(end+1) = moranI(cells, a0*dist);
    cells = cells_out;
    [cells_out, changed, H] = update_cells_continuum_mf(cells, Con, K, hill, fN, prec);
end
t_out_mf = t;

% play trajectory
hin = figure;
disp_mol = 1;
showI = 0;
h = reset_cell_figure(hin, pos, rcell);
cells = cells_hist_mf{1};
t=0;
update_figure_periodic_scatter(h, cells, t, disp_mol, showI, a0, dist);
for t=1:numel(cells_hist_mf)-1
    pause(0.01);
    cells = cells_hist_mf{t+1};
    update_figure_periodic_scatter(h, cells, t, disp_mol, showI, a0, dist);
end

%% Run uniform lattice simulation
x_uniform = zeros(tmax+1, 1);
x_uniform(1) = mean( cells_hist{1}, 1);
this_x = x_uniform(1);
for t=1:tmax
    this_x = update_function_uniform(this_x, hill, Con, K, fN);
    x_uniform(t+1) = this_x;
end

%% Plot muX(t)=<Xi>(t)
h=figure;
box on
hold on
mksize = 8;
lw = 1;
plot(0:tmax, xmean, 'b-', 'LineWidth', lw, 'MarkerSize', mksize);
%plot(0:tmax, xmean_mf, 'r--', 'LineWidth', lw, 'MarkerSize', mksize);
plot(0:tmax, x_uniform, 'r--', 'LineWidth', lw, 'MarkerSize', mksize);
set(gca, 'FontSize', 32);
xlabel('Time');
ylabel('\langle X_i \rangle'); % correct label
xlim([0 t]);
ylim([0 1]);
%legend({'exact', 'mean-field', 'uniform lattice'}, 'Location', 'se');
legend({'exact', 'uniform lattice'}, 'Location', 'se');

% Plot sigmaX(t)
h2=figure;
box on
hold on
mksize = 8;
plot(0:tmax, xstd, 'b-', 'LineWidth', lw,  'MarkerSize', mksize);
%plot(0:tmax, xstd_mf, 'r--', 'LineWidth', lw, 'MarkerSize', mksize);
set(gca, 'FontSize', 32);
xlabel('Time');
ylabel('Varability \sigma(X_i)'); % correct label
xlim([0 t]);
ylim([0 0.4]);
legend({'exact', 'mean-field'}, 'Location', 'ne');

%Save fig
qsave = 1;
save_folder = 'H:\My Documents\Thesis\Single gene finite Hill\Fig1';
fname_str = strrep(sprintf('trajectory_mean_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_size_12_8_v2',...
    N, a0, Con, K, hill, initialID), '.','p');
fname_str_2 = strrep(sprintf('trajectory_std_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_size_12_8_v2',...
    N, a0, Con, K, hill, initialID), '.','p');
fname = fullfile(save_folder, fname_str);
fname2 = fullfile(save_folder, fname_str_2);
save_figure(h, 12, 8, fname, '.pdf', qsave)
save_figure(h2, 12, 8, fname2, '.pdf', qsave)

%% Plot muX(t)=<Xi>(t) with std as spread
h=figure;
box on
hold on

% create fill between curves
curve1 = (xmean+xstd)';
curve2 = (xmean-xstd)';
x2 = [0:tmax, fliplr(0:tmax)];
inBetween = [curve1, fliplr(curve2)];
p=fill(x2, inBetween, [0 0 1]);
alpha(p, 0.25);

%{
curve1 = (xmean_mf+xstd_mf);
curve2 = (xmean_mf-xstd_mf);
x2 = [0:tmax, fliplr(0:tmax)];
inBetween = [curve1, fliplr(curve2)];
p=fill(x2, inBetween, [1 0 0]);
alpha(p, 0.25);
%}
% Plot mean
mksize = 8;
lw = 1;
p1 = plot(0:tmax, xmean, 'b-', 'LineWidth', lw, 'MarkerSize', mksize);
%p2 = plot(0:tmax, xmean_mf, 'r*-', 'LineWidth', lw, 'MarkerSize', mksize);
p3 = plot(0:tmax, x_uniform, 'r-', 'LineWidth', lw, 'MarkerSize', mksize);
set(gca, 'FontSize', 32);
xlabel('Time');
ylabel('Mean expression level'); % correct label
xlim([0 t]);
ylim([0 1]);
%legend([p1 p2 p3], {'exact', 'mean-field', 'uniform lattice'}, 'Location', 'se');
legend([p1 p3], {'exact', 'uniform lattice'}, 'Location', 'se');

%Save fig
qsave = 1;
save_folder = 'H:\My Documents\Thesis\Single gene finite Hill\Fig1';
fname_str = strrep(sprintf('trajectory_A2_mean_std_fill_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_size_12_8',...
    N, a0, Con, K, hill, initialID), '.','p');
fname = fullfile(save_folder, fname_str);
save_figure(h, 12, 8, fname, '.pdf', qsave)
%% Check whether trajectory fell within allowed range
%{
test_mean_all = zeros(t-len_test, 1);
test_std_all = zeros(t-len_test, 1);

for time=len_test+1:t+1
    idx = time-len_test;
    test_mean_all(idx) = std(xmean(time-len_test+1:time));
    test_std_all(idx) = std(xstd(time-len_test+1:time));   
end

figure();
hold on
plot(len_test+1:t+1, test_mean_all, 'b');
plot(len_test+1:t+1, test_std_all, 'r');
plot([len_test+1 t+1], [err_thr_mean err_thr_mean], 'b--');
plot([len_test+1 t+1], [err_thr_std err_thr_std], 'r--');
legend({'last 10 xmean', 'last 10 xstd', 'xmean threshold', 'xstd threshold'});
xlabel('time');
ylabel('value');
%}
%% Plot (p,I) trajectory
%{
figure(5);
hold on

% plot hamiltonian map
E = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
piv = (0:N)/N;
Iv = -0.05:0.05:1;
[p_i,Imesh] = meshgrid(piv, Iv);
contourf(piv, Iv, E(p_i, Imesh),'LineStyle', 'none');
colormap('summer')
c = colorbar;

% plot trajectories
p1=plot(xmean, I, 'b');
p2=plot(xmean, Theta, 'r');
%plot endpoints
plot(xmean(end), I(end),'bx');
plot(xmean(end), Theta(end),'rx');

% plot Theta(p) = (2p-1)^2 for uniform lattice
Thetav = (2*piv-1).^2;
plot(piv, Thetav, 'k--');

% add legend
legend([p1 p2], {'I', 'Theta'});
%}
%% Plot in p, C space
%{
% C: correlation, unnormalized I
C = Theta - (2*xmean-1).^2;

figure(6);
hold on
% plot trajectories
p1=plot(xmean, C, 'b');
%plot starting and endpoints
plot(xmean(1), C(1),'bo');
plot(xmean(end), C(end),'bx');
%{
figure();
plot(0:t, C);
%}
%}
%% Save whole trajectory if interesting
%{
fname_str = strrep(sprintf('N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
    N, a0, Con, K, hill, initialID, p, sigma), '.','p');
i = 1;
fname = fullfile(pwd, 'data', 'dynamics', 'finiteHill', 'single_trajectories',... 
    strcat(fname_str,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(pwd, 'data', 'dynamics', 'finiteHill', 'single_trajectories',...
        strcat(fname_str,'-v',int2str(i),'.mat'));
end
save(fname)
%}
%%
% Check movies
%figure(6);
%movie(frames, 1, 10)

%%
%{
% Export movies
% normal
fname_str = strrep(sprintf('N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
    N, a0, Con, K, hill, initialID, p, sigma), '.','p');
i = 1;
fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
    strcat(fname_str,'-v',int2str(i),'.avi'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
        strcat(fname_str,'-v',int2str(i),'.avi'));
end

myVideo = VideoWriter(fname, 'Uncompressed AVI');
myVideo.FrameRate = 3;  % Default 30
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);
%%
% mean field 
fname_str = strrep(sprintf('N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f_MF',...
    N, a0, Con, K, hill, initialID, p, sigma), '.','p');
i = 1;
fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
    strcat(fname_str,'-v',int2str(i),'.avi'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
        strcat(fname_str,'-v',int2str(i),'.avi'));
end

myVideo = VideoWriter(fname, 'Uncompressed AVI');
myVideo.FrameRate = 15;  % Default 30
open(myVideo);
writeVideo(myVideo, framesMF);
close(myVideo);
%}

%%(K-(1+fN))/(2*(1+fN)*(Con-1))

%% Plot individual trajectories together with single-cell fixed points and MF solution
%{
all_trajectories = cell2mat(cells_hist);
h7 = figure(7);
hold on
plot(0:t, all_trajectories, 'b', 'LineWidth', 0.5);
%plot(0:t, Xm, 'g--', 'LineWidth', 2, 'DisplayName', '$$X_m(t)$$');
plot(0:t, mean(all_trajectories), 'y--', 'LineWidth', 2, 'DisplayName','<$$X_i$$(t)>');
plot(0:t, repmat(fp, 1, t+1), 'r--', 'LineWidth', 2);
%legend('X_i(t)', 'X_m(t)', '<X_i(t)>', 'X*');
%legend('show')
set(gca, 'FontSize', 24);
xlabel('$$t$$');
ylabel('$$X_i(t)$$');
set(gcf, 'Units', 'Inches', 'Position', [0 0 10 8]);
%}
%% Plot distribution of Xi at fixed time
%{
h8 = figure(8);
hold on
time = 0;
%nbins = numel(binc);
%binedges = [binc-1/nbins/2 1];
histogram(cells_hist{time+1})
%histogram(cells_hist{time}, binedges, 'Normalization', 'probability');
%plot(binc, MEprob(:,11), 'r-');
xlim([0 1]);
xlabel('$$X_i$$');
ylabel('Prob');
%}
%% Plot Xi values of final nonuniform state on a line
%{
h9 = figure(9);
this_cells = cells_hist{end};
hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
             'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
             'XLim',[0 1],...               %#   set the x axis limit,
             'YLim',[0 eps],...             %#   set the y axis limit (tiny!),
             'Color','none');               %#   and don't use a background color
if ~isempty(cells(this_cells < 0.5))         
    plot(this_cells(this_cells < 0.5), 0, 'r*','MarkerSize',10);  %# Plot data set 1
end
if ~isempty(cells(this_cells > 0.5))
    plot(this_cells(this_cells > 0.5), 0, 'b*','MarkerSize',10);  %# Plot data set 2
end
xlabel('$$X_i$$');
set(gca,'FontSize', 24);    

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Sample_lattice_Xi_N%d_a0_%.2f_K%d_Con%d_hill%.2f_%s',...
        N, a0, K, Con, hill, initialID), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_transition', fname_str);
    save_figure_pdf(h8, 10, 2, out_file);
end
%}

%% Find uniform lattice fixed points
%{
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
hfunc = @update_function_uniform;
hfunc2 = @(x) hfunc(x, hill, Con, K, fN);
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end
%}
%% Plot phase diagram 
%{
% For single cell in uniform lattice
% Bistability example: K=8, Con=16, a0=1.5, hill=2
set(0,'defaulttextinterpreter', 'latex');
h1=figure(1);
hold on
S = linspace(1+fN,Con*(1+fN),200);
%fS = 1 + (Con-1)*S.^hill./(K^hill+S.^hill);
plot(linspace(0, 1, 200), linspace(0, 1, 200), 'b-', 'LineWidth', 1.5);
plot(linspace(0, 1, 200), S.^hill./(K^hill+S.^hill), 'r-', 'LineWidth', 1.5);
plot(fp, fp, 'ko', 'MarkerSize', 10);
%xlabel('Sensed concentration $$Y_i$$', 'FontSize', 20);
%ylabel('Secreted concentration $$f(Y_i)$$', 'FontSize', 20);
xlabel('$$X_i(t)$$');
ylabel('$$X_i(t+1)$$');
%legend({'X_i', 'f(X_i)'}, 'FontSize', 20, 'Location', 'nw');
legend({'X_i(t)', 'X_i(t+1)'}, 'FontSize', 20, 'Location', 'nw');

set(gca, 'FontSize', 24);
xlim([0 1]);
ylim([0 1]);

%Save fig
fname = strrep(sprintf('fixed_points_map_Con%d_K%d_hill%.2f_a0_%.2f_uniform_lat', Con, K, hill, a0), '.','p');
%out_file = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'figures', fname);
%out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_transition', fname);
out_file = fullfile('H:\My Documents\Multicellular automaton\finite_Hill\figures', fname);
save_figure_pdf(h1, 10, 8, out_file)
save_figure_eps(h1, 10, 8, out_file)
%}

%% Calculate minimum noise
%Kmin = ((1-fp(2))/fp(2))^(1/hill)*(1+fN).*((Con-1)*[fp(1) fp(3)]+1);
%alpha_min = abs(K-Kmin)/sqrt(N);
