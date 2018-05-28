% Time evolution of a system without noise and with visualization
% Finite Hill coefficient, mean field approach
close all
clear all
warning off
set(0, 'defaulttextinterpreter', 'latex');

% Lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 6;
Rcell = 0.2*a0;
% circuit parameters
K = 10;
Con = 20;
hill = 2; % Hill coefficient

% Initial configuration
initialID = 'montecarlo';  %see below
p = 0.5;
sigma = 0.27;
prec = 12; % precision
iniON = round(p*N);

% simulation parameters
qsave = 0;
%dir = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'fig_bistable_sim');
dir = fullfile('H:\My Documents\Multicellular automaton\latex', '2_finite_Hill_nonuniform', 'fig_trajectories');
tmax = 10000;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist = round(dist, 5);

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% Uniform lattice stationary solutions
syms x
S = @(x) (Con-1)*x + 1;
f = @(x) (S(x)*(1+fN))^hill/(K^hill + (S(x)*(1+fN))^hill) - x;
xMF = [fzero(f, [0 0.1]) fzero(f, [0.1 0.2]) fzero(f, [0.2 1])]; % manually set bounds

%% Plot phase diagram
% Bistability example: K=8, Con=16, a0=1.5, hill=2
set(0,'defaulttextinterpreter', 'latex');
h = figure(1);
hold on
X = linspace(0,1,100);
fX = ((Con-1)*X+1).^hill./(((Con-1)*X+1).^hill + K.^hill);
plot(X, fX, 'b-', 'LineWidth', 2);
plot(X, X, 'r-', 'LineWidth', 2);
xlabel('$$X_t$$', 'FontSize', 28);
ylabel('$$X_{t+1}$$', 'FontSize', 28);
set(gca,'FontSize', 24);
%legend({'f(Y_i)', 'Y_i'}, 'FontSize', 20, 'Location', 'nw');

%Save fig
qsave = 1;
if qsave
    basedir = 'H:\My Documents\Multicellular automaton';
    fname = strrep(sprintf('single_cell_dynamics_Con%d_K%d_hill%.2f', Con, K, hill), '.','p');
    out_file = fullfile(basedir, 'latex', '1_finite_Hill', 'figures', fname);
    save_figure_pdf(h, 10,8, out_file)
end
%

% Mean field / uniform lattice phase diagram
%{
Xi = linspace(0, 1, 100);
Yi = (1+fN)*((Con-1)*Xi + 1);
fYi = Yi.^hill./(K^hill+Yi.^hill);
h=figure();
hold on
plot(Xi, Xi, 'b-', 'LineWidth', 2);
plot(Xi, fYi, 'r-', 'LineWidth', 2);
xlabel('$$X_i$$');
ylabel('$$f(X_i)$$');
set(gca,'FontSize', 24);
%}
% Save figure
if qsave
    fname = strrep(sprintf('stability_diagram_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s',...
        N, a0, Con, K, hill, initialID), '.','p');
    out_file = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'fig_bistable_sim', fname);
    save_figure_pdf(h, 10, 8, out_file)
end
%}
%% initialize cells
cells_hist = {};
if strcmp(initialID, 'uniform')
    cells_hist{1} = rand(N, 1); %uniformly distributed
elseif strcmp(initialID, 'singleOFF')
    cells_hist{1} = ones(N, 1);
    x = randi(N);
    cells_hist{1}([x x+1 x+gridsize x+gridsize+1]) = 0;
elseif strcmp(initialID, 'fixedp')
    cells_hist{1} = p + sigma*randn(N, 1);
    cells_hist{1}(cells_hist{1} < 0) = 0;
    cells_hist{1}(cells_hist{1} > 1) = 1;
elseif strcmp(initialID, 'fixedp_unif')
    cells_hist{1} = p*ones(N, 1);
elseif strcmp(initialID, 'binary')
    cells_temp = zeros(N, 1);
    cells_temp(randperm(N, iniON)) = 1;
    cells_hist{1} = cells_temp;
elseif strcmp(initialID, 'binaryrand')
    cells_hist{1} = randi(2, N, 1) - 1;
elseif strcmp(initialID, 'montecarlo')
    cells_hist{1} = init_random_cells_montecarlo(N, p);
elseif strcmp(initialID, 'in_out_binary')
    cells_hist{1} = init_random_cells_binary_constraint(N, iniON, fp);
end
%cells = [0.62*ones(floor(N/2),1); zeros(ceil(N/2),1)]; % extreme case: separation of ON/OFF cells

% Load previous simulation
%fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
%    'N225_a01p5_Con18p00_K14p00_hill2p00_fixedp_pini0p27_sigma0p10-v1.mat');
%load(fname, 'cells_hist');
%% Run simulation
cells = cells_hist{1};
% video for saving
frames = struct('cdata',[],'colormap',[]);

% initialize figure
hin = figure(1);
t = 0;
H = [];
xmean = [];
xstd = [];
I = [];
%cells_hist = {};

%update_cell_figure_continuum(hin, pos, a0, cells, cell_type, t);
frames(1) = getframe(gcf);
xmean(end+1) = mean(cells);
xstd(end+1) = std(cells);
I(end+1) = moranI(cells, a0*dist);
[cells_out, changed, H(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);

while changed && t<tmax
    disp(t);
    t = t+1;
    %k = waitforbuttonpress;
    %pause(0.1);
    %update_cell_figure_continuum(hin, pos, a0, cells_out, cell_type, t);
    frames(end+1) = getframe(gcf);
    cells_hist{end+1} = cells_out;
    xmean(end+1) = mean(cells_out);
    xstd(end+1) = std(cells);
    I(end+1) = moranI(cells, a0*dist);
    cells = cells_out;
    [cells_out, changed, H(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
end
%% Plot final configuration
h = figure(1);
update_cell_figure_continuum(h, pos, a0, cells_hist{end}, cell_type, t);

%% Plot muX(t)=<Xi>(t) with MF steady states
h=figure(3);
% create fill between curves
curve1 = xmean+xstd;
curve2 = xmean-xstd;
x2 = [0:t, fliplr(0:t)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.6 0.7 1]);
hold on

% plot mean and std 
plot(0:t, xmean, 'b-', 'Linewidth', 3);
plot(0:t, curve1, 'b-', 'Linewidth', 1.5);
plot(0:t, curve2, 'b-', 'Linewidth', 1.5);

% Plot together with mean field steady state
plot([0 t], [xMF(1) xMF(1)], 'r--', 'Linewidth', 2)
plot([0 t], [xMF(2) xMF(2)], 'r--', 'Linewidth', 2)
plot([0 t], [xMF(3) xMF(3)], 'r--', 'Linewidth', 2)

set(gca, 'FontSize', 20);
xlabel('time');
ylabel('$$\langle X_i \rangle$$'); % correct label
xlim([0 t]);
ylim([-0.2 1]);
%legend({'\langle X_i \rangle', '\sigma_X'});
%legend({'\mu_X (Standard)', '\sigma_X (Standard)', '\mu_X (Mean field)', '\sigma_X (Mean field)'});

%Save fig
qsave = 1;
if qsave
    %fname_str = strrep(sprintf('trajectory_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f_ebar_MF_sol_style2',...
    %    N, a0, Con, K, hill, initialID, p, sigma), '.','p');
    fname_str = strrep(sprintf('trajectory_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_prec%d_w_error_unif_sol_style2',...
        N, a0, Con, K, hill, initialID, p, prec), '.','p');
    i = 1;
    fname = fullfile(dir,...
        strcat(fname_str,'-v',int2str(i),'.pdf'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(dir, strcat(fname_str,'-v',int2str(i),'.pdf'));
    end
    save_figure_pdf(h, 10, 8, fname)
end
%% Plot fraction of ON cells
fracON = zeros(t+1, 1);
for i=1:t+1
    cells = cells_hist{i};
    fracON(i) = sum(cells>xMF(2))/N;
end

h=figure(2);
plot(0:t, fracON, 'LineWidth', 1.5);
xlim([0 t]);
ylim([0 1]); 
xlabel('time');
ylabel('fraction ON cells');
set(gca, 'FontSize', 20);

qsave = 1;
if qsave
    %fname_str = strrep(sprintf('trajectory_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f_ebar_MF_sol_style2',...
    %    N, a0, Con, K, hill, initialID, p, sigma), '.','p');
    fname_str = strrep(sprintf('frac_ON_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_prec%d_w_error_unif_sol_style2',...
        N, a0, Con, K, hill, initialID, p, prec), '.','p');
    i = 1;
    fname = fullfile(dir,...
        strcat(fname_str,'-v',int2str(i),'.pdf'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(dir, strcat(fname_str,'-v',int2str(i),'.pdf'));
    end
    save_figure_pdf(h, 5, 4, fname)
end
%% Mean field simulation
%{
% mean field with dynamics
HMF = [];
xmeanMF = [];
xstdMF = [];
IMF = [];
%cells_histMF = {};
tMF = 0;
cellsMF = cells_hist{1};

% video for saving
framesMF = struct('cdata',[],'colormap',[]);

% initialize figure
hin = figure(2);
update_cell_figure_continuum(hin, pos, a0, cellsMF, cell_type, tMF);
framesMF(1) = getframe(gcf);

%cells_histMF{end+1} = cellsMF;
xmeanMF(end+1) = mean(cellsMF);
xstdMF(end+1) = std(cellsMF);
IMF(end+1) = moranI(cellsMF, a0*dist);
[cells_outMF, changed, HMF(end+1)] = update_cells_continuum_mf(cellsMF, Con, K, hill, fN);

while changed && tMF<tmax
    tMF = tMF+1;
    %k = waitforbuttonpress;
    %pause(0.2);
    update_cell_figure_continuum(hin, pos, a0, cells_outMF, cell_type, tMF);
    framesMF(end+1) = getframe(gcf);

    %cells_histMF{end+1} = cells_outMF;
    xmeanMF(end+1) = mean(cells_outMF);
    xstdMF(end+1) = std(cells_outMF);
    IMF(end+1) = moranI(cells, a0*dist);
    cellsMF = cells_outMF;
    [cells_outMF, changed, HMF(end+1)] = update_cells_continuum_mf(cellsMF, Con, K, hill, fN);
end
%% Uniform lattice simulation
% Calculate trajectory with all cells Xm
Xm_ini = mean(cells_hist{1}); % initial Xm
tsteps = max(t, tMF)+1;
Xm = zeros(1, tsteps+1);
Xm(1) = Xm_ini;
for i=1:tsteps
    Y = (1+fN)*(Xm(i)*(Con-1) + 1);
    Xm(i+1) = Y^hill/(K^hill + Y^hill);
end

%% Compare dynamics (sim, MF and unif)
h=figure(4);
hold on
%errorbar(0:t, xmean, xstd, 'bo-');
plot(0:t, xmean, 'bx-');
plot(0:t, xstd, 'ro-');

% Plot together with mean field dynamics
plot(0:tMF, xmeanMF, 'bx--');
plot(0:tMF, xstdMF, 'rx--');

% Plot uniform lattice 
h;
hold on
plot(0:tsteps-1, Xm(1:end-1), 'g-', 'LineWidth', 1.5);
%legend({'simulation', 'calculation'});
legend({'\mu_X', '\sigma_X',...
    '\mu_X (MF)', '\sigma_X (MF)',...
    '\mu_X (Unif)'}, 'Location', 'e');

set(gca, 'FontSize', 20);
xlabel('time');
ylabel('$$\langle X_i \rangle$$'); % correct label
xlim([0 40]);
%legend({'\langle X_i \rangle', '\sigma_X'});
%legend({'\mu_X (Standard)', '\sigma_X (Standard)', '\mu_X (Mean field)', '\sigma_X (Mean field)'});

%Save fig
%qsave=1;
if qsave
    %fname_str = strrep(sprintf('trajectory_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
    %    N, a0, Con, K, hill, initialID, p, sigma), '.','p');
    fname_str = strrep(sprintf('trajectory_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f',...
        N, a0, Con, K, hill, initialID, p), '.','p');

    i = 1;
    fname = fullfile(dir,...
        strcat(fname_str,'-v',int2str(i),'.pdf'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd, 'latex', '17-08-14_finite_Hill_nonuniform', 'fig_trajectories',...
            strcat(fname_str,'-v',int2str(i),'.pdf'));
    end
    save_figure_pdf(h, 10, 8, fname)
end
%}
%% Plot hamiltonian
%{
% calculate hamiltonian for uniform lattice
%hm = -Xm.*((1+fN)*(Xm*(Con-1)+1) - K);
hm = (Xm(2:end) - Xm(1:end-1)).^2;

h=figure(5);
hold on
tmax = 10;
plot(0:tmax, H(1:tmax+1)/N, 'b','LineWidth', 1.5); %h = H/N
plot(0:tmax, HMF(1:tmax+1)/N, 'r','LineWidth', 1.5); %h = H/N
plot(0:tmax, hm(1:tmax+1), 'g','LineWidth', 1.5); %h = H/N
set(gca, 'FontSize', 20);
xlabel('time');
ylabel('$$h(t)$$');
legend({'Standard', 'Mean field', 'Uniform lattice'});

%Save fig
if qsave
    fname_str = strrep(sprintf('hamiltonian_new_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
        N, a0, Con, K, hill, initialID, p, sigma), '.','p');
    i = 1;
    fname = fullfile(dir,...
        strcat(fname_str,'-v',int2str(i),'.pdf'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(dir,...
            strcat(fname_str,'-v',int2str(i),'.pdf'));
    end
    save_figure_pdf(h, 10, 8, fname)
end
%% Plot (p,I) trajectory
figure(6);
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
p2=plot(xmeanMF, IMF, 'r');
%plot endpoints
plot(xmean(end), I(end),'bx');
plot(xmeanMF(end), IMF(end),'rx');

% add legend
legend([p1 p2], {'normal', 'mean field'});
%}
%% Save whole trajectory if interesting
%{
interesting = 0;
if interesting
    fname_str = strrep(sprintf('N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
        N, a0, Con, K, hill, initialID, p, sigma), '.','p');
    i = 1;
    fname = fullfile(dir,...
        strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(dir,...
            strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    save(fname)
end
%}
%%
% Check movies
%figure(7);
%movie(framesMF, 1, 2)

%%
% Export movies
%{
% normal
fname_str = strrep(sprintf('N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
    N, a0, Con, K, hill, initialID, p, sigma), '.','p');
i = 1;
fname = fullfile(dir,...
    strcat(fname_str,'-v',int2str(i),'.avi'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(dir,...
        strcat(fname_str,'-v',int2str(i),'.avi'));
end

myVideo = VideoWriter(fname, 'Uncompressed AVI');
myVideo.FrameRate = 10;  % Default 30
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);
%%
%{
% mean field 
fname_str = strrep(sprintf('N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f_MF',...
    N, a0, Con, K, hill, initialID, p, sigma), '.','p');
i = 1;
fname = fullfile(pwd, dir,...
    strcat(fname_str,'-v',int2str(i),'.avi'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(pwd, dir,...
        strcat(fname_str,'-v',int2str(i),'.avi'));
end

myVideo = VideoWriter(fname, 'Uncompressed AVI');
myVideo.FrameRate = 15;  % Default 30
open(myVideo);
writeVideo(myVideo, framesMF);
close(myVideo);
%}
%}