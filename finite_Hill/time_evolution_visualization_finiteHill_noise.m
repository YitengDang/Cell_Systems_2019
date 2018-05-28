% Time evolution of a system without noise and with visualization
% Finite Hill coefficient
close all
clear all
warning off
%%
% Lattice parameters
gridsize = 15;
N = gridsize^2;
a0 = 5.2;
Rcell = 0.2*a0;

% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient
prec = 5; %precision
noise = 0; %10^(1);

% Initial configuration
initialID = 'montecarlo'; %'uniform' 'allON' 'singleOFF' 'fixedp'
% If fixedp, assume initial states ~ Normal(p,sigma)
p = 0.5;
sigma = 0.1; 
iniON = round(p*N);
qsave = 0;

% simulation parameters
tmax = 100;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist = round(dist, 5);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% Find single cell fixed points
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
hfunc = @update_function_uniform;
hfunc2 = @(x) hfunc(x, hill, Con, K, fN);
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end
%% Plot phase diagram 
%{
% For single cell in uniform lattice
% Bistability example: K=8, Con=16, a0=1.5, hill=2
set(0,'defaulttextinterpreter', 'latex');
h=figure();
hold on
S = linspace(1+fN,Con*(1+fN),200);
fS = 1 + (Con-1)*S.^hill./(K^hill+S.^hill);
plot(linspace(0, 1, 200), (S-1)/(Con-1), 'b-', 'LineWidth', 1.5);
plot(linspace(0, 1, 200), (fS-1)/(Con-1), 'r-', 'LineWidth', 1.5);
%xlabel('Sensed concentration $$Y_i$$', 'FontSize', 20);
%ylabel('Secreted concentration $$f(Y_i)$$', 'FontSize', 20);
xlabel('$$X_i(t)$$');
ylabel('$$X_i(t+1)$$');
legend({'X_i', 'f(X_i)'}, 'FontSize', 20, 'Location', 'nw');
set(gca, 'FontSize', 24);
xlim([0 1]);
ylim([0 1]);

%Save fig
fname = strrep(sprintf('single_cell_dynamics_Con%d_K%d_hill%.2f_a0_%.2f', Con, K, hill, a0), '.','p');
out_file = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'figures', fname);
%out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_transition', fname);
%save_figure_pdf(h, 10, 8, out_file)
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
elseif strcmp(initialID, 'montecarlo')
    cells_hist{1} = init_random_cells_montecarlo(N, p);
end
%cells = [0.62*ones(floor(N/2),1); zeros(ceil(N/2),1)]; % extreme case: separation of ON/OFF cells

% Load previous simulation
%fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
%    'N225_a01p5_Con18p00_K14p00_hill2p00_fixedp_pini0p27_sigma0p10-v1.mat');
%load(fname, 'cells_hist');%% Run simulation
cells = cells_hist{1};
% video for saving
frames = struct('cdata',[],'colormap',[]);

% initialize figure
hin = figure(1);
t = 0;
H = [];
xmean = [];
xstd = [];
xskew = []; %skewness
xkurt = []; %kurtosis
dY_mean = [];

I = [];
Theta = [];
%cells_hist = {};
update_cell_figure_continuum(hin, pos, a0, cells, cell_type, t);
%update_cell_figure_continuum_binary(hin, pos, a0, cells, cell_type, t, fp(2));
frames(1) = getframe(gcf);
xmean(end+1) = mean(cells);
xstd(end+1) = std(cells);
%xskew(end+1) = skewness(cells);
%xkurt(end+1) = kurtosis(cells);
[I(end+1), Theta(end+1)] = moranI(cells, a0*dist);
[cells_out, changed, dY_mean(end+1)] = update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, prec);

while changed && t<tmax
    t = t+1;
    %k = waitforbuttonpress;
    %pause(0.1);
    update_cell_figure_continuum(hin, pos, a0, cells_out, cell_type, t);
    %update_cell_figure_continuum_binary(hin, pos, a0, cells_out, cell_type, t, fp(2));
    frames(end+1) = getframe(gcf);
    cells_hist{end+1} = cells_out;
    xmean(end+1) = mean(cells_out);
    xstd(end+1) = std(cells_out);
    %xskew(end+1) = skewness(cells);
    %xkurt(end+1) = kurtosis(cells);
    [I(end+1), Theta(end+1)] = moranI(cells, a0*dist);
    cells = cells_out;
    [cells_out, changed, dY_mean(end+1)] = update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, prec);
end

% Check stability of final state
[stable, evals] = stability_fun(cells_out, a0, Rcell, K, Con, hill, dist);

%% Calculate uniform lattice solution
%{
Xm_ini = mean(cells_hist{1}); % initial Xm
Xm = zeros(1, t+1);
Xm(1) = Xm_ini;
for i=1:t
    Y = (1+fN)*(Xm(i)*(Con-1) + 1);
    Xm(i+1) = Y^hill/(K^hill + Y^hill);
end
%}
%% Plot muX(t)=<Xi>(t) and sigmaX(t)
%
h=figure(2);
hold on
mksize = 8;
%plot(0:t, Xm, 'gd', 'LineWidth', 1.5, 'MarkerSize', mksize); % mean field
%errorbar(0:t, xmean, xstd, 'bo-');
plot(0:t, xmean, 'bo', 'MarkerSize', mksize);
plot(0:t, xstd, 'ro', 'MarkerSize', mksize);
%plot(0:t, xskew, 'yo-');
%plot(0:t, xkurt, 'co-');

set(gca, 'FontSize', 20);
xlabel('time');
ylabel('$$\langle X_i \rangle$$'); % correct label
%xlim([0 40]);
%legend({'\langle X_i \rangle', '\sigma_X'});
%legend({'\mu_X (Standard)', '\sigma_X (Standard)', '\mu_X (Mean field)', '\sigma_X (Mean field)'});
%legend({'simulation', 'calculation'});
%legend({'\mu_X', '\sigma_X', 'Skew[X]', 'Kurt[X]',...
%    '\mu_X (Unif)'}, 'Location', 'e');
%legend({'\mu_X', '\sigma_X',...
%    '\mu_X (Unif)'}, 'Location', 'e');
xlim([0 t]);

% ----(Optional) include Master Equation result----------------------------
%{
tmax = min(t, numel(MExmean)-1);
plot(0:tmax, MExmean(1:tmax+1), 'x', 'Color', [0 0 1],...
    'LineWidth', 1.5, 'MarkerSize', mksize);
plot(0:tmax, sqrt(MExvar(1:tmax+1)), 'x', 'Color', [1 0 0],...
    'LineWidth', 1.5, 'MarkerSize', mksize);
%}
% -------------------------------------------------------------------------

%Save fig
qsave=0;
if qsave
    if strcmp(initialID,'uniform')
        fname_str = strrep(sprintf('trajectory_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s',...
            N, a0, Con, K, hill, initialID), '.','p');
    else
        fname_str = strrep(sprintf('trajectory_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
            N, a0, Con, K, hill, initialID, p, sigma), '.','p');
    end
    i = 1;
    fname = fullfile(pwd, 'figures', 'time_evolution_moments_ME',...
        strcat(fname_str,'-v',int2str(i)));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd, 'figures', 'time_evolution_moments_ME',...
            strcat(fname_str,'-v',int2str(i),'.pdf'));
    end
    save_figure_pdf(h, 10, 8, fname)
    save_figure_eps(h, 10, 8, fname)
end
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
%% Correlation
%{
C = Theta - (2*xmean-1).^2;
figure();
plot(0:t, C);
%}
%% Plot dY_mean
figure();
plot(0:t, dY_mean);
%% Plot in p, C space
%{
C = Theta - (2*xmean-1).^2;

figure(6);
hold on
% plot trajectories
p1=plot(xmean, C, 'b');
%plot starting and endpoints
plot(xmean(1), C(1),'bo');
plot(xmean(end), C(end),'bx');
%}
%% Plot hamiltonian
%{
% calculate hamiltonian for uniform lattice
%hm = -Xm.*((1+fN)*(Xm*(Con-1)+1) - K);
%hm = (Xm(2:end) - Xm(1:end-1)).^2;
hvals = E(xmean, I);

h=figure(4);
hold on
tmax = t;
%plot(0:tmax, H(1:tmax+1)/N, 'b','LineWidth', 1.5); %h = H/N
plot(0:tmax, hvals(1:tmax+1), 'b','LineWidth', 1.5); %h = H/N
%plot(0:tmax, HMF(1:tmax+1)/N, 'r','LineWidth', 1.5); %h = H/N
%plot(0:tmax, hm(1:tmax+1), 'g','LineWidth', 1.5); %h = H/N
set(gca, 'FontSize', 20);
xlabel('time');
ylabel('$$h(t)$$');
legend({'Standard', 'Mean field', 'Uniform lattice'});

%Save fig
if qsave
    fname_str = strrep(sprintf('hamiltonian_new_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
        N, a0, Con, K, hill, initialID, p, sigma), '.','p');
    i = 1;
    fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'fig_bistable_sim',...
        strcat(fname_str,'-v',int2str(i),'.pdf'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'fig_bistable_sim',...
            strcat(fname_str,'-v',int2str(i),'.pdf'));
    end
    save_figure_pdf(h, 10, 8, fname)
end
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
plot(0:t, Xm, 'g--', 'LineWidth', 2, 'DisplayName', '$$X_m(t)$$');
plot(0:t, mean(all_trajectories), 'y--', 'LineWidth', 2, 'DisplayName','<$$X_i$$(t)>');
plot(0:t, repmat(fp, 1, t+1), 'r--', 'LineWidth', 2);
%legend('X_i(t)', 'X_m(t)', '<X_i(t)>', 'X*');
%legend('show')
set(gca, 'FontSize', 24);
xlabel('$$t$$');
ylabel('$$X_i(t)$$');
set(gcf, 'Units', 'Inches', 'Position', [0 0 10 8]);
%% Plot distribution of Xi at fixed time
h8 = figure(8);
hold on
time = 56;
nbins = numel(binc);
binedges = [binc-1/nbins/2 1];
histogram(cells_hist{time}, binedges, 'Normalization', 'probability');
plot(binc, MEprob(:,11), 'r-');
xlim([0 1]);
%% Plot Xi values of final nonuniform state on a line
h9 = figure(9);
this_cells = cells_out;
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