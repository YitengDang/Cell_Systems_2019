% Time evolution of a system without noise and with visualization
% Finite Hill coefficient
% all: also plot mean-field and uniform lattice solutions
close all
clear variables
warning off
%% Parameters
% Lattice parameters
gridsize = 15;
N = gridsize^2;
a0 = 6;
Rcell = 0.2*a0;

% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient
noise = 0; %10^(-2); %10^(0.25);
prec = 8;

% Initial configuration
initialID = 'uniform'; %'fixedp_unif'; %'uniform' 'montecarlo' 'allON' 'singleOFF' 'fixedp'
% If fixedp, assume initial states ~ Normal(p,sigma)
p = 0.5;
sigma = 0.1;
iniON = round(p*N);

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

% Find uniform lattice fixed points
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
%hfunc = @update_function_uniform;
hfunc2 = @(x) update_function_uniform(x, hill, Con, K, fN) - x; % update function(x) = x
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end

%% Calculate minimum noise
%Kmin = ((1-fp(2))/fp(2))^(1/hill)*(1+fN).*((Con-1)*[fp(1) fp(3)]+1);
%alpha_min = abs(K-Kmin)/sqrt(N);

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

% video for saving
%frames = struct('cdata',[],'colormap',[]);
% initialize figure
hin = figure(1);
t = 0;
H = [];
xmean = [];
xstd = [];
%I = [];
%cells_hist = {};

cells = cells_hist{1};
update_cell_figure_continuum(hin, pos, a0, cells, cell_type, t);
%frames(1) = getframe(gcf);

xmean(end+1) = mean(cells);
xstd(end+1) = std(cells);
%I(end+1) = moranI(cells, a0*dist);
[cells_out, changed, H(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);

while changed && t<tmax
    disp(t);
    t = t+1;
    %k = waitforbuttonpress;
    %pause(0.01);
    update_cell_figure_continuum(hin, pos, a0, cells_out, cell_type, t);
    %frames(end+1) = getframe(gcf);
    cells_hist{end+1} = cells_out;
    xmean(end+1) = mean(cells_out);
    xstd(end+1) = std(cells);
    %I(end+1) = moranI(cells, a0*dist);
    cells = cells_out;
    [cells_out, changed, H(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec);
end
t_out = t;
%}
%% Alternative: load trajectory, run mean-field and uniform lattice with initial state of loaded trajectory

% Load previous simulation
fname = fullfile('N:\tnw\BN\HY\Shared\Yiteng\one_signal_finite_Hill\sample_trajectories',...
    'N121_a0_5p00_K8_Con16_hill2p00_t73_xmeanf_0p66_prec_8_uniform-v3_full_synchronization');
load(fname, 'cells_hist');%% Run simulation

%% Run mean-field simulation

xmean_mf = [];
xstd_mf = [];
%I_mf = [];

t = 0;
cells = cells_hist{1};
hin = figure;
update_cell_figure_continuum(hin, pos, a0, cells, cell_type, t);

xmean_mf(end+1) = mean(cells);
xstd_mf(end+1) = std(cells);
%I_mf(end+1) = moranI(cells, a0*dist);
[cells_out, changed, H] = update_cells_continuum_mf(cells, Con, K, hill, fN, prec);
while changed && t<tmax
    disp(t);
    t = t+1;
    %k = waitforbuttonpress;
    pause(0.01);
    update_cell_figure_continuum(hin, pos, a0, cells_out, cell_type, t);
    %frames(end+1) = getframe(gcf);
    cells_hist{end+1} = cells_out;
    xmean_mf(end+1) = mean(cells_out);
    xstd_mf(end+1) = std(cells);
    %I_mf(end+1) = moranI(cells, a0*dist);
    cells = cells_out;
    [cells_out, changed, H] = update_cells_continuum_mf(cells, Con, K, hill, fN, prec);
end
t_out_mf = t;


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
hold on
mksize = 8;
plot(0:t_out, xmean, 'bo-', 'LineWidth', 1.5, 'MarkerSize', mksize);
plot(0:t_out_mf, xmean_mf, 'b*-', 'LineWidth', 1.5, 'MarkerSize', mksize);
plot(0:tmax, x_uniform, 'b.-', 'LineWidth', 1.5, 'MarkerSize', mksize);

set(gca, 'FontSize', 20);
xlabel('Time');
ylabel('$$\langle X_i \rangle$$'); % correct label
xlim([0 t]);

% Plot sigmaX(t)
h=figure;
hold on
mksize = 8;
plot(0:t_out, xstd, 'r-', 'MarkerSize', mksize);
plot(0:t_out_mf, xstd_mf, 'r--', 'MarkerSize', mksize);
set(gca, 'FontSize', 20);
xlabel('Time');
ylabel('$$\sigma(X_i)$$'); % correct label
xlim([0 t]);

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