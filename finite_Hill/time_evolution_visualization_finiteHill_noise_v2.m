% Time evolution of a system without noise and with visualization
% Finite Hill coefficient
% v2: new equilibrium condition (based on threshold on fluctuations in mean
% and std)
close all
clear variables
warning off
%%
% Lattice parameters
gridsize = 15;
N = gridsize^2;
a0 = 8;
Rcell = 0.2*a0;

% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient
noise = 10^(-2); %10^(0.25);

% Initial configuration
initialID = 'montecarlo'; %'fixedp_unif'; %'uniform' 'allON' 'singleOFF' 'fixedp'
% If fixedp, assume initial states ~ Normal(p,sigma)
p = 0.2;
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

%% Find uniform lattice fixed points
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

%% Fix the error thresholds by running simulation starting in fixed point with noise
cm = 2;
cs = 2;
[err_thr_mean, err_thr_std, ~] = fix_error_thresholds_v2(dist, a0, Con, K, fN, hill, noise, cm, cs);
fprintf('Error threshold mean %d \n ', err_thr_mean);
fprintf('Error threshold std %d \n', err_thr_std);

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

% video for saving
%frames = struct('cdata',[],'colormap',[]);

% initialize fig & vars
hin = figure(2);
t = 0;
xmean = [];
xstd = [];
corr = [];
Theta = [];
%I = [];

% initialize cells
cells = cells_hist{1};

% store initial variables
[~, theta] = moranI(cells, a0*dist);
Theta(1) = theta;
xmean(1) = mean(cells);
xstd(1) = std(cells);
corr(1) = theta - (2*mean(cells)-1).^2;

% First run to populate the moving average test
len_test = 15;
ini_offset = 5; % start test at t=5
xmean_last = zeros(len_test, 1);
xstd_last = xmean_last;
corr_last = xmean_last;

% start populating moving average
xmean_last(end) = mean(cells);
xstd_last(end) = std(cells);
corr_last(end) = Theta(1) - (2*xmean(1)-1).^2;

for t=1:len_test
    ind_aux = len_test - t+1; % populate in reverse order (older last)
    [cells, ~] = ...
        update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, 10);
    [~, theta] = moranI(cells, a0*dist);
    Theta(end+1) = theta;
    % update test variables
    xmean_last(ind_aux) = mean(cells);
    xstd_last(ind_aux) = std(cells);
    corr_last(ind_aux) = theta - (2*mean(cells)-1).^2;
    % store variables
    cells_hist{t+1} = cells;
    xmean(t+1) = mean(cells);
    xstd(t+1) = std(cells);
    corr(t+1) = theta - (2*mean(cells)-1).^2;
    % update figure
    update_cell_figure_continuum(hin, pos, a0, cells, cell_type, t);
    %frames(end+1) = getframe(gcf);
end

%%
% Check if errors fall within the thresholds
cont = true;
while cont
    if std(xmean_last) < err_thr_mean && std(xstd_last) < err_thr_std
        cont = false;
    elseif t > tmax
    %if t > tmax
        disp('DNF; t_max reached');
        cont = false;
    else
        t = t+1;
        [cells, ~] = ...
        	update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, 10);
        % rotate the list to replace the last 
        xmean_last = circshift(xmean_last, 1);
        xstd_last = circshift(xstd_last, 1);
        % update test variables
        xmean_last(1) = mean(cells);
        xstd_last(1) = std(cells);
        % store variables
        cells_hist{t+1} = cells;
        xmean(t+1) = mean(cells);
        xstd(t+1) = std(cells);
        [~, Theta(end+1)] = moranI(cells, a0*dist);
        % update figure
        update_cell_figure_continuum(hin, pos, a0, cells, cell_type, t);
        %frames(end+1) = getframe(gcf);
    end
end
%}
%% Plot muX(t)=<Xi>(t) and sigmaX(t)
%{
h=figure(3);
hold on
mksize = 8;
l1a = mean(xmean)+err_thr_mean;
l1b = mean(xmean)-err_thr_mean;
l2a = mean(xstd)+err_thr_std;
l2b = mean(xstd)-err_thr_std;
plot(0:t, xmean, 'bo-', 'MarkerSize', mksize);
plot(0:t, xstd, 'ro-', 'MarkerSize', mksize);
%plot(0:t, corr, 'yo-', 'MarkerSize', mksize);
%plot([0 t], [l1a l1a] , 'b--');
%plot([0 t], [l1b l1b], 'b--');
%plot([0 t], [l2a l2a], 'r--');
%plot([0 t], [l2b l2b], 'r--');
set(gca, 'FontSize', 20);
xlabel('time');
ylabel('$$\langle X_i \rangle$$'); % correct label
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