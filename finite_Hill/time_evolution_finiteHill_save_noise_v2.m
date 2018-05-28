% Finite Hill coefficient
% v2: new equilibrium condition (based on threshold on fluctuations in mean
% and std)
close all
clear variables
warning off
%%
% Lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 6;
Rcell = 0.2*a0;

% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient
noise = 0.1;

% Initial configuration
initialID = 'uniform'; %'binary' 'uniform' 'allON' 'singleOFF' 'fixedp'
% If fixedp, assume initial states ~ Normal(p,sigma)
p = 0.2;
sigma = 0.1; 
iniON = round(p*N);

% simulation parameters
nsims = 10;
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
%% Find uniform lattice fixed points
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
hfunc = @update_function_uniform;
hfunc2 = @(x) hfunc(x, hill, Con, K, fN);
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end
%}
%% Error threshold method (appendix)
% Fix the error thresholds by running simulation starting in fixed point with noise
%{
t_test = 100;
cells_test_hist = cell(t_test+1, 1);
xmean = zeros(t_test+1, 1);
xstd = zeros(t_test+1, 1);
corr = xmean;
cells_test = fp(end)*ones(N, 1); % start at fixed point
cells_test_hist{1} = cells_test;
xmean(1) = mean(cells_test);
xstd(1) = std(cells_test);
[~, Theta] = moranI(cells_test, a0*dist);
corr(1) = Theta - (2*xmean(1)-1).^2;

for t=0:t_test
    [cells_test, ~] = ...
        update_cells_noise_hill(cells_test, dist, Con, K, a0, Rcell, noise, hill, 10);
    [~, Theta] = moranI(cells_test, a0*dist);
    cells_test_hist{t+1} = cells_test;
    xmean(t+1) = mean(cells_test);
    xstd(t+1) = std(cells_test);
    corr(t+1) = Theta - (2*xmean(1)-1).^2;
end

figure();
hold on
%plot(0:t_test, xmean);
%plot(0:t_test, xstd);
plot(0:t_test, corr);
offset = 10;
cm = 1;
cs = 1;
cc = 1;
err_thr_mean = cm*std(xmean(offset:end));
err_thr_std = cs*std(xstd(offset:end));
err_thr_corr = cc*std(corr(offset:end));
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

%% Loop over simulations

for sim=1:nsims

% initialize fig & vars
hin = figure(1);
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

% update figure
update_cell_figure_continuum(hin, pos, a0, cells, cell_type, t);
    
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
end
% Check if errors fall within the thresholds
cont = true;
while cont
    if std(xmean_last) < err_thr_mean && std(xstd_last) < err_thr_std %&& std(corr_last) < err_thr_corr
        cont = false;
    elseif t > tmax
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
        [~, theta] = moranI(cells, a0*dist);
        Theta(end+1) = theta;
        xmean(t+1) = mean(cells);
        xstd(t+1) = std(cells);
        corr(t+1) = theta - (2*mean(cells)-1).^2;
        % update figure
        update_cell_figure_continuum(hin, pos, a0, cells, cell_type, t);
        %frames(end+1) = getframe(gcf);
    end
end
%}

% Save trajectory
fname_str = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_noise_%.2f_cm%.1f_cs%.1f',...
    N, a0, K, Con, hill, noise, cm, cs), '.', 'p');
i = 1;
fname = fullfile('H:\My Documents\Multicellular automaton\temp',...
    strcat(fname_str,initialID,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile('H:\My Documents\Multicellular automaton\temp',...
        strcat(fname_str,initialID,'-v',int2str(i),'.mat'));
end
save(fname);

end

%% Plot muX(t)=<Xi>(t) and sigmaX(t)
%{
h=figure(2);
hold on
mksize = 8;
l1a = mean(xmean)+err_thr_mean;
l1b = mean(xmean)-err_thr_mean;
l2a = mean(xstd)+err_thr_std;
l2b = mean(xstd)-err_thr_std;
%plot(0:t, xmean, 'bo-', 'MarkerSize', mksize);
%plot(0:t, xstd, 'ro-', 'MarkerSize', mksize);
plot(0:t, corr, 'go-', 'MarkerSize', mksize);
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
test_corr_all = zeros(t-len_test, 1);

for time=len_test+1:t+1
    idx = time-len_test;
    test_mean_all(idx) = std(xmean(time-len_test+1:time));
    test_std_all(idx) = std(xstd(time-len_test+1:time));   
    test_corr_all(idx) = std(corr(time-len_test+1:time));   
end

figure();
hold on
%plot(len_test+1:t+1, test_mean_all, 'b');
%plot(len_test+1:t+1, test_std_all, 'r');
plot(len_test+1:t+1, test_corr_all, 'go-');
%plot([len_test+1 t+1], [err_thr_mean err_thr_mean], 'b--');
%plot([len_test+1 t+1], [err_thr_std err_thr_std], 'r--');
plot([len_test+1 t+1], [err_thr_corr err_thr_corr], 'g--');
legend({'last 10 xmean', 'last 10 xstd', 'xmean threshold', 'xstd threshold'});
xlabel('time');
ylabel('value');
%}