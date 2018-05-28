% Monte Carlo algorithm to generate random configuration with fixed N
% and fraction of ON cells (determined by basin of attraction, determined by fixed points)
% Works for Hill coefficient n=2 
close all
clear variables
clc

% Lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 5.4;
Rcell = 0.2*a0;

% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist = round(dist, 5);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% Find fixed points
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
%hfunc = @update_function_uniform;
hfunc2 = @(x) (((Con-1)*x + 1))^hill/(K^hill+ ((Con-1)*x + 1)^hill) - x;
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end

%%
% If p is not 0 or 1
p = 0.5;
mcsteps = 10^5;
ntrials = 10^4;
n = round(p*N);
Non = zeros(ntrials, 2); % 1: initial, 2: final
var_vals_ini = zeros(ntrials, 2); % (_,1): ON cells, (_,2): OFF cells
var_vals_fin = zeros(ntrials, 2);

for j=1:ntrials
    r = randi(N, mcsteps, 1);
    cells = fp(1)*ones(N, 1);
    cells( randperm(N, n) ) = fp(3);
    Non(j, 1) = sum(cells > fp(2));
    var_vals_ini(j, 1) = var( cells(cells>fp(2)) );
    var_vals_ini(j, 2) = var( cells(cells<fp(2)) );
    accepted = 0;
    for i=1:mcsteps
        %delta = (2*rand()-1)/2; % -0.5 <= delta <= 0.5
        alpha = min(fp(3)-fp(2), fp(2)-fp(1))/3; % arbitrary choice, tuned to give satisfactory fraction of accepted statess
        delta = alpha*randn(); % delta ~ N(0, minimum distance to unstable fixed point) 
        cell_new = cells(r(i))+delta; 
        state_old = cells(r(i)) > fp(2); 
        state_new = cell_new > fp(2);
        if all(state_old == state_new) && all(abs(2*cell_new-1) <= 1) % if both cells still in range [0, 1]
            accepted = accepted+1;
            cells(r(i)) = cell_new;
        end
    end
    fprintf('accepted = %d out of %d steps \n', accepted, mcsteps);
    %Non(j, 2) = sum(cells > fp(2));
    var_vals_fin(j, 1) = var( cells(cells>fp(2)) );
    var_vals_fin(j, 2) = var( cells(cells<fp(2)) );
end

dir = fullfile(pwd, 'figures\pin_pout\data_binary\init_random_cells_binary_constraint');
fname = strrep(sprintf('Trialdata_onOff_N%d_K%d_Con%d_a0%.1f_hill%.1f_p%.2f_mcsteps%d_ntrials%d', ...
    N, K, Con, a0, hill, p, mcsteps,  ntrials), '.', 'p');
save( fullfile(dir, fname) );

%%
%figure();
%histogram(cells, 0:0.05:1);
%% Load data
mcvals = 10.^(0:5);
var_ini_mean = zeros( numel(mcvals), 2); %(_,1): ON cells, (_,2): OFF cells
var_fin_mean = zeros( numel(mcvals), 2);
var_ini_err = zeros( numel(mcvals), 2);
var_fin_err = zeros( numel(mcvals), 2);

dir = fullfile(pwd, 'figures\pin_pout\data_binary\init_random_cells_binary_constraint');
for i=1:numel(mcvals)
    % clear("var_vals_ini", "var_vals_fin");
    fname = strrep(sprintf('Trialdata_onOff_N%d_K%d_Con%d_a0%.1f_hill%.1f_p%.2f_mcsteps%d_ntrials%d', ...
        N, K, Con, a0, hill, p, mcvals(i),  ntrials), '.', 'p');
    load(fullfile(dir, fname), 'var_vals_ini', 'var_vals_fin');
    var_ini_mean(i, :) = mean(var_vals_ini);
    var_fin_mean(i, :) = mean(var_vals_fin);
    var_ini_err(i, :) = std(var_vals_ini)/sqrt(ntrials);
    var_fin_err(i, :) = std(var_vals_fin)/sqrt(ntrials);
end

%% Plot figure
h1=figure(1);
errorbar(mcvals, var_fin_mean(:,1), var_fin_err(:,1), 'bo--', 'LineWidth', 2);
set(gca,'XScale','log');
hold on
errorbar(mcvals, var_fin_mean(:,2), var_fin_err(:,2), 'ro--', 'LineWidth', 2 );
%errorbar(mcvals, var_ini_mean(:,1), var_ini_err(:,1), 'bo--' );
%errorbar(mcvals, var_ini_mean(:,2), var_ini_err(:,2), 'ro--' );
%semilogx( mcvals, var_fin_mean, 'bo-' );
legend({'ON','OFF'}, 'Location', 'nw');
xlabel('MC steps');
ylabel('Variance in cell states');
set(gca,'FontSize', 24);
