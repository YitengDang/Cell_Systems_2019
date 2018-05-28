% This script calculates the map of p_in p_eq using exact simulation.
close all
clear variables

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 5.6;
Rcell = 0.2*a0;
K = 8;
Con = 16;
hill = 2;
%noise = 1;
noiselist = 10.^[-2:0.25:0.5];

% method for sampling initial states; binary=ON/OFF, normal = N(p, sigma)
sampmethod = 'montecarlo'; %'binary' 'normal' 'uniform'

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
dist_vec = dist(1, :);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength

%%
for run=1:numel(noiselist)
    noise = noiselist(run);
    
% Load file
fname = strrep(sprintf('pin_pout_noise_%.4f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_%s', ...
    noise, N, Con, K, a0, hill, sampmethod), '.', 'p');
path = 'H:\My Documents\Multicellular automaton\finite_Hill\figures\pin_pout\data\old_noise';
%[fname, path, ~] = uigetfile(path);
load(fullfile(path, strcat(fname, '.mat')));
close all;

%% redo p=1 simulation
count(N+1, :) = zeros(1,N+1);
t_av(N+1) = 0;
corr_av(N+1,:) = zeros(1,3);

% simulation parameters
n_smpl = 1000; % number of samples for each value
tmax = 1000;
len_test = 10;

% fix error thresholds
[err_thr_mean, err_thr_std, ~] = fix_error_thresholds(dist, a0, Con, K, fN, hill, noise);

k=N; % k = <Xi(t=0)> = p_initial
disp(k)
ntmax_aux = zeros(n_smpl,1);
corr_aux = zeros(n_smpl,1);
k_aux = zeros(n_smpl,1);
t_aux = zeros(n_smpl,1);
%-----------------
% Test n_smpl random initial configurations and update cells
% until they are in final configuration
parfor i = 1:n_smpl
    % initiate cells
    cells = init_random_cells_montecarlo(N, k/N);

    % first populate the test
    xmean_last = zeros(len_test, 1);
    xstd_last = xmean_last;
    for t=1:len_test
        ind_aux = len_test - t+1; % populate in reverse order (older last)
        [cells, ~] = ...
            update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, 10);
        % update test variables
        xmean_last(ind_aux) = mean(cells);
        xstd_last(ind_aux) = std(cells);
    end

    % run until errors below thresholds
    cont = true;
    while cont
        if std(xmean_last) < err_thr_mean && std(xstd_last) < err_thr_std
            cont = false;
        elseif t > tmax
            disp('DNF; t_max reached');
            cont = false;
            ntmax_aux(i) = 1;
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
        end
    end
    %-----------------
    t_aux(i) = t;
    kout = round(sum(cells)); % round to nearest integer value
    k_aux(i) = kout;
    [~, theta] = moranI(cells, a0*dist);
    corr_aux(i) = theta - (2*mean(cells)-1).^2;
end
%--- count ntmax---
ntmax(k+1) = sum(ntmax_aux);
%-----------------
% update count
for i = 0:N
    count(k+1,i+1) = count(k+1,i+1) + sum(sum(round(k_aux) == i));
end
% mean final t
t_av(k+1) = sum(sum(t_aux))/n_smpl;
% mean, min and max final I
corr_av(k+1,1) = sum(sum(corr_aux))/n_smpl;
corr_av(k+1,2) = min(corr_aux);
corr_av(k+1,3) = max(corr_aux);

prob = transpose(count./repmat(sum(count,2),1,N+1));
%% save new file
close all
fname_str = strrep(sprintf('pin_pout_noise_%.4f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_%s', ...
    noise, N, Con, K, a0, hill, sampmethod), '.', 'p');
file_out = fullfile(pwd, 'figures', 'pin_pout', 'data', strcat(fname_str,'.mat'));
save(file_out);

%% plot the map
set(0,'defaulttextinterpreter', 'latex');
h1 = figure(1);
im_fig = imagesc(p,p,prob);
% set title and font
%title(sprintf('$$N = %d, K = %d, S_{ON} = %.1f, a_0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18)
set(gca,'Ydir','normal','FontSize', 24)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', count' > 0);
% set colorbar and labels
c = colorbar;
c.Label.String = 'Probability';
xlabel('$$\langle X \rangle_{in}$$', 'FontSize', 24)
ylabel('$$\langle X \rangle_{out}$$', 'FontSize', 24)

% Organize and save
save_fig = 1;
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout', 'noise', strcat(fname_str, '_map'));
    save_figure_pdf(h1, 10, 8, out_file);
    save_figure_eps(h1, 10, 8, out_file);
end

% Solve for uniform lattice fixed points and include in plot
fp = zeros(3, 1);
x0 = [0.03 0.2 0.65]; %estimates based on previous graph
hfunc = @update_function_uniform;
hfunc2 = @(x) hfunc(x, hill, Con, K, fN);
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end
figure(h1);
hold on
plot([0 1], [fp(1) fp(1)], 'r--'); % stable fixed point
%plot([0 1], [fp(2) fp(2)], 'r--');
plot([0 1], [fp(3) fp(3)], 'r--'); % stable fixed point
plot([fp(2) fp(2)], [0 1], 'g--'); % unstable fixed point

% Organize and save
save_fig = 1;
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout', 'noise', strcat(fname_str, '_map_v2'));
    save_figure_pdf(h1, 10, 8, out_file);
    save_figure_eps(h1, 10, 8, out_file);
end
%%
% Plot the average number of steps it takes to reach equilibrium
%
h2 = figure(2);
plot(p, t_av, 'r-o')
set(gca,'FontSize', 24)
%title(sprintf('$$N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f$$', ...
%    N, K, Con, a0, Rcell),'FontSize', 18)
xlabel('$$\langle X \rangle_{in}$$', 'FontSize', 24)
%ylabel('Average # steps for eq.', 'FontSize', 24)
ylabel('$$\langle t_{eq} \rangle$$');
pmax = p(t_av == max(t_av));
title(sprintf('max $$t_{eq}$$ at $$p_{in}$$ = %.2f', pmax));

save_fig = 1;
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout', 'noise', strcat(fname_str, '_tav'));
    save_figure_pdf(h2, 10, 8, out_file);
    save_figure_eps(h2, 10, 8, out_file);
end
%}
%%
% Plot the average corr in equilibrium
h3 = figure(3);
plot(p, corr_av, 'r-o')
set(gca,'FontSize', 24)
%title(sprintf('$$N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f$$', ...
%s    N, K, Con, a0, Rcell),'FontSize', 18)
xlabel('$$\langle X \rangle_{in}$$', 'FontSize', 24)
ylabel('$$\langle C(X_{out}) \rangle$$', 'FontSize', 24)

save_fig = 1;
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout', 'noise', strcat(fname_str, '_corr'));
    save_figure_pdf(h3, 10, 8, out_file);
    save_figure_eps(h3, 10, 8, out_file);
end

end