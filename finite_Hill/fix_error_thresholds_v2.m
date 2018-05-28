function [err_thr_mean, err_thr_std, err_thr_corr] = fix_error_thresholds_v2(dist, a0, Con, K, fN, hill, noise, cm, cs)
N = size(dist, 1);
Rcell = 0.2*a0;

% Find uniform lattice fixed points
fp = zeros(3, 1);
x0 = [0.0 0.5 1]; %estimates based on previous graph
hfunc = @update_function_uniform;
hfunc2 = @(x) hfunc(x, hill, Con, K, fN);
for idx=1:3
    fp(idx) = fzero(hfunc2, x0(idx));
end
%}
%% Fix the error thresholds by running simulation starting in fixed point with noise
t_test = 1100;
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

for t=1:t_test
    [cells_test, ~] = ...
        update_cells_noise_hill(cells_test, dist, Con, K, a0, Rcell, noise, hill, 10);
    [~, Theta] = moranI(cells_test, a0*dist);
    cells_test_hist{t+1} = cells_test;
    xmean(t+1) = mean(cells_test);
    xstd(t+1) = std(cells_test);
    corr(t+1) = Theta - (2*xmean(1)-1).^2;
end

%figure();
%hold on
%plot(0:t_test, xmean);
%plot(0:t_test, xstd);
%plot(0:t_test, corr);

offset = 100;
%cm = 2;
%cs = 2;
err_thr_mean = cm*std(xmean(offset:end));
err_thr_std = cs*std(xstd(offset:end));
err_thr_corr = cs*std(corr(offset:end));
%fprintf('Error threshold mean %d \n ', err_thr_mean);
%fprintf('Error threshold std %d \n', err_thr_std);