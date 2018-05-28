function [t_mat] = transition_matrix_func(gridsize, Con, K, a0, Rcell)

% path to search for saved results
data_path = fullfile(pwd, 'data', 'transition_matrix');

% Parameters
%gridsize = 11;
[dist, ~] = init_dist_hex(gridsize, gridsize);
N = gridsize^2;
% Feedback
%Con = 15;
%K = 6;
%a0 = 1.5;
%Rcell = 0.2*a0;
% Spatial order
%I = 0.5;

filename = sprintf('t_mat_grid%d_Con%d_K%d_a0%d_R%d.mat', ...
    gridsize, round(Con), round(K), round(10*a0), ...
    round(100*Rcell));

fname = fullfile(data_path, filename);

if exist(filename, 'file') == 2
    load(fname)
else
    t_mat = zeros(N+1);
    for n = 0:N
        %disp(n)
        [ptsum, ~, ~] = transition_prob(dist, a0, Rcell, K, Con, n);
        t_mat(n+1, :) = ptsum;
    end
    save(fname)
end

end