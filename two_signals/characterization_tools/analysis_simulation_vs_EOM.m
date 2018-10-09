% analyze results to compare simulation and EOM
clear all
close all
%% Load simulation data
% NB. two simulations should have the same tmax

% Simulation
folder = 'H:\My Documents\Multicellular automaton\temp';
fname_str = 'plane_wave_formation_period_15';
fname = fullfile(folder, fname_str);
load(fname, 'cells_hist', 't_out');
cells_out_sim = cells_hist{end};
cells_idx = cells_out_sim*[1; 2];
n_out_sim = reshape(histcounts(cells_idx, -0.5:3.5), 2, 2);
t_out_sim = t_out;
N = size(cells_out_sim, 1);
n_all_sim = zeros(t_out_sim+1, 1);
for i=1:t_out_sim+1
    
end
%% Load EOM data (multiple)
n_sim = 5;
d_t_out = zeros(n_sim, 1);
d1 = zeros(n_sim, 1);
d2 = zeros(n_sim, 1);

for v=1:n_sim
    folder = 'H:\My Documents\Multicellular automaton\temp';
    fname_str_temp = sprintf('%s_EOM-v%d', fname_str, v);
    fname = fullfile(folder, fname_str_temp);
    load(fname, 'n_all', 't_out');
    n_out_EOM = n_all(:, :, end);
    t_out_EOM = t_out;
    
    % Compare t_out
    d_t_out(v) = t_out_EOM - t_out_sim;

    % Compare n_out
    % (1) Hamming dist
    d1(v) = sum(abs(n_out_EOM(:) - n_out_sim(:)))/4/N;

    % (2) Euclidean dist
    d2(v) = sum((n_out_EOM(:)/N - n_out_sim(:)/N).^2);

end

fprintf('<d(t_out)> = %.2f \n', mean(d_t_out));
fprintf('<d1(initial state, final state)> = %.2f \n', mean(d1));
fprintf('<d2(initial state, final state)> = %.2f \n', mean(d2));

%% Plot histogram of d(t_out)
figure;
histogram(d_t_out);

%% Plot histogram of final state differences

figure;
%hold on
histogram(d1, [0:0.1:1.5]);

figure;
histogram(d2, [0:0.1:1]);
