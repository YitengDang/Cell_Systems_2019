a0 = 1.5;
fname_str = 'two_signal_mult_N225_p1_0p50_p2_0p50_K12_22_t_out_66_period_15-v1';
p=[0.5 0.5];
parent_folder = 'L:\BN\HY\Shared\Yiteng\two_signals\travelling_wave_analysis\vs_p0';
subfolder = strrep(sprintf('ini_p1_%.2f_p2_%.2f', p(1), p(2)), '.', 'p');
folder = fullfile(parent_folder, subfolder);
load(fullfile(folder, fname_str));
%%
new_res = travelling_wave_test(cells_hist, a0, period, t_out);
old_res = travelling_wave_test_old(cells_hist, a0, period, t_out);
disp(new_res)
disp(old_res)
%%    
s = size(cells_hist{1}, 2);
gz = 15;
N = gz^2;
[dist, pos] = init_dist_hex(gz, gz);

p_last = zeros(period+1, 2);
I_last = zeros(period+1, 2);

for i1=1:period+1
    i2 = i1 + t_out-period;
    cells = cells_hist{i2};
    p_last(i1, :) = sum(cells, 1)/N;
    for j1=1:s
        [I_last(i1, j1), ~] = moranI(cells(:, j1), a0*dist);
    end
end

I_last = round(I_last, 5); % deal with rounding mistakes
%%
figure;
plot(p_last)
