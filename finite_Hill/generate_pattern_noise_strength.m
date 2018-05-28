%% Measures the deviations obtained from perturbing the system with a given noise strength
clear all
close all
clc
%%
gz = 11;
N = gz^2;
n_smpl = 100; 
n_noise = 100; 
p0 = 0.5;
%iniON = round(p0*N);
noisevals = 10.^[-3:0.2:1];

% (1) uniform initial config
%cells = rand(N, 1);
% (2) binary initial config
%cells = zeros(N, 1);
%cells(randperm(N, iniON)) = 1;
% (3) monte carlo initial config
%cells = init_random_cells_montecarlo(N, k/N);

%% Calculate absolute changes
d1vals = zeros(numel(noisevals), n_smpl, n_noise);
for i0=1:numel(noisevals)
    noise = noisevals(i0);
    for i1=1:n_smpl
        cells_ini = init_random_cells_montecarlo(N, p0);
        for i2=1:n_noise
            cells = generate_pattern_noise(cells_ini, noise);
            d1vals(i0,i1,i2) = sum(abs(cells-cells_ini));
        end
    end
end
%% mean
h=figure();
d1mean = mean(mean(d1vals, 3), 2);
d1std = sqrt(sum(std(d1vals/N, 0, 3).^2, 2))/n_smpl;
plot = errorbar(noisevals, d1mean/N, d1std, 'LineWidth', 2);
set(get(plot,'Parent'), 'XScale', 'log');
xlabel('noise');
ylabel('$$\langle d_1 \rangle$$', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [1 1 10 8]);
%errorbar(noisevals, d1mean, d1std);

qsave = 1;
if qsave
    fname_out = strrep(sprintf('analysis_noise_strength_mean_montecarlo_p0_%.2f', p0),...
        '.', 'p');
    out_file = fullfile(pwd, 'figures', 'sensitivity_initial_cond_continuum',...
        fname_out); %filename
    save_figure_pdf(h, 10, 8, out_file);
    save_figure_eps(h, 10, 8, out_file);
end

%%