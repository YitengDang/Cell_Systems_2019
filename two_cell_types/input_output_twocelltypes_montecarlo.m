% Batch simulations of trajectories with different inputs (p1, p2, (theta))
% vs outputs
% Save only the main statistics of the simulations, not detailed
% trajectories (too many data points)
close all
clear variables
warning on
set(0,'defaulttextinterpreter', 'latex');
%--------------------------------------------------------------------------
% Set parameters of the system
gridsize = 15;
N = gridsize^2;
frac_type1 = 0.5;
N1 = round(frac_type1*N);
N2 = N - N1;
a0 = 1.5;

% Parameters
K = [10 10];
Con = [20 20]; %column j = type j 
Coff = [1 1];
%p = [0.5 0.5];
rcell = [0.2 0.2];
Rcell = rcell*a0;
%K_loop = 10; %[2 5 10]; %2:15;
%Con_loop = 10; %5:5:40; %2:1:40;

% Simulation parameters
p1 = 0:0.1:1;
p2 = 0:0.1:1;
%[p1_all, p2_all] = meshgrid(p1, p2);
%theta_in = []
theta_str = 'theta_random';

nruns = 10;
tmax = 1000;
%save_trajectory = 0; %whether to save individual trajectories 

% Communication matrix
Mcomm = [1 1; 1 1]; % type i reacts to type j iff M_ij=1
Mcomm_str = sprintf('Mcomm_%d_%d_%d_%d', Mcomm(1,1), Mcomm(1,2), Mcomm(2,1), Mcomm(2,2));
Mcomm = logical(Mcomm);

% Check input
if ~all(size(Mcomm)==[2 2])
    disp('Wrong communication matrix size!');
    return
end
for idx=1:numel(Mcomm)
    if ~(Mcomm(idx)==0 || Mcomm(idx)==1)
        disp('Wrong communication matrix entries!');
        return
    end
end

%% output
fname_str = strrep(sprintf('IO_MC_N1_%d_N2_%d_a0_%.1f_K1_%d_K2_%d_Con1_%d_Con2_%d_nruns_%d_%s_%s_v3_limit_theta_range',...
    N1, N2, a0, K(1), K(2), Con(1), Con(2), nruns, Mcomm_str, theta_str), '.', 'p');

%savepath = fullfile(pwd, 'data', 'input-output_maps', 'montecarlo');
%savepath_fig = fullfile(pwd, 'figures', 'input-output_maps', 'montecarlo', Mcomm_str);
savepath = 'H:\My Documents\Multicellular automaton\temp';
savepath_fig = savepath;
%--------------------------------------------------------------------------
%% Get interaction matrices
% Distance and position
[dist, pos] = init_dist_hex(gridsize, gridsize);
[f_mat, g_mat, f_std, g_std] = calc_interaction_strengths(N, frac_type1, a0, rcell, dist, Mcomm);
%f_mat = [f11 f12; f21 f22];
%g_mat = [g11 g12; g21 g22];
if sum(sum(f_std > 10^(-2)))
    disp('Warning! Interaction strengths have large standard deviation');
end
%% Run trajectories
%
[p1_out, p2_out, theta11_out, theta12_out, theta22_out, t_final_out]...
    = input_output_twocelltypes_montecarlo_run(p1, p2, N1, N2, f_mat, g_mat, Coff, Con, K,...
    tmax, nruns);
%}
% Load trajectories
%load(fullfile(savepath, fname_str));

%%
% save statistics
qsave = 0;
if qsave
    fname = fullfile(savepath, strcat(fname_str, '.mat'));
    clear dist
    save(fname);
end
%% Plot statistics
% plot p_out
h1 = figure(1);
imagesc(p1, p2, (N1*p1_out+N2*p2_out)'/N);
set(gca, 'YDir', 'normal', 'FontSize', 24);
c = colorbar;
xlabel('$$p_1$$');
ylabel('$$p_2$$');
ylabel(c, '$$\langle p_{out} \rangle$$', 'Interpreter', 'latex');
caxis([0 1]);
set(h1, 'Units', 'Inches', 'Position', [1 1 10 8]);
%% plot <p_out> - p_in
[N1_x, N2_x] = meshgrid(p1*N1, p2*N2);
px = (N1_x+N2_x)/N;

% h101: unnormalized
h101 = figure(101);
imagesc(p1, p2, (N1*p1_out+N2*p2_out)'/N - px);
set(gca, 'YDir', 'normal', 'FontSize', 24);
c = colorbar;
xlabel('$$p_1$$');
ylabel('$$p_2$$');
ylabel(c, '$$\langle p_{out} \rangle - p_{in}$$', 'Interpreter', 'latex');
set(h101, 'Units', 'Inches', 'Position', [1 1 10 8]);

% h102: normalized
h102 = figure(102);
imagesc(p1, p2, (N1*p1_out+N2*p2_out)'/N - px);
set(gca, 'YDir', 'normal', 'FontSize', 24);
c = colorbar;
xlabel('$$p_1$$');
ylabel('$$p_2$$');
ylabel(c, '$$\langle p_{out} \rangle - p_{in}$$', 'Interpreter', 'latex');
caxis([-1 1]);
set(h102, 'Units', 'Inches', 'Position', [1 1 10 8]);
colormap('jet');

%% plot I_out
%{
h2 = figure(2);
imagesc(p1, p2, I_out');
set(gca, 'YDir', 'normal', 'FontSize', 24);
c = colorbar;
xlabel('$$p_1$$');
ylabel('$$p_2$$');
ylabel(c, '$$\langle I_{out} \rangle$$', 'Interpreter', 'latex');
caxis([0 0.5]);
set(h2, 'Units', 'Inches', 'Position', [1 1 10 8]);
%}
%% plot p1_out
h3 = figure(3);
imagesc(p1, p2, p1_out');
set(gca, 'YDir', 'normal', 'FontSize', 24);
c = colorbar;
xlabel('$$p_1$$');
ylabel('$$p_2$$');
ylabel(c, '$$\langle p_{1, out} \rangle$$', 'Interpreter', 'latex');
caxis([0 1]);
set(h3, 'Units', 'Inches', 'Position', [1 1 10 8]);

% h301: unnormalized
h301 = figure(301);
imagesc(p1, p2, p1_out' - N1_x/N1);
set(gca, 'YDir', 'normal', 'FontSize', 24);
c = colorbar;
xlabel('$$p_1$$');
ylabel('$$p_2$$');
ylabel(c, '$$\langle p_{1, out} \rangle - p_{1, in}$$', 'Interpreter', 'latex');
set(h301, 'Units', 'Inches', 'Position', [1 1 10 8]);

%% plot p2_out
h4 = figure(4);
imagesc(p1, p2, p2_out');
set(gca, 'YDir', 'normal', 'FontSize', 24);
c = colorbar;
xlabel('$$p_1$$');
ylabel('$$p_2$$');
ylabel(c, '$$\langle p_{2, out} \rangle$$', 'Interpreter', 'latex');
caxis([0 1]);
set(h4, 'Units', 'Inches', 'Position', [1 1 10 8]);

% h401: unnormalized
h401 = figure(401);
imagesc(p1, p2, p2_out' - N2_x/N2);
set(gca, 'YDir', 'normal', 'FontSize', 24);
c = colorbar;
xlabel('$$p_1$$');
ylabel('$$p_2$$');
ylabel(c, '$$\langle p_{2, out} \rangle - p_{2, in}$$', 'Interpreter', 'latex');
set(h401, 'Units', 'Inches', 'Position', [1 1 10 8]);

%% plot t_final_out
h5 = figure(5);
imagesc(p1, p2, t_final_out');
set(gca, 'YDir', 'normal', 'FontSize', 24);
c = colorbar;
xlabel('$$p_1$$');
ylabel('$$p_2$$');
ylabel(c, '$$\langle t_{final, out} \rangle$$', 'Interpreter', 'latex');
%caxis([0 tmax]);
set(h5, 'Units', 'Inches', 'Position', [1 1 10 8]);

%% plot theta11_out
h6 = figure(6);
imagesc(p1, p2, theta11_out');
set(gca, 'YDir', 'normal', 'FontSize', 24);
c = colorbar;
xlabel('$$p_1$$');
ylabel('$$p_2$$');
ylabel(c, '$$\langle \theta^{11}_{out} \rangle$$', 'Interpreter', 'latex');
%caxis([0 tmax]);
set(h6, 'Units', 'Inches', 'Position', [1 1 10 8]);

%% Save figures
qsave = 1;
if qsave
    fname_out1 = fullfile(savepath_fig, strcat(fname_str, '_p_out'));
    %fname_out2 = fullfile(savepath_fig, strcat(fname_str, '_I_out'));
    fname_out3 = fullfile(savepath_fig, strcat(fname_str, '_p1_out'));
    fname_out4 = fullfile(savepath_fig, strcat(fname_str, '_p2_out'));
    fname_out5 = fullfile(savepath_fig, strcat(fname_str, sprintf('_tf_out_tmax_%d', tmax)));

    fname_out101 = fullfile(savepath_fig, strcat(fname_str, '_p_out_diff'));
    fname_out102 = fullfile(savepath_fig, strcat(fname_str, '_p_out_diff_norm'));
    fname_out301 = fullfile(savepath_fig, strcat(fname_str, '_p1_out_diff'));
    fname_out401 = fullfile(savepath_fig, strcat(fname_str, '_p2_out_diff'));

    save_figure_pdf(h1, 10, 8, fname_out1);
    %save_figure_pdf(h2, 10, 8, fname_out2);
    save_figure_pdf(h3, 10, 8, fname_out3);
    save_figure_pdf(h4, 10, 8, fname_out4);
    save_figure_pdf(h5, 10, 8, fname_out5);

    save_figure_pdf(h101, 10, 8, fname_out101);
    save_figure_pdf(h102, 10, 8, fname_out102);
    save_figure_pdf(h301, 10, 8, fname_out301);
    save_figure_pdf(h401, 10, 8, fname_out401);
end