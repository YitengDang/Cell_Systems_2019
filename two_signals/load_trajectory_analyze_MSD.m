% Load trajectory, analyze its p1-p2 dynamics using the Moment Scaling
% Spectrum
clear all
close all
clc

%% Load trajectory
%folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\Fig_2_Dropbox';
%folder = 'L:\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\patterns\Network 33';
%folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\travelling waves';
%folder = 'H:\My Documents\Multicellular automaton\temp';
%folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\Fig5-sample complex trajectory';
%folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2\selected\filmstrip_selection';
%folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\data\Fig_S6_types_of_TW';
%folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\FigS3_class_III_oscillating';
%folder = 'H:\My Documents\Multicellular automaton\app\data\time_evolution\patterns\Network 12';
%folder = 'H:\My Documents\Multicellular automaton\paper_2\figures\originals\Fig6-extensions\wave_examples';
folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity paper 2\movies\Fig_6_model_extension_patterns';
%fname_str = 'travelling_pulse_horizontal_M_int1_1_-1_0_t_out_243_period_30-network_19';
%fname_str = 'spiral_single_period_172_M_int_1_-1_1_0_network_15';
%fname_str = 'sample_complex_trajectory';
%fname_str = 'two_signal_mult_N225_ini_state_TW_params_5_mcsteps_400_t_out_15_period_15-v1';
fname_str = 'Fig_6_C_TW_formation_hill_10_ini_state_rand';

load(fullfile(folder, fname_str), 'cells_hist', 'positions', 'distances', 'save_consts_struct');
%load(fullfile(folder, fname_str), 'cells_hist', 'positions', 'distances', 'save_consts_struct', 'positions_all');
rcell = save_consts_struct.rcell;
a0 = save_consts_struct.a0;
t_out = numel(cells_hist)-1;
%% Replay trajectory, calculate p and I

p_all = zeros(numel(cells_hist), 2);
I_all = zeros(numel(cells_hist), 2);

h = figure;
disp_mol = 12;
[h_cells, h_borders]  = reset_cell_figure_minimal(h, positions, rcell);  
for i=0:t_out
    cells = cells_hist{i+1};
    %update_cell_figure_continuum(hin, pos, cells, cell_type, i, disp_mol);
    %update_figure_periodic_scatter(plot_handle, cells, time, disp_mol, showI, a0, distances)
    %update_cell_figure_external(h_cells, h_borders, cells, i, disp_mol, positions);    
    
    % save p, I
    p_all(i+1, :) = mean(cells, 1);
    I_all(i+1, 1) = moranI(cells(:,2), a0*distances);
    I_all(i+1, 2) = moranI(cells(:,2), a0*distances);
    
    %pause(0.5);
end

%% Plot p1-p2 dynamics
% plot frames
h=figure;
box on
%clf(h, 'reset');
p00 = scatter( p_all(1, 1), p_all(1, 2), 100, 'gs', 'filled');
hold on
p01 = scatter( p_all(end,1), p_all(end,2), 100, 'rs', 'filled');
set(h, 'Units', 'Inches', 'Position', [1 1 10 8]);
xlim([0 1])
ylim([0 1])
xlabel('p^{(1)}', 'Interpreter', 'tex');
ylabel('p^{(2)}', 'Interpreter', 'tex');
legend([p00 p01], {'Initial state', 'Final state'}, 'AutoUpdate','off');
set(gca, 'FontSize', 24);
for t_final=0:t_out
    p1 = plot( p_all(1:(t_final+1),1), p_all(1:(t_final+1),2), 'b-' );
    p1.Color(4) = 0.6;
    p2 = scatter( p_all(t_final+1,1), p_all(t_final+1,2), 100, 'ko', 'filled');
    uistack(p00,'top'); % move plots to top
    uistack(p01,'top');
    
    title(sprintf('Time = %d', t_final));
    
    pause(0.001);
    delete(p1);
    delete(p2);
end
hold off 

%% Calculate gamma_P

p_pow = 2;
tau_all = 1:20;
xp_tau = zeros(numel(tau_all), 1); % <x^p>(tau)
for i=1:numel(tau_all)
    tau = tau_all(i);
    dx_all = (p_all(1+tau:end, :)-p_all(1:end-tau, :));
    dx_norm = sqrt( dx_all(:,1).^2 + dx_all(:,2).^2 );
    xp_tau(i) = sum( dx_norm.^p_pow  )/length(dx_norm);    
end

h = figure;
plot( tau_all, xp_tau, 'bo-' );

%plot( log(tau_all), log(xp_tau), 'o-' )
%xlabel('log(\tau)');
%ylabel('log(\langle |x|^p \rangle)');
%% Calculate 