clear all
close all
clc
%% Load trajectory
%folder = 'H:\My Documents\Multicellular automaton\data\two_signals\parameter_gradient\vertical_step_function';
%fname_str = 'Parameter_gradient_K_2_1_square_wave_t_out_1000_period_16-v1_trav_wave_vertical';

folder = 'H:\My Documents\Multicellular automaton\data\two_signals\parameter_gradient\horizontal_step_function_Ax_0p5';
fname_str = 'Parameter_gradient_K_2_1_square_wave_t_out_5000_period_24-v1';
load( fullfile(folder, fname_str) );

%%
[period_ub, ~] = periodicity_test_short(cells_hist);
[period, t_onset] = periodicity_test_detailed(cells_hist, t_out,...
    period_ub);
t_wave = t_onset;

%% plot snapshot
t_wave = t_onset + 1;

cells = cells_hist{t_wave+1};
disp_mol = 12;
showI = 0;
a0 = save_consts_struct.a0;
rcell = save_consts_struct.rcell;
h = figure;
plot_handle = reset_cell_figure(h, positions, rcell);
update_figure_periodic_scatter(plot_handle, cells, t_wave, disp_mol, showI, a0, distances);

%%
orientation_all = {};
for t_wave=t_onset:t_onset+period
    try 
        wave_state = -1;
        [orientation, bands_in_wave, number_of_waves, diag_vec, band_vec,...
            wave_state, bended] = determine_wave_properties(cells_hist, t_wave, wave_state);
        [wave_state_translated] = translate_states(wave_state);

        orientation_all{end+1} = orientation;
    catch
        disp('Not a travelling wave');
    end
end
%%
h=figure;
histogram(categorical(orientation_all));