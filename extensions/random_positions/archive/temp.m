clear all;
close all
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\trav_wave_stability_general\run2_net_parameters_TW_sim';
fname_str = 'Wave_type_1_network_36_states_F2_M4_B3_E1_Con_K_values_waves_sim';

load(fullfile(folder, fname_str));

folder2 = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\batch_sim_all_topologies_run2';
fname_str2 = 'batch_sim_analyzed_data_batch2';
load(fullfile(folder2, fname_str2), 'M_int_all_reduced');

network = 19;
M_int = M_int_all_reduced{network};

save(fullfile(folder, fname_str), 'a0', 'hill', 'lambda', 'M_int', 'N',...
    'network', 'noise', 'rcell', 'K_wave_sim', 'Con_wave_sim');