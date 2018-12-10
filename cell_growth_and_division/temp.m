clear all
folder = 'H:\My Documents\Multicellular automaton\temp\sustained_inhomogeneity';
for ii=1:38
    fname_str = strcat('two_signals_rcell_sigma_0p1_K_growth_1p5_sigma_D_0p100_t_out_1000-v', num2str(ii), '.mat');
    load( fullfile(folder, fname_str) );
    save_consts_struct.sigma_D = 0.1;
    
    save( fullfile(folder, fname_str), 'save_consts_struct', 'cells_hist', 't_out',...
        'changed', 'positions', 'distances', 'positions_all', 'rcell_hist');
end

