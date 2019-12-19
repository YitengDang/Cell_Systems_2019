% Construct and simulate a cyclic cellular automaton
clear all
close all
clc

% Parameters
l = 300; % size of system, L = l x l
kappa = 4; % number of colors
r = 1; % range of cells
% T = 1; % threshold
%% Single simulation
%{
% Data to store
cells_hist = {};

% Initiate lattice
cells = randi(kappa, l, l)-1;
cells_hist{1} = cells;

% Plot lattice
%{
h=figure;
set(h, 'Position', [100 100 800 700]);
p=imagesc(cells);
pbaspect([1 1 1])
colormap(parula(kappa))
c=colorbar;
title(sprintf('t=%d', 0));
set(c, 'YTick', []);
pause(0.01);
%}
% Run simulation
t_max = 100;
for t=1:t_max
    disp(t);
    cells = update_cells(cells, kappa, r, T);
    cells_hist{end+1} = cells;
    %{
    p.CData = cells;
    title(sprintf('t=%d', t));
    pause(0.01);
    %}
end
%% Replay trajectory
% Plot lattice
h=figure;
set(h, 'Position', [100 100 800 700]);
p=imagesc(cells_hist{1});
pbaspect([1 1 1])
colormap(parula(kappa))
c=colorbar;
title(sprintf('t=%d', 0));
set(c, 'YTick', []);
pause(0.1);

for t=1:t_max
    p.CData = cells_hist{t+1};
    title(sprintf('t=%d', t));
    pause(0.1);
end
%}
%% Batch simulations
sim_count = 500;
T_all = 1:4;
for T_idx = 1:numel(T_all)
    T = T_all(T_idx);
    
    save_folder = 'N:\tnw\BN\HY\Shared\Yiteng\cyclic_ca_2D\batch_2019_05_22';

    fpattern = sprintf('Cyclic_CA_l%d_kappa%d_r%d_T%d_sim_%s.mat', l, kappa, r, T, '(\d+)');
    files = dir(save_folder);
    count = 0;
    for i=3:numel(files)
        filename = files(i).name;
        out = regexp(filename, fpattern, 'match');
        if ~isempty(out)
            count = count+1;
        end
    end
    fprintf('Simulations to do: %d \n', sim_count-count)

    for sim=count+1:sim_count
        disp(sim);
        % Data to store
        cells_hist = {};

        % Initiate lattice
        cells = randi(kappa, l, l)-1;
        cells_hist{1} = cells;

        % Run simulation
        t_max = 100;
        for t=1:t_max
            cells = update_cells(cells, kappa, r, T);
            cells_hist{end+1} = cells;
        end

        % Save simulation
        fname_str = sprintf('Cyclic_CA_l%d_kappa%d_r%d_T%d_sim_%d',...
            l, kappa, r, T, sim);
        fname = fullfile(save_folder, fname_str);
        save(fname, 'cells_hist', 'kappa', 'r', 'T')
    end

end