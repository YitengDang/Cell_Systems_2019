% Computes the spectral density of a trajectory
clear all
close all
set(0, 'defaulttextinterpreter', 'latex');

%% Load trajectory
folder = 'K:\bn\hy\Shared\Yiteng\Multicellularity\videos\selected\data';
%folder = 'K:\bn\hy\Shared\Yiteng\Multicellularity\videos\synchronization';

% Input specific filename 
fname_str = 'Chaos_in_growing_domain';
load(fullfile(folder, fname_str), 'cells_hist');

%{
% Get all file names in the directory
listing = dir(folder);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
    end
end
%}
%%
%for file_idx=1:count
    %{
    fname_str = names{file_idx};
    fprintf('%s\n', fname_str);
    load(fullfile(folder, fname_str), 'cells_hist');
    %}

    % Process loaded data
    t_out = numel(cells_hist)-1;

    % Take a slice
    t_sel = 1:t_out;
    cells_hist = cells_hist(t_sel+1);

    % Variables
    A = cell2mat(cells_hist);
    X1_all = A(:, 1:2:end); % X^(1)(t) for all cells, N x t_out+1
    X2_all = A(:, 2:2:end); %

    % normalize data by mean
    X1_all = X1_all - mean(X1_all, 1);
    X2_all = X2_all - mean(X2_all, 1);

    % other parameters
    n = t_out+1;       % number of samples
    f = (0:n-1)/n;     % frequency range of FFT
    f = fftshift(f); f(f>=1/2) = f(f>=1/2)-1; % shifted frequency range
    %% Single cell trajectory
    %{
    cell_idx = 1;
    Xf1 = fft(X1_all(cell_idx,:));
    %Xf2 = fft(X2_all(cell_idx,:));

    n = t_out+1;          % number of samples
    f = (0:n-1)/n;     % frequency range
    power1 = abs(Xf1).^2/n;    % power of the DFT

    figure;
    plot(f, power1)
    xlabel('Frequency')
    ylabel('Power')

    figure;
    hold on
    plot(1./f, power1)
    plot([4 4], [0 5], 'r--');
    xlabel('Period')
    ylabel('Power')
    set(gca, 'XScale', 'log');
    %}
    %% Averaged over all cells
    Xf1 = fftshift(fft(X1_all, n, 2));
    Xf2 = fftshift(fft(X2_all, n, 2));

    power1 = mean(abs(Xf1).^2/n, 1);    % power of the DFT
    power2 = mean(abs(Xf2).^2/n, 1);    % power of the DFT

    h1=figure;
    hold on
    plot(f, power1)
    plot(f, power2)
    legend({'1', '2'});
    xlabel('Frequency')
    ylabel('Power')
    set(gca, 'FontSize', 20);
    %%
    idx = f>0; % choose only frequencies within range 
    
    h2=figure;
    hold on
    p1=plot(1./f(idx), power1(idx));
    p2=plot(1./f(idx), power2(idx));
    ymax = ceil(max(max(power1(idx)), max(power2(idx))));
    plot([2 2], [0 ymax], 'k--');
    plot([4 4], [0 ymax], 'r--');
    legend([p1 p2], {'1', '2'});
    xlabel('Period')
    ylabel('Power')
    set(gca, 'XScale', 'log');
    set(gca, 'FontSize', 20);

    %% Save spectral densities
    qsave = 1;
    save_folder = 'K:\bn\hy\Shared\Yiteng\Multicellularity\videos\selected\spectral densities';
    %fname_str = '';
    t_s = sprintf('_t%d_to_%d', t_sel(1), t_sel(end));
    fname_out = fullfile(save_folder, strcat(fname_str, '_spectral_density')); %, t_s));
    save_figure(h1, 10, 8, fname_out, '.pdf', 0)

    fname_out = fullfile(save_folder, strcat(fname_str, '_spectral_density_vs_period')); %, t_s));
    save_figure(h2, 10, 8, fname_out, '.pdf', qsave)
    
    close all
%end