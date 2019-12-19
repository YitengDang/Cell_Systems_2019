% Compares trajectories of the finite and infinite Hill coefficient models
% with the same initial configuration
% New version: revised metrics
clear all
close all
set(0, 'defaulttextinterpreter', 'tex');
%warning off
%% Parameters
% lattice parameters
gridsize = 11;
N = gridsize^2;
% circuit parameters
Con = 16;
K = 8;
hill = 2;
% initial conditions
initialID = 'binaryrand'; %'uniform'; %'binaryrand'; %'binary'
prec = 8;

% Loop over a0
a0list_run = [5 5.4 5.6 6]; % 5:0.05:5.7];
a0_all = a0list_run; %a0list_run; 

nruns = 100;
tmax = 10000;

%load_folder = fullfile(pwd, 'data', 'dynamics_transition 2017-10-25'); 
parent_folder = 'N:\tnw\BN\HY\Shared\Yiteng\one_signal_finite_Hill\dynamics';
%subfolder = '2017-10-16_vs_a0_4p00_to_6p00_K8_Con16_hill2p00';
%subfolder = '2017-10-25_vs_a0_5p00_to_6p00_K8_Con16_hill2p00_binaryrand';
subfolder = '2018-02-19_a0_5p00_to_6p00_K8_Con16_hill2p00_binaryrand';
load_folder = fullfile(parent_folder, subfolder); 

save_fig_folder = fullfile('H:\My Documents\Thesis\Single gene finite Hill\Fig5');

%% Variables to store
% ---------- static (end) data --------------------
% overall spread
sigma_cells_all = zeros(numel(a0_all), nruns);

% uniform/non-uniform
unif_count_all = zeros(numel(a0_all), nruns);
unif_count_all_2 = unif_count_all;
unif_count_all_3 = unif_count_all;

% d variables
d_all = zeros(numel(a0_all), nruns);
%dmean = zeros(numel(a0_all), 1);
%dstd = dmean;
%unif_frac = zeros(numel(a0_all), 1);

% I2 variables
I2_final_all = zeros(numel(a0_all), nruns);
%{
I2mean = zeros(numel(a0_all), 1);
I2std = I2mean;
I2mean_ord = I2mean;
I2std_ord = I2mean;
I2ordcount = I2mean;
%}

% fractions of ON/OFF cells
p_all = zeros(numel(a0_all), nruns);

% spread within ON/OFF cells 
sigma_cells_ON_OFF_all = zeros(numel(a0_all), nruns, 2); % 1: ON cells, 2: OFF cells
%{
sigma_mean = zeros(numel(a0_all), 2); % column 1: ON cells, column 2: OFF cells
sigma_std = sigma_mean;
p_mean = zeros(numel(a0_all), 1);
p_std = p_mean;
%}

t_out_all = zeros(numel(a0_all), nruns);

% --------- dynamics -----------------
p_dynamics_all = zeros(numel(a0_all), nruns, tmax+1); % fraction ON cells
I2_dynamics_all = zeros(numel(a0_all), nruns, tmax+1); % spatial index
sigma_dynamics_all = zeros(numel(a0_all), nruns, tmax+1, 3); %variability 1: overall std, 2: std OFF cells, 3: std ON cells
d_dynamics_all = zeros(numel(a0_all), nruns, tmax+1); % Hamming distance between initial and current state (rounded)

% histograms
%{
d_edges = 0:5:(floor(N/5)+1)*5;
d_hist = zeros(numel(a0_all), numel(d_edges)-1);

I2_edges = -0.08:0.02:0.18;
I2_hist = zeros(numel(a0_all), numel(I2_edges)-1);

% histogram
sigma_edges = 0:0.02:0.2;
sigma_hist_ON = zeros(numel(a0_all), numel(sigma_edges)-1);
sigma_hist_OFF = sigma_hist_ON;
%}

% Loop over a0
for idx = 1:numel(a0_all)
    % setup
    a0 = a0_all(idx);
    Rcell = 0.2*a0;
    fprintf('a0 = %.2f \n', a0);
    
    % ---------------- Calc single cell / uniform lattice fixed points -----
    % get fN
    [dist, pos] = init_dist_hex(gridsize, gridsize);
    dist = round(dist, 5);
    dist_vec = a0*dist(1,:);
    r = dist_vec(dist_vec>0); % exclude self influence
    fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
    
    % Calc fixed points
    fp = zeros(3, 1);
    x0 = [0.03 0.2 0.65]; %estimates based on previous graph
    %hfunc = @update_function_uniform;
    hfunc2 = @(x) (((Con-1)*x + 1)*(1+fN))^hill/(K^hill+ (((Con-1)*x + 1)*(1+fN))^hill ) - x;
    for i=1:3
        fp(i) = fzero(hfunc2, x0(i));
    end
    
    %% Load data, analyze and plot
    % Path to search for the saved data. It searchs by the name, defined by the
    % parameters chosen
    straux = '(\d+)';
    straux2 = '(\w+)';
    
    if ismember(a0, a0list_run)
        % binaryrand pattern (old)
        %fpattern = strrep(sprintf('N%d_a0_%.2f_Con%.2f_K%.2f_hill%.2f_%s-v%s', ...
        %    N, a0, Con, K, hill, initialID, straux), '.', 'p');
        % uniform pattern / binaryrand new 
        fpattern = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_prec_%d_%s-v%s', ...
            N, a0, K, Con, hill, straux,straux2, prec, initialID, straux), '.', 'p');
    else
        fpattern = strrep(sprintf('N%d_a0_%.3f_Con%.2f_K%.2f_hill%.2f_%s-v%s', ...
            N, a0, Con, K, hill, initialID, straux), '.', 'p');
    end

    % Get all file names in the directory
    listing = dir(load_folder);
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
    
    % Load data
    % Calculate and store Hamming distances between states
    count = 0;
    for i = 1:numel(names)
        [tokens, ~] = regexp(names{i}, fpattern, 'tokens','match');
        if numel(tokens) > 0 && count<nruns
            disp(names{i}) % displays the file name
            count = count+1;
            
            % load the data
            % -- binaryrand (only initial and final cells) --
            %{
            load(fullfile(load_folder,strcat(names{i},'.mat')), 'cells_ini',...
                'cells_final');      
            %}
            % -- uniform (all cells) --
            %
            load(fullfile(load_folder,strcat(names{i},'.mat')), 'cells_hist');
            cells_final = cells_hist{end};
            cells_ini = cells_hist{1};
            t_out_all(idx, count) = numel(cells_hist)-1;
            %}
            
            
            % --------- Dynamic data (all) --------------------------------
            for t=1:numel(cells_hist)
                cells = cells_hist{t};
                
                % Hamming d
                d_dynamics_all(idx, count, t) = sum(abs(round(cells_ini) - round(cells)));
                
                % Spatial index I2
                [~, theta] = moranI(cells, a0*dist);
                I2_dynamics_all(idx, count, t) = theta - (2*sum(cells)/N-1)^2;
                
                % Spread in values (overall)
                sigma_dynamics_all(idx, count, t, 1) = std(cells);
                
                % Binary classification of cells
                cells_ON = (cells > fp(2));
                
                % Store final p and std of ON/OFF cells
                this_p = mean(cells_ON);
                p_dynamics_all(idx, count, t) = this_p;
                if this_p > 0 % only register if there are ON cells
                    sigma_dynamics_all(idx, count, t, 2) = std(cells(cells_ON));
                end
                if this_p < 1
                    sigma_dynamics_all(idx, count, t, 3) = std(cells(~cells_ON));
                end
            end
            
            % fill the remaining values with NaN
            if numel(cells_hist)<tmax+1
                d_dynamics_all(idx, count, numel(cells_hist)+1:tmax+1) = NaN;
                I2_dynamics_all(idx, count, numel(cells_hist)+1:tmax+1) = NaN;
                sigma_dynamics_all(idx, count, numel(cells_hist)+1:tmax+1, :) = NaN;
                p_dynamics_all(idx, count, numel(cells_hist)+1:tmax+1) = NaN;
            end
            
            
            % -------- Static data (start and end) ------------------------
            % d
            %d = sum(abs(cells_ini - round(cells_final)));
            d = sum(abs(round(cells_ini) - round(cells_final)));
            d_all(idx, count) = d;

            % I2
            [~, theta] = moranI(cells_final, a0*dist);
            p_temp = sum(cells_final)/N;
            %I2fin(end+1) = theta - (2*p_temp-1)^2;
            I2_final_all(idx, count) = theta - (2*p_temp-1)^2;
            
            % overall spread
            sigma_cells_all(idx, count) = std(cells_final);
           
            % count nonuniform
            threshold = 0.01;
            unif_count_all(idx, count) = std(cells_final) < threshold ;
            
            % Alternative method
            unif_count_all_2(idx, count) = sum(round(cells_final))==0 || sum(round(cells_final))==N;
            
            % Altenrative 3
            threshold_3 = 0.01;
            unif_count_all_3(idx, count) = (max(cells_final) - min(cells_final)) < threshold;
            
            % get ON and OFF cells
            cells_ON = (cells_final > fp(2));
            this_p = mean(cells_ON);
            p_all(idx, count) = this_p;
            if this_p > 0 % only register if there are ON cells
                sigma_cells_ON_OFF_all(idx, count, 1) = std(cells_final(cells_ON));
            end
            if this_p < 1
                sigma_cells_ON_OFF_all(idx, count, 2) = std(cells_final(~cells_ON));
            end
            %}
        end
    end
    
    % 
    
    
    % fraction uniform
    %unif_frac(idx) = unif_count/count;
    
    % store data
    %{
    % d
    d_hist(idx,:) = histcounts(d_all, d_edges);
    dmean(idx) = mean(d_all);
    dstd(idx) = std(d_all);

    % I2
    I2_hist(idx, :) = histcounts(I2fin, I2_edges); % Heat map of I2 values
    I2mean(idx) = mean(I2fin); 
    I2std(idx) = std(I2fin);
    %}
    
    % I2: Ordered states: define as having sigma < eps
    %{
    eps = 10^(-9);
    I2mean_ord(idx) = mean(I2fin(abs(I2fin)>eps));
    I2std_ord(idx) = std(I2fin(abs(I2fin)>eps));
    I2ordcount(idx) = sum(abs(I2fin)>eps);
    %}
    
    % spread
    %{
    %[sigma_hist_ON(idx, :), sigma_edges] = histcounts(sigma_ON_all);
    sigma_mean(idx, 1) = mean(sigma_ON_all);
    sigma_mean(idx, 2) = mean(sigma_OFF_all);
    sigma_std(idx, 1) = std(sigma_ON_all)/sqrt(numel(sigma_ON_all));
    sigma_std(idx, 2) = std(sigma_OFF_all)/sqrt(numel(sigma_OFF_all));
    p_mean(idx) = mean(p);
    p_std(idx) = std(p)/sqrt(numel(p));
    %}
end

%% Save analyzed data
%{
fname_str = sprintf('analyzed_data_sigma_cells_all_unif_%s_N%d_a0_%.2f_K%d_Con%d_hill%.2f_v2',...
    N, a0, K, Con, hill);
% Save static data
save( fullfile(parent_folder, fname_str), 'a0_all', 'sigma_cells_all',...
    'unif_count_all', 'unif_count_all_2', 'unif_count_all_3', 'd_all', ...
    'I2_final_all', 'p_all', 'sigma_cells_ON_OFF_all');
%}
% Save dynamic data
fname_str = strrep(sprintf('analyzed_data_%s_all_data_a0_%.1f_to_%.1f_4_values_dynamics',...
    subfolder, a0_all(1), a0_all(end)), '.', 'p');
save( fullfile(parent_folder, fname_str), 'a0_all', 't_out_all',...
    'p_dynamics_all', 'I2_dynamics_all', 'sigma_dynamics_all', 'd_dynamics_all');
%}
%% Load saved data
fname_str = strrep(sprintf('analyzed_data_%s_all_data_a0_%.1f_to_%.1f_4_values_dynamics',...
    subfolder, a0_all(1), a0_all(end)), '.', 'p');
load( fullfile(parent_folder, fname_str), 'a0_all', 't_out_all',...
    'p_dynamics_all', 'I2_dynamics_all', 'sigma_dynamics_all', 'd_dynamics_all');

%% Plot spread (lattice uniformity) vs a0
% (1) scatter plot
x_data = repmat(a0_all', 1, nruns);
y_data = sigma_cells_all;
h = figure;
scatter(x_data(:), y_data(:) );
%}
%% (2) density plot
%idx_a0_sel = 2:numel(a0_all)-1;
idx_a0_sel = 1:numel(a0_all);
c_data = sigma_cells_all(idx_a0_sel, :);
bin_edges = 0:0.025:0.30;
bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2;
density_all = zeros(size(c_data, 1), numel(bin_edges)-1 );
for i=1:size(c_data, 1)
    density_all(i, :) = histcounts( c_data(i, :), bin_edges )/nruns;
end

h=figure;
imagesc(a0_all(idx_a0_sel), bin_centers, density_all')
set(gca, 'YDir', 'normal')
xlabel('a_0');
ylabel('Standard deviation \sigma[X_k]'); % \sigma(X_k(t_{eq}))');
c = colorbar;
%ylabel(c, 'P(\sigma[X_k])');
ylabel(c, 'Probability');
set(gca, 'FontSize', 24);

qsave = 0;
fname_str = strrep(sprintf('Sigma_cells_density_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s_%druns',...
    a0_all(1), a0_all(end), N, K, Con, hill, initialID, nruns), '.', 'p');
fname = fullfile(save_fig_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot fraction of NON-UNIFORM lattices
unif_frac = sum(unif_count_all, 2)/nruns;
unif_frac_2 = sum(unif_count_all_2, 2)/nruns;
unif_frac_3 = sum(unif_count_all_3, 2)/nruns;

%idx_a0_sel = 2:numel(a0_all);
idx_a0_sel = 1:numel(a0_all);

h=figure;
hold on
plot(a0_all(idx_a0_sel), 1-unif_frac(idx_a0_sel), 'o-', 'LineWidth', 2);
%plot(a0_all(idx_a0_sel), 1-unif_frac_2(idx_a0_sel), 'o--', 'LineWidth', 2);
%plot(a0_all(idx_a0_sel), 1-unif_frac_3(idx_a0_sel), 'o--', 'LineWidth', 2);
xlabel('a_0');
ylabel('Fraction non-uniform lattices');
set(gca,'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
%legend({'\sigma[X_k]<0.01', 'Round(X_k)'});
%legend({'\sigma[X_k]<0.01', 'max(X_k)-min(X_k) < 0.01'}, 'Location', 'se');
%xlim([a0_all(1)-0.005 a0_all(end)+0.005]);
%set(gca, 'XTick', [4.5 5.0 5.5]);
%ylim([0 70]);
ylim([0 1]);

qsave = 0;
fname_str = strrep(sprintf('Frac_non_unif_lattices_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s_%druns',...
    a0_all(1), a0_all(end), N, K, Con, hill, initialID, nruns), '.', 'p');
fname = fullfile(save_fig_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot mean Hamming distance against a0
%
dmean = mean(d_all, 2);
dstd = std(d_all, 1, 2); %/sqrt(nruns); 

h3=figure;
errorbar(a0_all, dmean, dstd, 'LineWidth', 2);
%plot(a0_all, dmean, 'o-', 'LineWidth', 2);
xlabel('a_0');
ylabel('\langle d_H(X_{ini}, X_{eq}) \rangle');
set(gca,'FontSize', 32);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
%xlim([a0_all(1)-0.05 a0_all(end)+0.05]);
ylim([0 80]);

qsave = 0;
fname_str = strrep(sprintf('Hamming_dist_avg_a0_%.1fto%.1f_N%d_K%d_Con%d_hill_%.1f_%s_%druns_v2_errorbar_std',...
    a0_all(1), a0_all(end), N, K, Con, hill, initialID, nruns), '.', 'p');
fname = fullfile(save_fig_folder, fname_str);
save_figure(h3, 10, 8, fname, '.pdf', qsave);

%% Plot heat map of d_H
%{
h5 = figure;
bincenters = (d_edges(1:end-1)+d_edges(2:end))/2;
nsim = numel(d_all); % number of simulations
imagesc(a0list, bincenters, d_hist'/nsim);
c = colorbar;
xlabel('a_0');
ylabel('d(X_{in}, X_{out})');
ylabel(c, 'Probability');
set(gca,'YDir', 'Normal', 'FontSize', 24);
xticks(5.0:0.1:5.9);
caxis([0 1]);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('d_heat_map_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s_500runs',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_Hamming_dist', fname_str);
    save_figure_pdf(h5, 12, 8, fname);
    save_figure_eps(h5, 12, 8, fname);
end
%}

%% Plot mean I2 against a0
I2mean = mean(I2_final_all, 2);
I2std = std(I2_final_all, 1, 2); %/sqrt(nruns);

h6=figure;
errorbar(a0_all, I2mean, I2std, 'LineWidth', 2);
%plot(a0_all, I2mean, 'o-', 'LineWidth', 2);
xlabel('a_0');
ylabel('\langle I_2(t_{eq}) \rangle');
set(gca,'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
%xlim([a0list(1)-0.05 a0list(end)+0.05]);
xlim([a0_all(1)-0.005 a0_all(end)+0.005]);
ylim([-0.02 0.15]);

qsave = 0;
fname_str = strrep(sprintf('I2_final_avg_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s_%druns_errorbar_std',...
    a0_all(1), a0_all(end), N, K, Con, hill, initialID, nruns), '.', 'p');
fname = fullfile(save_fig_folder, fname_str);
save_figure(h6, 10, 8, fname, '.pdf', qsave);

%% Number of ordered final states
h7 = figure(7);
plot(a0list, I2ordcount, '-o');
xlabel('a0');
ylabel('Number of ordered final states');


%% Heat map
%{
h5 = figure(8);
bincenters = (I2_edges(1:end-1)+I2_edges(2:end))/2;
nsim = numel(I2fin); % number of simulations
imagesc(a0_all, bincenters, I2_hist'/nsim);
c = colorbar;
xlabel('$$a_0$$');
ylabel('$$I_2$$');
ylabel(c, 'Probability');
set(gca,'YDir', 'Normal', 'FontSize', 24);
caxis([0 1]);
xticks(5.0:0.1:5.9);

qsave = 0;
if qsave
    fname_str = strrep(sprintf('I2_final_heat_map_a0_%.2fto%.2f_N%d_K%d_Con%d_hill_%.1f_%s_500runs',...
        a0list(1), a0list(end), N, K, Con, hill, initialID), '.', 'p');
    fname = fullfile(pwd, 'figures', 'autonomy_transition_I2_final', fname_str);
    save_figure_pdf(h5, 12, 8, fname);
    save_figure_eps(h5, 12, 8, fname);
end
%}
%% Plot fraction of ON cells
%
p_mean = mean(p_all, 2);
p_std = std(p_all, 1, 2);

h9 = figure;
hold on
errorbar(a0_all, p_mean, p_std, 'LineWidth', 2)
%plot(a0_all, p_mean, 'o-', 'LineWidth', 2)
plot([a0_all(1) a0_all(end)], [0.5 0.5], 'r--');
xlabel('a_0');
ylabel('Mean fraction ON cells');
set(gca,'FontSize', 32);
set(gcf, 'Units', 'Inches', 'Position', [0.5 0.5 10 8]);
%xlim([a0_all(1)-0.005 a0_all(end)+0.005]);
ylim([0 1.05]); 

qsave = 0;
fname_str = strrep(sprintf('Frac_ON_cells_avg_a0_%.1fto%.1f_N%d_K%d_Con%d_hill_%.1f_%s_%druns_errorbar_std_detailed',...
    a0_all(1), a0_all(end), N, K, Con, hill, initialID, nruns), '.', 'p');
fname = fullfile(save_fig_folder, fname_str);
save_figure(h9, 10, 8, fname, '.pdf', qsave);

%% --------------------Plot dynamics---------------------------------------
for a0_idx=4 %1:4
    
% Plot p(t) 
%a0_idx = 1;
this_a0 = a0_all(a0_idx);
this_tmax = ceil(max(t_out_all(a0_idx, :))/100)*100;

p_data = squeeze(p_dynamics_all(a0_idx, :, 1:this_tmax));
h=figure;
plot(0:this_tmax-1, p_data', 'LineWidth', 1)
ylim([0 1]);
xlabel('t');
ylabel('Fraction ON cells');
set(gca, 'FontSize', 32);

qsave = 0;
fname_str = strrep(sprintf('Frac_ON_cells_vs_t_a0_%.1f_N%d_K%d_Con%d_hill_%.1f_%s_%druns',...
    this_a0, N, K, Con, hill, initialID, nruns), '.', 'p');
fname = fullfile(save_fig_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot I2(t) 
%a0_idx = 1;
this_a0 = a0_all(a0_idx);
this_tmax = ceil(max(t_out_all(a0_idx, :))/100)*100;

I2_data = squeeze(I2_dynamics_all(a0_idx, :, 1:this_tmax));
h=figure;
plot(0:this_tmax-1, I2_data', 'LineWidth', 1)
ylim([-0.2 0.2]);
xlabel('t');
ylabel('I_2');
set(gca, 'FontSize', 32);

qsave = 0;
fname_str = strrep(sprintf('I2_vs_t_a0_%.1f_N%d_K%d_Con%d_hill_%.1f_%s_%druns',...
    this_a0, N, K, Con, hill, initialID, nruns), '.', 'p');
fname = fullfile(save_fig_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot Hamming d
%a0_idx = 1;
this_a0 = a0_all(a0_idx);
this_tmax = ceil(max(t_out_all(a0_idx, :))/100)*100;

d_data = squeeze(d_dynamics_all(a0_idx, :, 1:this_tmax));
h=figure;
plot(0:this_tmax-1, d_data', 'LineWidth', 1)
xlim([0 this_tmax])
ylim([0 80]);
xlabel('t');
ylabel('Hamming distance');
set(gca, 'FontSize', 32);

qsave = 0;
fname_str = strrep(sprintf('Hamming_d_vs_t_a0_%.1f_N%d_K%d_Con%d_hill_%.1f_%s_%druns',...
    this_a0, N, K, Con, hill, initialID, nruns), '.', 'p');
fname = fullfile(save_fig_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot sigma(t) (overall)
%a0_idx = 1;
this_a0 = a0_all(a0_idx);
this_tmax = ceil(max(t_out_all(a0_idx, :))/100)*100;

sigma_data = squeeze(sigma_dynamics_all(a0_idx, :, 1:this_tmax, 1));
h=figure;
plot(0:this_tmax-1, sigma_data', 'LineWidth', 1)
ylim([0 1]);
xlabel('t');
ylabel('Variability \sigma(X_k)');
set(gca, 'FontSize', 32);

qsave = 1;
fname_str = strrep(sprintf('Sigma_cells_vs_t_a0_%.1f_N%d_K%d_Con%d_hill_%.1f_%s_%druns',...
    this_a0, N, K, Con, hill, initialID, nruns), '.', 'p');
fname = fullfile(save_fig_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%% Plot sigma(t) (ON/OFF)
%a0_idx = 1;
this_a0 = a0_all(a0_idx);
this_tmax = ceil(max(t_out_all(a0_idx, :))/100)*100;

sigma_data = squeeze(sigma_dynamics_all(a0_idx, :, 1:this_tmax, 2:3));
h=figure;
box on
hold on
p1 = plot(0:this_tmax-1, sigma_data(:,:,1)', 'r-', 'LineWidth', 1); %OFF cells
p2 = plot(0:this_tmax-1, sigma_data(:,:,2)', 'g-', 'LineWidth', 1); %ON cells
ylim([0 0.2]);
xlabel('t');
ylabel('Variability \sigma(X_k)');
set(gca, 'FontSize', 32);
legend([p1(1) p2(1)], {'OFF cells', 'ON cells'});

qsave = 0;
fname_str = strrep(sprintf('Sigma_ON_OFF_cells_vs_t_a0_%.1f_N%d_K%d_Con%d_hill_%.1f_%s_%druns',...
    this_a0, N, K, Con, hill, initialID, nruns), '.', 'p');
fname = fullfile(save_fig_folder, fname_str);
save_figure(h, 10, 8, fname, '.pdf', qsave);

%close all

end