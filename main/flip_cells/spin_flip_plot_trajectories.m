%% Plots trajectories with almost similar starting conditions and their Hamming distance
% (1) Hamiltonian map
% (2) Hamming distance
% lattice parameters
clear all
close all
clc

gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
K = 10;
Con = 5;
p0 = 0.7;
iniON = round(p0*N);
qsave = 0;

% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = 'H:\My Documents\Multicellular automaton\data\dynamics\single_spin_flip'; 
straux = '(\d+)';
%fpattern = sprintf('N%d_n%d_neq_%s_a0%d_K_%d_Son_%d_t_%s_flip%s-v%s', ...
%    N, iniON, straux, 10*a0, K, Con, straux, straux, straux);
fpattern = sprintf('spin_flip_N%d_n%d_flip%s_neq_%s_a0%d_K_%d_Son_%d_t_%s-v%s', ...
            N, iniON, straux, straux, 10*a0, K, Con, straux, straux);

% Get all file names in the directory
listing = dir(path);
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
%% Load base trajectory (used for flipping)
for i = 1:numel(names)
    % first get the filename given by the student
    fpattern0 = sprintf('spin_flip_N%d_n%d_flip0_neq_%s_a0%d_K_%d_Son_%d_t_%s-v%s',...
        N, iniON, straux, 10*a0, K, Con, straux, straux);
    [tokens, ~] = regexp(names{i},fpattern0,'tokens','match');
    if numel(tokens)>0
        disp('yes');
        load(fullfile(path,strcat(names{i},'.mat')), 'cells_hist', 'fN');
        cells_base = cells_hist;
    end
end
%% Load single other trajectory and check Hamming distance
%{
for i = 1:numel(names)
    % first get the filename given by the student
    fpattern0 = sprintf('spin_flip_N%d_n%d_flip1_neq_%s_a0%d_K_%d_Son_%d_t_%s-v%s',...
        N, iniON, straux, 10*a0, K, Con, straux, straux);
    [tokens, ~] = regexp(names{i},fpattern0,'tokens','match');
    if numel(tokens)>0
        disp('yes');
        load(fullfile(path,strcat(names{i},'.mat')), 'cells_hist', 'fN');
    end
end
%%
% calculate distance from unflipped trajectory
t_this = numel(cells_hist);
t_base =  numel(cells_base);
t_max = max(t_this, t_base);
dist_temp = zeros(t_max, 1);
if t_this < t_base
    for i1=1:t_base-t_this
        cells_hist{end+1} = cells_hist{t_this};
    end
else 
    for i1=1:t_this-t_base
        cells_base{end+1} = cells_base{t_this};
    end
end
for i=1:t_max
    dist_temp(i) = sum(abs(cells_hist{i}-cells_base{i}));
end

figure();
plot(1:t_max, dist_temp)
%}
%% Load data and make plots
% Initialize all variables and figures
h1 = figure(1);
hold on
h2 = figure(2);
hold on
h3 = figure(3);
hold on
first = true;
pini = []; %zeros(numel(names),1);
pend = []; %pini;
Iini = []; %pini;
Iend = []; %pini;
hini = []; %pini;
hend = []; %pini;
tend = []; %pini;
hdist_tend = [];
hdist_end = [];
cells_final = {}; %cell(numel(names), 1);
dist_trajectories = {}; % Hamming distance with initial state as a function of time

for i = 1:numel(names)
    % first get the filename given by the student
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        % load the data
        load(fullfile(path,strcat(names{i},'.mat')));
        % store final cells
        cells_final{end+1} = cells_hist{end};        
        % Get the sequence of p
        p = zeros(numel(cells_hist),1);
        for k = 1:numel(cells_hist)
            p(k) = sum(cells_hist{k})/N;
        end
        figure(h1)
        % plot energy contour before plotting the lines
        if first
            first = false;
            E = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
                -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
            piv = (0:N)/N;
            Iv = -0.05:0.05:1;
            [p_i,Imesh] = meshgrid(piv, Iv);
            contourf(piv, Iv, E(p_i, Imesh),'LineStyle', 'none')
            colormap('summer')
            c = colorbar;
        end
        theta = 4*I'.*p.*(1-p) + (2*p-1).^2;
        % Plot the line
        plot(p, I, 'Color', 'r')
        % Save the initial and end points for plotting the marks later
        pend(end+1) = p(end);
        Iend(end+1) = I(end);
        pini(end+1) = p(1);
        Iini(end+1) = I(1);
        hini(end+1) = H(1)/N;
        hend(end+1) = H(end)/N;
        tend(end+1) = t;
        % plot the hamiltonian history
        figure(h2)
        plot(0:t, H/N, 'Color', 'r')
        
        % calculate distance from unflipped trajectory
        t_this = numel(cells_hist);
        t_base =  numel(cells_base);
        t_max = max(t_this, t_base);
        dist_temp = zeros(t_max, 1);
        if t_this < t_base
            for i1=1:t_base-t_this
                cells_hist{end+1} = cells_hist{t_this};
            end
        else 
            for i1=1:t_this-t_base
                cells_base{end+1} = cells_base{t_base};
            end
        end
        for i=1:t_max
            dist_temp(i) = sum(abs(cells_hist{i}-cells_base{i}));
        end
        % Store final variables
        hdist_tend(end+1) = t_max;
        hdist_end(end+1) = dist_temp(end);
        % Plot trajectories
        figure(h3);
        plot(1:t_max, dist_temp, 'Color', 'b');
        
    end
end

% Plot the initial and final points in both figures
figure(h1)
scatter(pend, Iend, 'kx')
scatter(pini, Iini, 'ko')

figure(h2)
scatter(tend, hend, 'kx')
scatter(zeros(size(hini)), hini, 'ko')

figure(h3)
scatter(hdist_tend, hdist_end, 'kx')

% Plot the maximum I estimated with first neighbour approx.
%{
figure(h1)
fa0 = sinh(Rcell)*exp(Rcell-a0)/a0;
Imax = @(p) (6 - 4./sqrt(p*N) - 6./(p*N) - 6*p)*fa0./(1-p)/fN;
Imax2 = @(p) (6*p - 4./sqrt((1-p)*N) - 6./((1-p)*N))*fa0./p/fN;
plot(piv, Imax(piv), 'b')
plot(piv, max(Imax(piv), Imax2(piv)), '--b', 'Linewidth', 1.5)
ylim([-0.05 1])
%}

%--Set fonts and labels of the map--
figure(h1)
hold off
set(gca,'FontSize', 20)
xlabel('p', 'FontSize', 24)
ylabel('I', 'FontSize', 24)
ylabel(c, 'h=H/N', 'FontSize', 24)
%ylim([-0.05 0.3])
set(gcf, 'Units', 'Inches', 'Position', [1 1 10 8]);

% save the map
qsave = 1;
if qsave
    name = sprintf('N%d_n%d_a0%d_K_%d_Con_%d', ...
        N, iniON, 10*a0, K, Con);
    out_file = fullfile(pwd, 'figures', 'single_spin_flip',... 
        strcat(name,'_hmap_pI')); %filename
    save_figure_pdf(h1, 10, 6, out_file);
    save_figure_eps(h1, 10, 6, out_file);
end

% Set fonts and labels of the hamiltonian graph
figure(h2)
hold off
set(gca,'FontSize', 20)
xlabel('Time (steps)', 'FontSize', 24)
ylabel('h', 'FontSize', 24)
%ylim([-0.05 0.3])
box on

% save the pdf
qsave = 1;
if qsave
    name = sprintf('N%d_n%d_a0%d_K_%d_Con_%d', ...
        N, iniON, 10*a0, K, Con);
    out_file = fullfile(pwd, 'figures', 'single_spin_flip',... 
        strcat(name,'_h_traj')); % filename
    save_figure_pdf(h2, 10, 8, out_file);
    save_figure_eps(h2, 10, 8, out_file);
end

% Set fonts and labels of the hamiltonian graph
figure(h3)
hold off
set(gca,'FontSize', 24)
xlabel('Time (steps)', 'FontSize', 24)
ylabel('$$d_H$$', 'FontSize', 24, 'Interpreter', 'latex')
%ylim([-0.05 0.3])
box on

% save the pdf
qsave = 1;
if qsave
    name = sprintf('N%d_n%d_a0%d_K_%d_Con_%d', ...
        N, iniON, 10*a0, K, Con);
    out_file = fullfile(pwd, 'figures', 'single_spin_flip',...
        strcat(name,'_hamming_dist_traj')); % filename
    save_figure_pdf(h3, 10, 8, out_file);
    save_figure_eps(h3, 10, 8, out_file);
end
%% Plot Hamming distance
hamming = zeros(N,1);
for i=1:N
    cells_noflip = cells_base{end};
    hamming(i) = sum(abs(cells_final{i+1}-cells_noflip));
end

h4 = figure(4);
histogram(hamming, 'Normalization', 'pdf');
xlim([0 max(hamming)+1]);
xlabel('$$d_H(X_{ini}, X_{eq})$$', 'Interpreter', 'latex');
ylabel('Probability');
set(gca, 'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [1 1 10 8]);

% save the pdf
qsave = 1;
if qsave
    name = sprintf('N%d_n%d_a0%d_K_%d_Con_%d', ...
        N, iniON, 10*a0, K, Con);
    out_file = fullfile(pwd, 'figures', 'single_spin_flip',... 'no_noise', ...
        strcat(name,'_eq_hamming_dist_hist')); % filename
    save_figure_pdf(h4, 10, 8, out_file);
    save_figure_eps(h4, 10, 8, out_file);
end
%% view configurations
%{
h1 = figure();
update_cell_figure(h1, pos, a0, cells_noflip, cell_type, t);

k = find(hamming == max(hamming))+1; %greatest difference
h2 = figure();
update_cell_figure(h2, pos, a0, cells_final{k}, cell_type, t);
%}
