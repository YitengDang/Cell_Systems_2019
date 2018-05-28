% Plot all parameters of a saved dynamics data for all the files saved that
% have the set of parameters established.

% This plots the lines in p vs I space for several saved runs on top of the
% hamiltonian map. It also plots the history of the hamiltonian
clear variables
close all
warning off

%--------------------------------------------------------------------------
gridsize = 15;
N = gridsize^2;
frac_type1 = 0.5;
N1 = round(frac_type1*N);
N2 = N - N1;
a0 = 1.5;

% Parameters 
Con = [35 5]; %column j = type j 
K = [5 5];
p = [0.5 0.5];
Rcell = [0.2 0.2]*a0;

Mcomm = [1 0; 0 1]; % Communication matrix, type i reacts to type j iff M_ij=1
Mcomm = logical(Mcomm);
tmax = 100;

%--------------------------------------------------------------------------

% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes\temp';
%path = 'H:\My Documents\Multicellular automaton\two_cell_types\data\time_evolution';
straux = '(\d+)';

fpattern = strrep(sprintf(...
    'N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2f_Rcell2_%.2f_t%s-v%s', ...
    N1, N2, p(1), p(2), a0, K(1), K(2), Con(1), Con(2), Rcell(1), Rcell(2), straux, straux), '.', 'p');

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

% Initialize all variables and figures
h1 = figure(1); % trajectories
hold on
h2 = figure(2); % h map
hold on
first = true;
ii = 0;
for i = 1:numel(names)
    % first get the filename given by the student
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        ii = ii + 1;
        %figure(ii);
        %hold on;
        disp(names{i}) % displays the file name
        % load the data
        load(fullfile(path,strcat(names{i},'.mat')), ...
            'cells_hist', 'Non1', 'Non2', 'I', 'h');
        % Get the sequence of p
        p1 = Non1/N1;
        p2 = Non2/N2;
        p = (Non1+Non2)/N;
        %p1 = zeros(numel(cells_hist),1);
        %p2 = zeros(numel(cells_hist),1);
        %for k = 1:numel(cells_hist)
        %    p1(k) = sum(cells_hist{k}(idx1))/numel(idx1);
        %    p2(k) = sum(cells_hist{k}(idx2))/numel(idx2);
        %end
        figure(h1)
        % Plot the p trajectories
        plot(0:numel(cells_hist)-1, p1, '-r')
        plot(0:numel(cells_hist)-1, p2, '-b')
        % Plot trajectories in p,I space 
        figure(h2)
        plot(p, I, 'Color', 'r')
        % Save the initial and end points for plotting the marks later
        p1end(ii) = p1(end);
        p1ini(ii) = p1(1);
        p2end(ii) = p2(end);
        p2ini(ii) = p2(1);
        Iend(ii) = I(end);
        Iini(ii) = I(1);
        tend(ii) = numel(cells_hist)-1;

    end
end

% Set fonts and labels of the map
figure(h1)
scatter(tend, p1end, 'rx')
scatter(zeros(size(p1ini)), p1ini, 'ro') 
scatter(tend, p2end, 'bx')
scatter(zeros(size(p2ini)), p2ini, 'bo')
hold off
set(gca,'FontSize', 20)
xlabel('Time', 'FontSize', 24)
ylabel('p', 'FontSize', 24)

figure(h2)
scatter((p1end*N1+p2end*N2)/N, Iend, 'rx');
scatter((p1ini*N1+p2ini*N2)/N, Iini, 'ro');
set(gca,'FontSize', 20)
xlabel('p=p1+p2', 'FontSize', 24)
ylabel('I', 'FontSize', 24)
xlim([0 1]);
ylim([-0.2 1]);

% save the map
qsave = 0;
if qsave
    name = sprintf('N1-%d_n1-%d_N2-%d_n2-%d_a0-%d_K1-%d_K2-%d_Son1-%d_Son2-%d', ...
        N1, n1, N2, n2, round(10*a0), round(K1), round(K2), round(Son1), round(Son2));
    out_file = fullfile( ...
        pwd,'figures', 'twocelltypes', strcat(name,'_dynamics_p')); %filename
    save_figure_pdf(h1, 10, 6, out_file);
end
