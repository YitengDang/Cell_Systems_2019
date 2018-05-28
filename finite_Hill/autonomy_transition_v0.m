% Study the transition to autonomy for the finite Hill lattice in more
% detail
% At fixed K, CON, plot fraction of nonuniform lattices against a0
close all
clear all
warning off

% Lattice parameters
gridsize = 11;
N = gridsize^2;
initialID = 'uniform';
prec = 8; %precision

% circuit parameters
Klist = [8]; %[7 8 8 9 9 9 9 10 10 10 10 10 11 11 11 11 11 12 12 12 13];
Conlist = [16]; %[13 15 16 17 18 19 20 19 20 21 22 23 21 22 23 24 25 23 24 25 25];

%K = 7;
%Con = 13;
hill = 2; % Hill coefficient
all_a0 = 5.2:0.05:5.7; % all a0 to compare
nruns = 100; % number of runs for each a0

% Variables to store
tfinal = zeros(numel(all_a0), nruns); % equilibration times
count_aut = zeros(numel(all_a0), 1); % number of nonunif final states
cells_nonunif = {}; % final configurations which are non-uniform
Noff = zeros(numel(all_a0), nruns); % number of `ON' cells

% Load data
path = fullfile('H:\My Documents\Multicellular automaton', 'data', 'finiteHill', 'dynamics');
straux = '(\d+)'; % patterns for later
straux2 = '(\w+)';

% Get all file names in the directory
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;

for i_loop = 1:num_files
    filename = listing(i_loop+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
    end
end
%%
for i=1:numel(Klist) % loop over selected K, Con
    K=Klist(i);
    Con=Conlist(i);
    
    % Variables to store
    tfinal = zeros(numel(all_a0), nruns); % equilibration times
    count_aut = zeros(numel(all_a0), 1); % number of nonunif final states
    cells_nonunif = {}; % final configurations which are non-uniform
    Noff = zeros(numel(all_a0), nruns); % number of `ON' cells

% Loop over all a0
for j_loop=1:numel(all_a0)
    count2 = 0; % count number of loaded files
    a0 = all_a0(j_loop);
    fpattern = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_prec_%d_%s-v%s',...
            N, a0, K, Con, hill, straux, straux2, prec, initialID, straux), '.', 'p');
    for i_loop = 1:numel(names)
        % first get the filename given by the student
        [tokens, ~] = regexp(names{i_loop},fpattern,'tokens','match');
        if numel(tokens) > 0
            disp(names{i_loop}) % displays the file name
            %disp(tokens{1});
            count2 = count2+1;
            % load the data
            load(fullfile(path,strcat(names{i_loop},'.mat')));
            % store relevant data
            tfinal(j_loop, count2) = t; % eq. time
            epsilon = 0.01; % threshold to compare final standard deviation with
            if xstd(t+1) > epsilon
                cells_nonunif{end+1} = cells_hist{end}; % store non-unif. cells
                count_aut(j_loop) = count_aut(j_loop) + 1;
            end
            %cells_binary = round(cells_hist{end}); % round to get ON/OFF cells
            %Noff(j, count2) = N-sum(cells_binary); 
        end
    end
end

%%
% Fraction of autonomous lattices
frac = count_aut/nruns;
h1 = figure(1);
plot(all_a0, frac, 'o-', 'LineWidth', 1.5);
set(gca,'FontSize', 24);
xlabel('$$a_0$$', 'Interpreter', 'latex');
ylabel('Fraction non-uniform lattices');
xlim([a0list(1) a0list(end)]);
qsave = 1;
if qsave
    fname_str = strrep(sprintf('Transition_f_aut_N%d_a0_%.2fto%.2f_K%d_Con%d_hill%.2f_%s',...
        N, all_a0(1), all_a0(end), K, Con, hill, initialID), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_transition', fname_str);
    save_figure_pdf(h1, 10, 8, out_file);
    save_figure_eps(h1, 10, 8, out_file);

end
%%
%{
% Equilibration time
qsave = 0;

this_i = 2; %plot for this index of a0
%for this_i=1:numel(all_a0)
    close all
    fprintf('a0=%.2f\n', all_a0(this_i));

    h2 = figure();
    hold on
    histogram(tfinal(this_i, :), 'Normalization', 'pdf');
    [Pt, ti] = ksdensity(tfinal(this_i, :));
    plot(ti, Pt, 'LineWidth', 1.5);
    set(gca,'FontSize', 24);
    xlabel('$$t_{eq}$$', 'Interpreter', 'latex');
    ylabel('$$P(t_{eq})$$', 'Interpreter', 'latex');
    hold off
    
    if qsave
        fname_str = strrep(sprintf('Eq_time_N%d_a0_%.2f_K%d_Con%d_hill%.2f_%s',...
            N, all_a0(this_i), K, Con, hill, initialID), '.', 'p');
        out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_transition', fname_str);
        save_figure_pdf(h2, 10, 8, out_file);
    end
%end
%}
%% Mean equilibration times
qsave = 1;

h3=figure(3);
plot(all_a0, mean(tfinal, 2), 'o-', 'LineWidth', 1.5);
set(gca,'FontSize', 24);
xlabel('$$a_0$$', 'Interpreter', 'latex');
ylabel('$$\langle t_{eq} \rangle$$', 'Interpreter', 'latex');

if qsave
    fname_str = strrep(sprintf('Transition_eqtimes_N%d_a0_%.2fto%.2f_K%d_Con%d_hill%.2f_%s',...
        N, all_a0(1), all_a0(end), K, Con, hill, initialID), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_transition', fname_str);
    save_figure_pdf(h3, 10, 8, out_file);
    save_figure_eps(h3, 10, 8, out_file);
end
fprintf('K = %d, Con = %d \n', K, Con);
fprintf('max eq. time = %.2f, at a0 = %.1f \n', max(mean(tfinal, 2)), all_a0(max(mean(tfinal, 2)) == mean(tfinal, 2)));
end
%{
%% Number of ON/OFF cells
qsave = 0;

this_i = 3; %plot for this index of a0
%for this_i=1:numel(all_a0)
    close all;
    fprintf('a0=%.2f\n', all_a0(this_i)); 
    
    h4 = figure(4);
    hold on
    histogram(Noff(this_i, :)/N, (0:N)/N, 'Normalization', 'pdf');
    xlim([0 (max(max(Noff))+1)/N]);
    set(gca,'FontSize', 20);
    xlabel('$$f_{OFF}$$', 'Interpreter', 'latex');
    ylabel('$$P(f_{OFF})$$');
    
    if qsave
        fname_str = strrep(sprintf('F_OFF_N%d_a0_%.2f_K%d_Con%d_hill%.2f_%s',...
            N, all_a0(this_i), K, Con, hill, initialID), '.', 'p');
        out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_transition', fname_str);
        save_figure_pdf(h4, 10, 8, out_file);
    end
%end
%% Plot lattice of a final non-uniform state
this_cells = cells_nonunif{1};
hin = figure();
update_cell_figure_continuum(hin, pos, a0, this_cells, cell_type, t);

qsave = 0;
a0 = 0.5;
if qsave
    fname_str = strrep(sprintf('Final_lattice_N%d_a0_%.2f_K%d_Con%d_hill%.2f_%s',...
        N, a0, K, Con, hill, initialID), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_transition', fname_str);
    save_figure_pdf(hin, 10, 8, out_file);
end

%% Plot Xi values of final nonuniform state on a line
h5 = figure(5);
hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
             'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
             'XLim',[0 1],...               %#   set the x axis limit,
             'YLim',[0 eps],...             %#   set the y axis limit (tiny!),
             'Color','none');               %#   and don't use a background color
plot(this_cells(this_cells < 0.5),0,'r*','MarkerSize',10);  %# Plot data set 1
plot(this_cells(this_cells > 0.5),0,'b*','MarkerSize',10);  %# Plot data set 2
xlabel('$$X_i$$');
set(gca,'FontSize', 24);    

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Sample_lattice_Xi_N%d_a0_%.2f_K%d_Con%d_hill%.2f_%s',...
        N, a0, K, Con, hill, initialID), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_transition', fname_str);
    save_figure_pdf(h5, 10, 2, out_file);
end
%}