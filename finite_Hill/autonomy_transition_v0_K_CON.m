% Study the transition to autonomy for the finite Hill lattice in more
% detail
% Plot a map with K, C_{ON}
close all
clear all
warning off

% Lattice parameters
gridsize = 11;
N = gridsize^2;
initialID = 'uniform';

% circuit parameters
Klist = 5:15;
Conlist = 10:25;

hill = 2; % Hill coefficient
a0 = 6; % all a0 to compare
nruns = 50; % number of runs for each a0

% Variables to store
tfinal = zeros( numel(Klist), numel(Conlist), nruns ); % equilibration times
count_aut = zeros( numel(Klist), numel(Conlist) ); % number of autonomous lattices
cells_nonunif = {}; % final configurations which are non-uniform
%Noff = zeros( numel(Klist), numel(Conlist) , nruns); % number of `ON' cells

% Load data
path = fullfile(pwd, 'data', 'dynamics', 'finiteHill', '2017-09-08'); 
straux = '(\d+)';
straux2 = '(\w+)';
fpattern = strrep(sprintf('N%d_a0_%s_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_%s-v%s',...
            N, straux2, straux, straux, hill, straux, straux2, initialID, straux), '.', 'p');

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
%%
% Loop over all K, Con
for j=1:numel(Klist)
    K = Klist(j);
    for k=1:numel(Conlist)
    Con = Conlist(k);
    count2 = 0; % count number of loaded files
    fpattern = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_%s-v%s',...
            N, a0, K, Con, hill, straux, straux2, initialID, straux), '.', 'p');
        for i = 1:numel(names)
            % first get the filename given by the student
            [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
            if numel(tokens) > 0
                disp(names{i}) % displays the file name
                %disp(tokens{1});
                count2 = count2+1;
                % load the data
                load(fullfile(path,strcat(names{i},'.mat')));
                % store relevant data
                tfinal(j, k, count2) = t; % eq. time
                epsilon = 0.01; % threshold to compare final standard deviation with
                if xstd(t+1) > epsilon
                    cells_nonunif{end+1} = cells_hist{end}; % store non-unif. cells
                    count_aut(j, k) = count_aut(j, k) + 1;
                end
                %cells_binary = round(cells_hist{end}); % round to get ON/OFF cells
                %Noff(j, count2) = N-sum(cells_binary); 
            end
        end
    end
end

%% Plot map with fraction of non-uniform final states
h1=figure();
set(0, 'defaulttextinterpreter', 'latex');
frac = count_aut/nruns;
imagesc(Klist, Conlist, frac');
c = colorbar;
colormap('summer')
xlabel('$$K$$', 'FontSize', 24)
ylabel('$$C_{ON}$$', 'FontSize', 24)
ylabel(c, 'Fraction nonuniform lattices');
set(gca, 'ydir', 'normal','FontSize', 24)

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Frac_autonomy_N%d_a0_%.1f_K%.2fto%.2f_Con%.2fto%.2f_hill%.2f_%s',...
        N, a0, Klist(1), Klist(end), Conlist(1), Conlist(end), hill, initialID), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_KCon_map', fname_str);
    save_figure_pdf(h1, 10, 8, out_file);
end

%% Plot vs. Con/K
figure();
hold on
for i=1:numel(Conlist)
    for j=1:numel(Klist)
        plot(Conlist(i)/Klist(j), frac(j,i), 'x');
    end
end


%% Plot equilibration time
h2 = figure();

eqtime = mean(tfinal, 3);
imagesc(Klist, Conlist, eqtime');
c = colorbar;
colormap('hot');
xlabel('$$K$$','FontSize', 24);
ylabel('$$C_{ON}$$','FontSize', 24);
ylabel(c, '$$\langle t_{eq} \rangle$$', 'Interpreter', 'latex');
set(gca, 'ydir', 'normal', 'FontSize', 24)

qsave = 0;
if qsave
    fname_str = strrep(sprintf('Eq_time_avg_N%d_a0_%.1f_K%.2fto%.2f_Con%.2fto%.2f_hill%.2f_%s',...
        N, a0, Klist(1), Klist(end), Conlist(1), Conlist(end), hill, initialID), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_KCon_map', fname_str);
    save_figure_pdf(h2, 10, 8, out_file);
end

%% Plot lattice of a final non-uniform state
this_cells = cells_nonunif{651};
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

% Plot Xi values of final nonuniform state on a line
h5 = figure();
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

%% Save data in the same folder as figures
% generate map with parameter values at which nonuniform lattices were
% encountered
idx_pattern = find(count_aut>0 == 1); 
[I, J] = ind2sub([numel(Klist), numel(Conlist)], idx_pattern);
t = zeros(numel(Klist), numel(Conlist) );
for i=1:numel(I)
    t(I(i), J(i)) = 1;
    Klist_aut(i) = Klist( I(i) );
    Conlist_aut(i) = Conlist( J(i) );
end

figure();
imagesc(Klist, Conlist, t)
scatter(Klist_aut, Conlist_aut);

% save data
qsave = 0;
if qsave
    fname_str = strrep(sprintf('Frac_autonomy_N%d_a0_%.1f_K%.2fto%.2f_Con%.2fto%.2f_hill%.2f_%s',...
        N, a0, Klist(1), Klist(end), Conlist(1), Conlist(end), hill,...
        initialID), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'finite_Hill_autonomy_KCon_map',...
        strcat(fname_str, '.mat'));
    save(out_file, 'N', 'a0', 'Klist', 'Conlist', 'nruns', 'count_aut',...
        'cells_nonunif', 'Klist_aut', 'Conlist_aut');
end
