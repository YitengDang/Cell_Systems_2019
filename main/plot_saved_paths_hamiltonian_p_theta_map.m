% Plot all parameters of a saved dynamics data for all the files saved that
% have the set of parameters established.

% This plots the lines in p vs Theta space for several saved runs on top of the
% hamiltonian map. It also plots the history of the hamiltonian
% Finite Hill coefficient

clear variables
close all
warning off

% lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 5.4;
Rcell = 0.2*a0;
% circuit parameters
Con = 16;
K = 8;
hill = 2;
% initial conditions
initialID = 'uniform';
p0 = 0.3;
iniON = round(p0*N);
% save figure?
qsave = 1;

% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = 'H:\My Documents\Multicellular automaton\data\dynamics\finiteHill\2017-08-07'; 
straux = '(\d+)';
straux2 = '(\w+)';
fpattern = sprintf('N%d_a0_%s_K%d_Con%d_hill%s_t%s_xmeanf_%s_%s-v%s', ...
    N, '5p40', K, Con, '2p00', straux, straux2, ...
    initialID, straux);
%%
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
% Initialize all variables and figures
h1 = figure(1);
hold on
h2 = figure(2);
hold on
first = true;
pini = zeros(numel(names),1);
pend = pini;
th_ini = pini;
th_end = pini;
hini = pini;
hend = pini;
tend = pini;
for i = 1:numel(names)
    % first get the filename given by the student
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        % load the data
        load(fullfile(path,strcat(names{i},'.mat')), 'cells_hist', 'fN', 'I', 'mom', 't');
        % Get the sequence of p
        p = zeros(numel(cells_hist),1);
        for k = 1:numel(cells_hist)
            p(k) = sum(cells_hist{k})/N;
        end
        figure(h1)
        % plot energy contour before plotting the lines
        if first
            first = false;
            E = @(p, theta) -0.5*(Con-1)*(1+theta) - ((Con+1)/2*(1+fN) - K)*p;
            piv = (0:N)/N;
            thetav = -0.3:0.05:1;
            [p_i, Thetamesh] = meshgrid(piv, thetav);
            contourf(piv, thetav, E(p_i, Thetamesh),'LineStyle', 'none')
            colormap('summer')
            c = colorbar;
        end
        theta = 4*I'.*p.*(1-p) + (2*p-1).^2;
        % Plot the line
        plot(p, theta, 'Color', 'r')
        % Save the initial and end points for plotting the marks later
        pend(i) = p(end);
        th_end(i) = theta(end);
        pini(i) = p(1);
        th_ini(i) = theta(1);
        hini(i) = mom(1)/N;
        hend(i) = mom(end)/N;
        tend(i) = t;
        % plot the hamiltonian history
        figure(h2)
        plot(0:t, mom/N, 'Color', 'r')
    end
end

% Plot the initial and final points in both figures
figure(h1)
scatter(pend, th_end, 'kx')
scatter(pini, th_ini, 'ko')

figure(h2)
scatter(tend, hend, 'kx')
scatter(zeros(size(hini)), hini, 'ko')

% Plot the maximum I estimated with first neighbour approx.
%figure(h1)
%fa0 = sinh(Rcell)*exp(Rcell-a0)/a0;
%Imax = @(p) (6 - 4./sqrt(p*N) - 6./(p*N) - 6*p)*fa0./(1-p)/fN;
%Imax2 = @(p) (6*p - 4./sqrt((1-p)*N) - 6./((1-p)*N))*fa0./p/fN;
%plot(piv, Imax(piv), 'b')
%plot(piv, max(Imax(piv), Imax2(piv)), '--b', 'Linewidth', 1.5)
%ylim([-0.05 1])

% Set fonts and labels of the hamiltonian ap
figure(h1)
hold off
set(gca,'FontSize', 20)
xlabel('p', 'FontSize', 24)
ylabel('\Theta', 'FontSize', 24)
ylabel(c, 'h=H/N', 'FontSize', 24)
%ylim([-0.05 0.3])

% save the map
if qsave
name = strrep(sprintf('N%d_a0_%.2f_K_%d_Con_%d_hill_%.2f', ...
    N, a0, K, Con, hill), '.', 'p');
out_file = fullfile(pwd, 'figures', 'hamiltonian_map', 'finiteHill',... 
    strcat(name,'_hmap_p_theta')); %filename
save_figure_pdf(h1, 10, 6, out_file);
end

% Set fonts and labels of the hamiltonian trajectories
figure(h2)
hold off

set(gca,'FontSize', 20)
xlabel('Time (steps)', 'FontSize', 24)
ylabel('h', 'FontSize', 24)
%ylim([-0.05 0.3])
box on

% save the pdf
if qsave
name = strrep(sprintf('N%d_a0_%.2f_K_%d_Con_%d_hill_%.2f', ...
    N, a0, K, Con, hill), '.', 'p');
out_file = fullfile(pwd, 'figures', 'hamiltonian_map', 'finiteHill', ...
    strcat(name,'_h_traj')); % filename
save_figure_pdf(h2, 10, 6, out_file);
end