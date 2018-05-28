% Plot all parameters of a saved dynamics data for all the files saved that
% have the set of parameters established.

% This plots the lines in p vs I space for several saved runs on top of the
% hamiltonian map. It also plots the history of the hamiltonian
clear variables
close all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 5.2;
Rcell = 0.2*a0;
initialID = 'uniform';

% parameters of the circuit
K = 8;
Con = 16;
hill = 2;

% save?
qsave = 0;

% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = 'H:\My Documents\Multicellular automaton\data\dynamics\finiteHill\2017-08-07'; 
straux = '(\d+)';
straux2 = '(\w+)';
%fpattern = strrep(sprintf('N%d_a0%d_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_%s-v%s',...
%            N, 10*a0, K, Con, hill, straux, straux2, initialID, straux), '.','p');
fpattern = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_%s-v%s',...
            N, a0, K, Con, hill, straux, straux2, initialID, straux), '.','p');
        
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

% homogenization condition
cells_temp = fp(1)*ones(N,1);
cells_temp2 = fp(3)*ones(N,1);
cells_temp(1) = fp(3);
cells_temp2(1) = fp(1);

sigma_c_OFF = std(cells_temp);
sigma_c_ON = std(cells_temp2);

%hin = figure();
%update_cell_figure_continuum_with_I2(hin, pos, a0, cells_temp2, cell_type, 0, 0);


%% Initialize all variables and figures
h1 = figure(1);
hold on
h2 = figure(2);
hold on
first = true;
pini = zeros(numel(names),1);
pend = pini;
Iini = pini;
Iend = pini;
hini = pini;
hend = pini;
tend = pini;
nruns = 0;
for i = 1:numel(names)
    % first get the filename given by the student
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        nruns = nruns + 1;
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
        pend(i) = p(end);
        Iend(i) = I(end);
        pini(i) = p(1);
        Iini(i) = I(1);
        hini(i) = -mom(1)/N;
        hend(i) = -mom(end)/N;
        tend(i) = t;
        % plot the hamiltonian history
        figure(h2)
        plot(0:t, -mom/N, 'Color', 'r')
    end
end

% Plot the initial and final points in both figures
figure(h1)
scatter(pend, Iend, 'kx')
scatter(pini, Iini, 'ko')

figure(h2)
scatter(tend, hend, 'kx')
scatter(zeros(size(hini)), hini, 'ko')

% Plot the maximum I estimated with first neighbour approx.
%{
figure(h1)
fa0 = sinh(Rcell)*exp(Rcell-a0)/a0;
Imax = @(p) (6 - 4./sqrt(p*N) - 6./(p*N) - 6*p)*fa0./(1-p)/fN;
Imax2 = @(p) (6*p - 4./sqrt((1-p)*N) - 6./((1-p)*N))*fa0./p/fN;
%plot(piv, Imax(piv), 'b')
%plot(piv, max(Imax(piv), Imax2(piv)), '--b', 'Linewidth', 1.5)
ylim([-0.05 1])
%}

% Set fonts and labels of the map
figure(h1)
hold off

set(gca,'FontSize', 20)
xlabel('p', 'FontSize', 24)
ylabel('I', 'FontSize', 24)
ylabel(c, 'h=H/N', 'FontSize', 24)
%ylim([-0.05 0.3])

% save the map
if qsave
name = strrep(sprintf('N%d_a0%d_K%d_Con%d_hill%.2f_%s_runs%d',...
            N, 10*a0, K, Con, hill, initialID, nruns), '.','p');
out_file = fullfile(pwd, 'figures', 'hamiltonian_map', 'finiteHill',... 
    strcat(name,'_h_path')); %filename
save_figure_pdf(h1, 10, 6, out_file);
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
if qsave
name = strrep(sprintf('N%d_a0%d_K%d_Con%d_hill%.2f_%s_runs%d',...
            N, 10*a0, K, Con, hill, initialID, nruns), '.','p');
out_file = fullfile(pwd, 'figures', 'hamiltonian_map', 'finiteHill', ...
    strcat(name,'_allH')); % filename
save_figure_pdf(h2, 10, 6, out_file);
end

%% Plot distribution of final p or I
h3 = figure();
idx = (pend~=0); %throw away zero entries
histogram(pend(idx), 0:0.05:1, 'Normalization', 'count');
xlim([0 1]);
xlabel('$$p_{final}$$', 'FontSize', 24);
ylabel('Count', 'FontSize', 24);
set(gca,'FontSize', 20);

%
h4=figure();
histogram(Iend(idx), -0.2:0.05:1, 'Normalization', 'count');
xlim([-0.2 1]);
xlabel('$$I_{final}$$', 'FontSize', 24);
ylabel('Count', 'FontSize', 24);
set(gca,'FontSize', 20);

% Save figures
qsave = 1;
if qsave
    name = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_%s_runs%d',...
                N, a0, K, Con, hill, initialID, nruns), '.','p');
    out_file = fullfile(pwd, 'figures',...
    'finite_Hill_final_pI_distribution', strcat(name,'_p_distribution'));
    out_file2 = fullfile(pwd, 'figures',...
        'finite_Hill_final_pI_distribution', strcat(name,'_I_distribution'));
    save_figure_pdf(h3, 10, 6, out_file);
    save_figure_pdf(h4, 10, 6, out_file2);
end