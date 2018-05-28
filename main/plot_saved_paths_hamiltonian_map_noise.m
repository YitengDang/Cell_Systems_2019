% Plot all parameters of a saved dynamics data for all the files saved that
% have the set of parameters established. This is the same as the file
% energy_figure_compare_saveddata.m but for files that run with noise

% This plots the lines in p vs I space for several saved runs on top of the
% hamiltonian map. It also plots the history of the hamiltonian
clear variables
close all
warning off

% Parameters of the system
gridsize = 21;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
pini = 0.8;
iniON = round(pini*N);

% parameters of the circuit
Son = 10;
K = 10;
noise = 0.1*K;

% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = 'H:\My Documents\Multicellular automaton\data\dynamics\noise';
straux = '(\d+)';
fpattern = sprintf('N%d_n%d_neq_%s_a0%d_K_%d_Son_%d_noise_%d_t_%s-v%s', ...
    N, iniON, straux, round(10*a0), round(K), round(Son), round(10*noise), straux, straux);

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
h1 = figure(1);
hold on
h2 = figure(2);
hold on
first = true;
pini = [];
pend = [];
Iini = [];
Iend = [];
hini = [];
hend = [];
tend = [];
mom_all = [];
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
            E = @(p, I) -0.5*(Son-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
                -(2*p-1).*(0.5*(Son+1)*(1+fN) - K);
            piv = (0:N)/N;
            Iv = -0.05:0.05:1;
            [p_i,Imesh] = meshgrid(piv, Iv);
            contourf(piv, Iv, E(p_i, Imesh),'LineStyle', 'none')
            colormap('summer')
            c = colorbar;
        end
        % Plot the line
        plot(p, I, 'Color', 'r')
        % Save the initial and end points for plotting the marks later
        pend(end+1) = p(end);
        Iend(end+1) = I(end);
        pini(end+1) = p(1);
        Iini(end+1) = I(1);
        hini(end+1) = -mom(1)/N;
        hend(end+1) = -mom(end)/N;
        tend(end+1) = t;
        mom_all(end+1:end+numel(mom)) = -mom/N;
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

% Set fonts and labels of the map
figure(h1)

set(gca,'FontSize', 20)
xlabel('p', 'FontSize', 24)
ylabel('I', 'FontSize', 24)
ylabel(c, 'h=H/N', 'FontSize', 24)
ylim([-0.05 1])

% Plot the maximum I estimated with first neighbour approx.
fa0 = sinh(Rcell)*exp(Rcell-a0)/a0;
G = 2*pi*exp(Rcell)*sinh(Rcell);
Acell = 0.5*sqrt(3)*a0^2;
Lon = @(p) real(sqrt((p*N - 7)*Acell/pi + 1.5*a0^2));
Loff = @(p) real(sqrt(((1-p)*N - 7)*Acell/pi + 1.5*a0^2));
Imax = @(p) (p.*(6*p*fa0 + G/Acell*(exp(-1.5*a0) - exp(-Lon(p)))) + ...
    (1-p).*(6*(1-p)*fa0 + G/Acell*(exp(-1.5*a0) - exp(-Loff(p)))) - ...
    fN*(2*p-1).^2)./p./(1-p)/4/fN;

plot(piv, Imax(piv), 'b')
hold off
%ylim([-0.05 1])

% save the map
name = sprintf('N%d_n%d_a0%d_K_%d_Son_%d_noise_%d_nruns_%d', ...
    N, iniON, round(10*a0), round(K), round(Son), round(10*noise), num_files);
out_file = fullfile(pwd, 'figures','hamiltonian_map','noise', strcat(name,'_h_path'));
save_figure_pdf(h1, 10, 6, out_file);

figure(h2)
hold off

set(gca,'FontSize', 20)
xlabel('Time (steps)', 'FontSize', 24)
ylabel('h', 'FontSize', 24)
%ylim([-0.05 0.3])
box on
name = sprintf('N%d_n%d_a0%d_K_%d_Son_%d_noise_%d_nruns_%d', ...
    N, iniON, 10*a0, K, Son, round(10*noise), num_files);
out_file = fullfile(pwd, 'figures','hamiltonian_map','noise', strcat(name,'_allH'));
save_figure_pdf(h2, 10, 6, out_file);

%% plot all energy differences
%{
h3 = figure(3);
scatter(mom_all(1:end-1), diff(mom_all))
xlabel('h_t');
%ylabel('h_{t+1}-h_t');
%}
