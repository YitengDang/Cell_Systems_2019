% Plot all parameters of a saved dynamics data for all the files saved that
% have the set of parameters established.

% This plots the lines in p vs I space for several saved runs on top of the
% hamiltonian map. It also plots the history of the hamiltonian

clear variables
close all
warning off

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
% circuit parameters
Con = 11;
K = 18;
% initial conditions
p0 = 0.5;
iniON = round(p0*N);

% save figure?
qsave = 0;

% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = 'H:\My Documents\Multicellular automaton\data\dynamics\no_noise'; 
%path = 'H:\My Documents\Multicellular automaton\rebuttal\trajectories\N121_a0_1p5_figS4\K3_Con24';
%path = 'H:\My Documents\Multicellular automaton\data\dynamics\single_spin_flip';

straux = '(\d+)';
fpattern = sprintf('N%d_n%d_neq_%s_a0%d_K_%d_Son_%d_t_%s-v%s', ...
    N, iniON, straux, 10*a0, K, Con, straux, straux);
%fpattern = sprintf('N%d_n%s_neq_%s_a0%d_K_%d_Son_%d_t_%s-v%s', ...
%    N, straux, straux, 10*a0, K, Con, straux, straux);
%fpattern = sprintf('spin_flip_N%d_n%d_flip%s_neq_%s_a0%d_K_%d_Son_%d_t_%s-v%s',...
%    N, iniON, straux, straux, 10*a0, K, Con, straux, straux);

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

% Only for plotting initial configuration
%[dist, pos] = init_dist_hex(gridsize, gridsize);
%dist_vec = a0*dist(1,:);
%r = dist_vec(dist_vec>0); % exclude self influence
%fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
%cell_type = zeros(N,1); % all the same here

% Initialize all variables and figures
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
for i = 1:numel(names)
    % first get the filename given by the student
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        % load the data
        load(fullfile(path,strcat(names{i},'.mat')), 'cells_hist', 'fN', 'I', 'h', 't');
        
        % Plot initial configuration
        %if mod(i, 10)==0
        %    cells = cells_hist{1};
        %    htemp = figure(3);
        %    update_cell_figure(htemp, pos, a0, cells, cell_type, 0)
        %end
        
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
            Iv = -0.15:0.05:1;
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
        hini(i) = h(1);
        hend(i) = h(end);
        tend(i) = t;
        % plot the hamiltonian history
        figure(h2)
        plot(0:t, h, 'Color', 'r')
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
figure(h1)
fa0 = sinh(Rcell)*exp(Rcell-a0)/a0;
Imax = @(p) (6 - 4./sqrt(p*N) - 6./(p*N) - 6*p)*fa0./(1-p)/fN;
Imax2 = @(p) (6*p - 4./sqrt((1-p)*N) - 6./((1-p)*N))*fa0./p/fN;
%plot(piv, Imax(piv), 'b')
%plot(piv, max(Imax(piv), Imax2(piv)), '--b', 'Linewidth', 1.5)
ylim([-0.05 1])

% Set fonts and labels of the map
figure(h1)
hold off

set(gca,'FontSize', 20)
xlabel('p', 'FontSize', 24)
ylabel('I', 'FontSize', 24)
ylabel(c, 'h=H/N', 'FontSize', 24)
%ylim([-0.05 0.3])
%%
% save the map
if qsave
%name = sprintf('N%d_n%d_a0%d_K_%d_Con_%d', ...
%    N, iniON, 10*a0, K, Con);
name = sprintf('N%d_all_n_a0%d_K_%d_Con_%d',...
    N, 10*a0, K, Con);
%out_file = fullfile(pwd, 'figures', 'hamiltonian_map',... 'no_noise',... 
%    strcat(name,'_hmap_pI')); %filename
out_file = fullfile(pwd, 'rebuttal', 'trajectories', 'h_maps',... 'no_noise',... 
    strcat(name,'_hmap_pI')); %filename
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
%if qsave
%name = sprintf('N%d_n%d_a0%d_K_%d_Con_%d', ...
%    N, iniON, 10*a0, K, Con);
%out_file = fullfile(pwd, 'figures', 'hamiltonian_map',... 'no_noise', ...
 %   strcat(name,'_h_traj')); % filename
%save_figure_pdf(h2, 10, 6, out_file);
%end