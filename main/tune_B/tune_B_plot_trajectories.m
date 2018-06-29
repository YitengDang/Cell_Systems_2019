% Plots individual simulated trajectories and saves work and energy data 

% --Use tune_B_analyze_work_distributions.m to analyze work data and
% genereate distributions
clear variables
close all
warning off
clc

% Parameters
gridsize = 11;
N = gridsize^2;
qsave = 1;
protocol = 6; % protocol ID
rev = 0; % reverse protocol?
if rev
    protid = strcat(num2str(protocol), 'R');
else
    protid = num2str(protocol);
end

% fix initial p
fixedp = 1;
pini = 0.2;
iniON = round(pini*N);

% Tuning parameters
% (1)---Constant a0, variable Son, K---
a0 = 0.5;
Rcell = 0.2*a0;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% --(1) Constant a0-------------
% protocol 1, 3
% 
Son_0 = 5; 
K_ini = (Son_0+1)/2*(1+fN);
tsteps = 20; % complete procedure in tsteps steps
K = linspace(K_ini, K_ini, tsteps+1);
if protocol == 3
    load('tuneB_protocol3_K4p5to4p5_Son5to11p5_B0to5.mat');
else
    Son = linspace(5, 15, tsteps+1);
end
if rev
    Son = Son(end:-1:1);
end
%}

% protocol 4% protocol 4 & 6
%
tsteps = 20;
K_ini = 12;
Son_ini = 2*K_ini/(1+fN) - 1;
%K = linspace(K_ini,2,tsteps+1); %protocol 4
K = linspace(K_ini,4,tsteps+1); %protocol 6
Son = linspace(Son_ini, Son_ini, tsteps+1);
if rev
    K = K(end:-1:1);
end
%}

%--- (2) constant Son, K, variable a0---
%{
tsteps = 10;
a0 = linspace(1.5, 0.5, tsteps+1);

[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
fN = zeros(1, tsteps+1);
for i=1:tsteps+1
    r = a0(i)*dist_vec(dist_vec>0); % exclude self influence
    Rcell = 0.2*a0(i);
    fN(i) = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
    %[out, map] = get_phenotype_map(a0, dist, Rcell, Son, K);
end
%
% fixed circuit parameters
Son = 5;
K = (Son+1)/2*(1+fN(1));
if rev
    fN = fN(end:-1:1);
    a0 = a0(end:-1:1);
end
%}
% calculate and show all B
B_all = (Son+1)/2*(1+fN) - K;
disp('B='); disp(B_all);
%}
%%
%----------------------------------------------
% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = strcat('H:\My Documents\Multicellular automaton\data\dynamics\B_field\protocol',protid); 
%path = strcat('H:\My Documents\Multicellular automaton\data\dynamics\B_field\protocol_test'); 
straux = '(\d+)';
if find(protocol == [1 3 4 5 6], 1)
    fpattern = strrep(sprintf('N%d_n%d_neq%s_a0_%.1f_K_%.fto%.f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d-v%s', ...
                N, iniON, straux, a0, K(1), K(end), Son(1), Son(end),...
                B_all(1), round(B_all(end),2), tsteps, straux), '.', 'p');
elseif protocol == 2
    fpattern = strrep(sprintf('N%d_n%s_neq%s_a0_%.1fto%.1f_K_%.1f_Son_%.f_B_%.2fto%.2f_tsteps%d-v%s', ...
            N, straux, straux, a0(1), a0(end), K, Son, B_all(1), B_all(end), tsteps, straux), '.', 'p');
else 
    disp('ERROR');
end
%
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
set(0, 'defaultTextInterpreter', 'latex'); 

% Plot trajectories in (p,I) plane
h1 = figure(1);
hold on
set(gca,'FontSize', 20)
xlabel('$$p$$', 'FontSize', 24)
ylabel('$$I$$', 'FontSize', 24)
xlim([0 1]);
ylim([-0.1 0.55]);

% Plot dq, dw against t
h2 = figure(2);
hold on
set(gca,'FontSize', 20)
xlabel('$$t$$', 'FontSize', 24)
ylabel('energy change', 'FontSize', 24)
xlim([0 tsteps+1]);
yyaxis right
p2a = plot((1:tsteps+1)-0.5, B_all, 'g', 'LineWidth', 2);
ylabel('$$B$$', 'FontSize', 24)

h3 = figure(3);
hold on
set(gca,'FontSize', 20)
xlabel('$$t$$', 'FontSize', 24)
ylabel('$$dh$$', 'FontSize', 24)
yyaxis right
plot((1:tsteps+1)-0.5, B_all, 'g', 'LineWidth', 2);
ylabel('$$B$$', 'FontSize', 24)

h4 = figure(4);
hold on
xlabel('$$dw+dq$$', 'FontSize', 24)
ylabel('$$dh$$', 'FontSize', 24)

%
w = [];
q = [];
delh = [];
pend = [];
Iend = [];
pfin = [];
Ifin = [];
pini = [];
Iini = [];
%test1 = []; 
%test2 = []; 
ntrials = 0;
for i = 1:numel(names) 
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0 %&& mod(i+5, 10)==0
        disp(names{i}) % displays the file name
        ntrials = ntrials + 1;
        % load the data
        load(fullfile(path,strcat(names{i},'.mat')));
        % Get the sequence of p
        p = Non/N;
        % plot trajectories in (p,I) space without hamiltonian map
        
        % Plot the line
        figure(h1)
        plot(p, I, 'Color', 'r')
        
        % Save initial and end points
        pend(end+1) = p(end);
        Iend(end+1) = I(end);
        pfin(end+1) = p(tsteps+1);
        Ifin(end+1) = I(tsteps+1);
        pini(end+1) = p(1);
        Iini(end+1) = I(1);
                       
        % Plot dw and dq trajectories, data from files
        figure(h2)
        yyaxis left
        p2b = plot(1:tsteps, dw, 'b-', 'LineWidth', 0.5);
        p2c = plot(1:tsteps, dq, 'r-', 'LineWidth', 0.5);
        w(end+1) = sum(dw);
        q(end+1) = sum(dq);
        
        % Plot dh trajectories, data from files
        figure(h3)
        dh = (H(2:end) - H(1:end-1))/N;
        yyaxis left
        plot(1:length(dh), dh, 'b-', 'LineWidth', 0.5);
        delh(end+1) = sum(dh);

        % Check dh == dw+dq;
        figure(h4);
        hold on
        plot(linspace(min(dw+dq), max(dw+dq), 100), linspace(min(dh), max(dh), 100), '--');
        scatter(dw+dq, dh(1:tsteps));
        
        % tests
        %test1(end+1) = all(dq<=0);
        %test2(end+1) = all(dh==dq+dw);
    end
end

% Plot the initial and final points
figure(h1)
scatter(pend, Iend, 'kx')
scatter(pini, Iini, 'ko')

% Add legend to figure 2
figure(h2)
legend([p2a p2b p2c], {'B', 'dw', 'dq'}, 'Location', 'eastoutside');
%%
% Save figures
% make file id for saving
if find(protocol == [1 3 4 5 6], 1)
    name = strrep(sprintf('N%d_a0_%.1f_K_%.1fto%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d_ntrials%d', ...
                N, a0, K(1), K(end), Son(1), Son(end), B_all(1), B_all(end),...
                tsteps, ntrials), '.', 'p');
elseif protocol == 2
    name = strrep(sprintf('N%d_a0_%.1fto%.1f_K_%.1f_Son_%.f_B_%.2fto%.2f_tsteps%d_ntrials%d', ...
            N, a0(1), a0(end), K, Son, B_all(1), B_all(end),...
            tsteps, ntrials), '.', 'p');
else
    disp('ERROR');
end

if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution',...
        strcat('protocol ', num2str(protocol) ),...
        strcat('trajectories_pI', name)); % filename
    save_figure_pdf(h1, 10, 8, out_file);
    
    out_file = fullfile(pwd, 'figures', 'work_distribution',...
        strcat('protocol ', num2str(protocol) ),...
        strcat('trajectories_dw_dq_', name)); % filename
    save_figure_pdf(h2, 12, 8, out_file);
    
    out_file = fullfile(pwd, 'figures', 'work_distribution',...
        strcat('protocol ', num2str(protocol) ),...
        strcat('trajectories_dh_', name)); % filename
    save_figure_pdf(h3, 10, 8, out_file);
end
%}
%% Save data
out_file = fullfile(pwd, 'data', 'work_distribution',...
    strcat('prot_', protid, '_', name, '.mat'));
save(out_file, 'q', 'w', 'delh');

%% Plot work against magnetic energy change Delta h_B
%{
DeltahB = - B_all(end)*(2*pend-1) + B_all(1)*(2*pini-1);
h3 = figure(3);
line = linspace(-4,2,100);
hold on
plot(line, line, 'LineWidth', 1.5); 
scatter(W, DeltahB, 100);
%xlim([-6 2]);
%ylim([-6 0]);
set(gca,'FontSize', 20)
xlabel('$$W$$','FontSize', 24);
ylabel('$$\Delta h_B$$', 'FontSize', 24);
legend('W = \Deltah_B', 'Location','nw');
if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution', ...
        strcat('W_vs_DeltahB_prot',protocol, name)); % filename
    save_figure_pdf(h3, 10, 8, out_file);
end
%}

%% Plot work distribution histogram
%--------Code moved to
%tune_B_analyze_work_distributions.m--------------------