% Plots individual simulated trajectories and saves work and energy data 

% --Use tune_B_analyze_work_distributions.m to analyze work data and
% genereate distributions
clear variables
close all
warning off

% Parameters
gridsize = 11;
N = gridsize^2;
qsave = 0;
protocol = 4; % protocol ID
rev = 0; % reverse protocol?
if rev
    protid = strcat(num2str(protocol), 'rev');
else
    protid = num2str(protocol);
end

% Tuning parameters
% (1)---Constant a0, variable Son, K---
%{
a0 = 1.5;
Rcell = 0.2*a0;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

Son_0 = 5; 
K_ini = (Son_0+1)/2*(1+fN);
tsteps = 10; % complete procedure in tsteps steps
%Son=[5 5.24528812001061 5.53124349331728 5.86746718837490 6.26660291113958 6.74553754586907 7.32716227705940 8.04299585757536 8.93715627322097 10.0724797646909 11.5401266664718];
Son = linspace(5, 15, tsteps+1);
if rev
    Son = Son(end:-1:1);
end

K = linspace(K_ini, K_ini, tsteps+1);
B_all = (Son+1)/2*(1+fN) - K;
disp('B='); disp(B_all)
%}
%
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

% fixed circuit parameters
Son_ini = 5;
K = (Son_ini+1)/2*(1+fN(1));
if rev
    fN = fN(end:-1:1);
    a0 = a0(end:-1:1);
end
K = K*ones(1,tsteps+1);
Son = Son*ones(1,tsteps+1);

% list of B
B_all = (Son+1)/2*(1+fN)- K;
disp(B_all)
%}

% protocol 4
a0 = 1.5;
Rcell = 0.2*a0;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

tsteps = 10;
K_ini = 12;
Son_ini = 2*K_ini/(1+fN) - 1;
K = linspace(K_ini,2,tsteps+1);
Son = linspace(Son_ini, Son_ini, tsteps+1);
if rev
    K = K(end:-1:1);
end

B_all = (Son+1)/2*(1+fN) - K;
disp('B='); disp(B_all);
%%
%----------------------------------------------
% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = strcat('H:\My Documents\Multicellular automaton\data\dynamics\B_field\protocol',protid); 
%path = strcat('H:\My Documents\Multicellular automaton\data\dynamics\B_field\protocol_test'); 
straux = '(\d+)';
if find(protocol == [1 3], 1)
    fpattern = strrep(sprintf('N%d_n%s_neq%s_a0_%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d-v%s', ...
                N, straux, straux, a0, Son(1), Son(end), B_all(1), round(B_all(end),2), tsteps, straux), '.', 'p');
elseif protocol == 4
    fpattern = strrep(sprintf('N%d_n%s_neq%s_a0_%.1f_K_%.fto%.f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d-v%s', ...
                N, straux, straux, a0, K(1), K(end), Son(1), Son(end), B_all(1), round(B_all(end),2), tsteps, straux), '.', 'p');
elseif protocol == 2
    fpattern = strrep(sprintf('N%d_n%s_neq%s_a0_%.1fto%.1f_K_%.1f_Son_%.f_B_%.2fto%.2f_tsteps%d-v%s', ...
            N, straux, straux, a0(1), a0(end), K, Son, B_all(1),...
            B_all(end), tsteps, straux), '.', 'p');
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

% make file id for saving later
if find(protocol == [1 3 4], 1)
    name = strrep(sprintf('N%d_a0_%.1f_K_%.1fto%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d_ntrials%d', ...
                N, a0, K(1), K(end), Son(1), Son(end), B_all(1), B_all(end), tsteps, count), '.', 'p');
elseif protocol
    name = strrep(sprintf('N%d_a0_%.1fto%.1f_K_%.1f_Son_%.f_B_%.2fto%.2f_tsteps%d_ntrials%d', ...
            N, a0(1), a0(end), K, Son, B_all(1), B_all(end),...
            tsteps, count), '.', 'p');
else 
    disp('ERROR');
end
%%
% Plot trajectories in (p,I) plane
set(0, 'defaultTextInterpreter', 'latex'); 
h1 = figure(1);
hold on

% Plot B field change with dw trajectories
h2 = figure(2);
hold on
yyaxis right
plot((1:tsteps+1)-0.5, B_all, 'r', 'LineWidth', 2);
ylabel('$$B$$', 'FontSize', 24)

%
w = [];
delh = [];
pend = [];
Iend = [];
pfin = [];
Ifin = [];
pini = [];
Iini = [];
test1 = []; 
test2 = []; 
test3 = [];
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
        
        % Plot dw trajectories, data from files
        figure(h2)
        yyaxis left
        plot(1:tsteps, dw, 'b-', 'LineWidth', 0.5);
        w(end+1) = sum(dw);
        
        % Plot dW2 trajectories, directly calculate dW2
        figure(h2)
        yyaxis left
        dB = B_all(2:end) - B_all(1:end-1);
        dW2 = -dB.*(2*p(1:tsteps)-1);
        plot(1:tsteps, dW2, 'c-', 'LineWidth', 0.5);
        
        % Plot dH trajectories
        figure(h2)
        yyaxis left
        % only valid for protocol 1 or 1v2
        hvals = -0.5*(Son-1).*(1+4*fN.*p(1:tsteps+1).*(1-p(1:tsteps+1)).*I(1:tsteps+1)...
            + fN*((2*p(1:tsteps+1)-1).^2))...
                -(2*p(1:tsteps+1)-1).*(0.5*(Son+1)*(1+fN) - K);
        dh = hvals(2:end) - hvals(1:end-1);
        %plot(1:tsteps, dh, 'r-', 'LineWidth', 0.5);  

        % Plot dQ trajectories
        figure(h2)
        yyaxis left
        %plot(1:tsteps, dh - dw, 'g-', 'LineWidth', 0.5);
        
        % Calculate dQ explicitly
        dq = -(2*B_all(2:end)+(Son(2:end)-1)*fN).*(p(2:tsteps+1)-p(1:tsteps))...
            -2*(Son(2:end)-1).*fN.*(p(2:tsteps+1).*(1-p(2:tsteps+1)).*I(2:tsteps+1)...
            - p(1:tsteps).*(1-p(1:tsteps)).*I(1:tsteps));
        figure(h2)
        yyaxis left
        %plot(1:tsteps, dq, 'm-', 'LineWidth', 0.5);
        
        % delH
        delh(end+1) = sum(dh);
        % tests
        test1(end+1) = all(dq<=0);
        test2(end+1) = all(dq==dh-dW2);
        test3(end+1) = all(dw==dW2);
    end
end

% Plot the initial and final points
figure(h1)
scatter(pend, Iend, 'kx')
scatter(pini, Iini, 'ko')

% Set remaining fonts and labels 
figure(h1)
hold off
set(gca,'FontSize', 20)
xlabel('$$p$$', 'FontSize', 24)
ylabel('$$I$$', 'FontSize', 24)
xlim([0 1]);
ylim([-0.1 0.55]);

figure(h2)
hold off
yyaxis left
set(gca,'FontSize', 20)
xlabel('$$t$$', 'FontSize', 24)
ylabel('$$dw$$', 'FontSize', 24)
xlim([0 tsteps+1]);
%legend({'B', 'dw'}, 'Location', 'ne');
%
% Save figures
if qsave
    out_file = fullfile(pwd, 'figures', 'work_distribution', ...
        strcat(name,'_pI_trajectories')); % filename
    save_figure_pdf(h1, 10, 8, out_file);

    out_file = fullfile(pwd, 'figures', 'work_distribution', ...
        strcat(name,'_dW_trajectories')); % filename
    save_figure_pdf(h2, 10, 8, out_file);
end
%%
% Calculate initial and final Hamiltonian
hfin = -0.5*(Son(end)-1)*(1+4*fN(end).*pfin.*(1-pfin).*Ifin + fN(end)*(2*pfin-1).^2) ...
                -(2*pfin-1).*(0.5*(Son(end)+1)*(1+fN(end)) - K(end));
hend = -0.5*(Son(end)-1)*(1+4*fN(end).*pend.*(1-pend).*Iend + fN(end)*(2*pend-1).^2) ...
                -(2*pend-1).*(0.5*(Son(end)+1)*(1+fN(end)) - K(end));
hini = -0.5*(Son(1)-1)*(1+4*fN(1).*pini.*(1-pini).*Iini + fN(1)*(2*pini-1).^2) ...
                -(2*pini-1).*(0.5*(Son(1)+1)*(1+fN(1)) - K(1));
test4 = (delh == (hfin-hini));

%% Plot work against magnetic energy change Delta h_B
tune_B_plot_work_vs_deltahB(B_all, pini, pend, w, qsave, protocol, name);

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
%% Save data
out_file = fullfile(pwd, 'data', 'work_distributions', ...
    strcat('prot_', protocol, '_', name, '.mat'));
save(out_file, 'W', 'hend', 'hini', 'hfin');

%% Plot work distribution histogram
%--------Code moved to
%tune_B_analyze_work_distributions.m--------------------
%{
set(0, 'defaultTextInterpreter', 'latex'); 
%load(fullfile(pwd, 'data', 'work_distributions',...
%    'N121_n73_a0_1p5_Son_ini5_B_0to7p59_tsteps10.mat'));
h4 = figure(4);
hold on
nbins = 30;
%histogram(W, nbins, 'Normalization', 'pdf');
histfit(W, nbins, 'kernel')
set(gca,'FontSize', 20)
xlabel('$$W$$', 'FontSize', 24)
ylabel('$$P(W)$$', 'FontSize', 24)
%xlim([25 max(W)+1]);
%mW = mean(W);
%plot([mW mW], [0 150], 'g-','LineWidth',2);
%
% save the pdf
if qsave
    name = strrep(sprintf('N%d_a0_%.1f_K_%.1fto%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d_ntrials%d', ...
                N, a0, K(1), K(end), Son(1), Son(end), B_all(1), B_all(end), tsteps, count), '.', 'p');
    %name = strrep(sprintf('N%d_a0_%.1fto%.1f_K_%.1f_Son_%.f_B_%.2fto%.2f_tsteps%d_ntrials%d', ...
    %        N, a0(1), a0(end), K, Son, B_all(1), B_all(end),...
    %        tsteps, count), '.', 'p');
    out_file = fullfile(pwd, 'figures', 'work_distribution', ...
        strcat(name,'_hist_kernel')); % filename
    save_figure_pdf(h4, 10, 8, out_file);
end


%}

