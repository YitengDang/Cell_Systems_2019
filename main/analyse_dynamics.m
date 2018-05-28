% Analyze in detail the dynamics of cells in time for saved data. 
% Output: histograms of concentrations sensed by ON and OFF cells evolving over time.

% Legend: look up meaning of transition lines and label in legend

close all
clear variables

% Path where the data of the runs are
data_path = 'H:\My Documents\Multicellular automaton\data';

% Load a saved data from a run
[fname, path, ~] = uigetfile(fullfile(data_path));
load(fullfile(path,fname))
[~,name,~] = fileparts(fname);

% initialize cells in grid
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Matrix used to calculate the influence of cells on each other
M = zeros(size(dist));
M(dist > 0) = exp(Rcell - a0*dist(dist > 0))*sinh(Rcell)./(a0*dist(dist > 0));

% Calculate the parameters of the squared matrix to estimate the standard
% deviation of the concentration sensed by the cells
aux = M*M';
gN = aux(1,1);
xiN = sum(aux(1,2:end));

% Initialize variables
cells_ini = cells_hist{1}; % get the first state of the cells
p_ini = sum(cells_ini)/N; % get the initial fraction of ON cells
h = figure(1); % initialize the figure to plot

% For all the history of the cells
for i = 2:numel(cells_hist)
    cells = cells_hist{i}; % Get the cells state
    Y_self = (Son - 1)*cells + 1; % self sensing
    Y_nei = M*Y_self; % sensing from neighbors
    p = sum(cells)/N; % fraction of ON cells
    X = 2*cells - 1; % Change to +1 if ON -1 if OFF
    % Subtract the mean from the X vector and calculate I
    X = X - mean(X);
    I = (X'*M*X)/(X'*X)/fN;
    % Average concentration from neighbors for ON cells
    muon = fN*(Son*p + 1 - p + (Son-1)*(1-p)*I);
    % Average concentration from neighbors for OFF cells
    muoff = fN*(Son*p + 1 - p - (Son-1)*p*I);
    % Std of concentration from neighbors for ON & OFF cells
    sigmaon = sqrt((Son-1)^2*(gN)*p*(1-p));
    sigmaoff = sqrt((Son-1)^2*(gN)*p*(1-p));
    % Clear the figure and plot the histogram for ON and OFF cells
    clf(h,'reset');
    hold on
    histogram(Y_nei(cells == 0), 'Normalization', 'probability');
    histogram(Y_nei(cells == 1), 'Normalization', 'probability');
    xlabel('Sensed concentration');
    ylabel('Probability');
    %histogram(Y_nei, 'Normalization', 'probability');
    x = min(Y_nei):0.01:max(Y_nei);
    % Approximate as a normal distribution and compare
    pdfon = normpdf(x,muon,sigmaon);
    pdfoff = normpdf(x,muoff,sigmaoff);
    plot(x, pdfon, '-b')
    plot(x, pdfoff, '-r')
    % approximate as the composition of two normal distributions and
    % compare
    mup = fN*(Son*p_ini + 1 - p_ini);
    mu_test = 0.5*(Son*fN + K - 1);
    w1 = (muon - mu_test)/(mup-mu_test);
    %plot(x, (1-p)*pdfoff + p*pdfon, '--k')
    plot(x, w1*normpdf(x, mup,sigmaon) + ...
        (1-w1)*normpdf(x, mu_test,sigmaon),'--k')
    % Plot the transition lines for ON and OFF cells
    plot([K-1 K-1], [0 0.5], 'b')
    plot([K-Son K-Son], [0 0.5], 'k')
    plot([mup mup], [0 0.5], 'k')
    plot([mu_test mu_test], [0 0.5], 'g')
    % Calculate the real std for ON and OFF cells
    sigma = std(Y_nei);
    sigmaoff_calc = std(Y_nei(cells == 0));
    sigmaon_calc = std(Y_nei(cells == 1));
    % Legend
    legend('ON cells', 'OFF cells', 'ON estimate', 'OFF estimate', 'Composition', 'Location', 'bestoutside');
    hold off
    drawnow;
    pause(1) % wait cefore next time step
end