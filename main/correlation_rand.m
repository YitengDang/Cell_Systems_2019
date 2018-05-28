% Computes the correlation function <Xi Xj> as a function of distance
% Averages over random configurations sampled from a certain distribution
clear variables
close all

% Parameters of the system
gridsize = 8;
N = gridsize^2;
%a0 = 0.5;
%Rcell = 0.2*a0;

% use hexagonal lattice
[dist, ~] = init_dist_hex(gridsize, gridsize);
%dist_vec = a0*dist(1,:);
%r = dist_vec(dist_vec>0); % exclude self influence
%fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% fix fraction of ON cells
p = 0.8;
iniON = round(p*N);

% number of configurations to average over
ntrials = 10000;

%% 
distvals = unique(round(dist, 5)); %unique distances; need to round to neglect numerical errors
correlation = zeros(length(distvals), ntrials); % Correlation function
nterms = zeros(length(distvals), 1); % number of terms to average over for each value

% get indices
rows = {};
cols = {};
for i=1:length(distvals)
    [row, col] = find(round(dist, 5) == distvals(i)); 
    rows{end+1} = row;
    cols{end+1} = col;
    nterms(i) = numel(row);
end

%%
for j=1:ntrials
    % initialize cells
    % (1) ON/OFF cells
    cells = 2*randi(2,N,1)-3;
    
    % (2) ON/OFF at fixed p
    %cells = zeros(N, 1)-1; % OFF cells: -1
    %cells(randperm(N, iniON)) = 1; %ON cells: 1
    
    % (3) Xi in [-1, 1], uniform
    %cells = 2*rand(N, 1)-1;
    
    % (4) Xi in [-1, 1], normal(p, min(p,1-p)) on [0,1];
    %cells = 2*(min(p,1-p)*randn(N, 1)+p) - 1;
    %idx = (cells < -1); idx2 = (cells > 1);
    %cells(idx) = -1; cells(idx2) = 1;
    
    for i=1:length(distvals)
        row = rows{i}; col = cols{i};
        correlation(i, j) = sum(cells(row).*cells(col))/numel(row)-mean(cells)^2; 
    end
end


corr_func = mean(correlation, 2);
figure();
plot(distvals, corr_func, 'x-');
xlabel('distance');
ylabel('correlation');

figure();
plot(distvals, nterms/N, 'x-');
xlabel('distance');
ylabel('number of neighbours');
