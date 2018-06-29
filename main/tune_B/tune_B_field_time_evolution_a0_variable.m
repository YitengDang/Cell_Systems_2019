% Time evolution of a system without noise, with changing B field
% Take a0 constant (saves calculating fN at each step)
% Takes only trajectories that start in autonomy

% Results to be analyzed with tune_B_plot_work_distribution.m
close all
clear all
clc
warning off
rng('shuffle');

% protocol
protocol = 2;
rev = 0;
if rev
    protid = strcat(num2str(protocol), 'R');
else
    protid = num2str(protocol);
end
% simulation parameters
nruns = 900;
qsave = 1;

% system size
gridsize = 11;
N = gridsize^2;

% variable a0, calculate fN
tsteps = 10;
a0 = linspace(1.5, 0.5, tsteps+1);
Rcellf = 0.2; %Rcell = Rcellf * a0
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
fN = zeros(1, tsteps+1);
for i=1:tsteps+1
    r = a0(i)*dist_vec(dist_vec>0); % exclude self influence
    Rcell = Rcellf*a0(i);
    fN(i) = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
    %[out, map] = get_phenotype_map(a0, dist, Rcell, Son, K);
end

% fixed circuit parameters
Son = 5;
K = (Son+1)/2*(1+fN(1));
if rev
    fN = fN(end:-1:1);
    a0 = a0(end:-1:1);
end

% list of B
B_all = (Son+1)/2*(1+fN)- K;
disp('B='); disp(B_all);
%%
for tests = 1:nruns
    disp(strcat('run ', int2str(tests)));
    
    % initialize cells
    cells = randi(2, N, 1)-1; % choose completely random configurations
    iniON = sum(cells);
    %cells = zeros(N,1);
    %iniON = 110;
    %cells(randperm(N, iniON)) = 1;
    
    % initialize variables
    t = 0;
    I = [];
    Non = [];
    dw = zeros(1, tsteps);
    dq = zeros(1, tsteps);
    dh = [];
    H = [];
    cells_hist = {};
    % store initial configuration
    cells_hist{end+1} = cells;
    Non(end+1) = sum(cells);
    I(end+1) = moranI(cells, a0(1)*dist);
    [~, ~, mom] = update_cells(cells, dist, Son, K, a0(1), Rcell);
    H(end+1) = -mom;
    % initial step: check whether cells are in equilibrium
    %[cells_out, changed, ~] = update_cells(cells, dist, Son_all(1), K_all(1), a0, Rcell);
    %if changed %take only configs that start in steady state
    %    disp('Initial state not stationary!');
    %    continue;
    %else
        % if in equilibrium, apply protocol to update cells
        for i=1:tsteps
            % update cells, calculate dw, dq
            [cells_out, changed, dq(i), dw(i), H(i+1)] =...
                update_cells_tune_B_field_varya0(cells, dist, a0(i+1), Rcellf, Son, K, a0(i));
            % update variables
            t = t+1;
            cells_hist{end+1} = cells_out;
            I(end+1) = moranI(cells_out, a0(i)*dist);
            Non(end+1) = sum(cells_out);      
            % calculate heat
            %dq(i) = -(2*B_all(i+1)+(Son(i+1)-1)*fN).*(Non(end) - Non(end-1))/N...
            %    -2*(Son(i+1)-1).*fN.*(Non(end)/N.*(1-Non(end)/N).*I(end)...
            %    - Non(end-1)/N.*(1-Non(end-1)/N).*I(end-1));
            cells = cells_out;
        end
    %end
    % if not in equilibrium after final step, continue until equilibrium is
    % reached
    while changed
        [cells_out, changed, mom] = update_cells(cells, dist, Son, K, a0(end), Rcellf*a0(end));
        %t = t+1;
        cells_hist{end+1} = cells_out;
        I(end+1) = moranI(cells_out, a0(end)*dist);
        Non(end+1) = sum(cells_out);
        H(end+1) = -mom;
        cells = cells_out;
    end
    
    % check dW and dQ agree with other calculations
    dh = (H(2:end) - H(1:end-1))/N;
    %disp('dq<=0?'); disp(all(dq<=0));
    %disp('dh='); disp(dq+dw);
    %disp('dq+dw==dh?'); disp(dh == dw+dq);
    
    % save file
    if qsave
        fname_str = strrep(sprintf('N%d_n%d_neq%d_a0_%.1fto%.1f_K_%.1f_Son_%.f_B_%.2fto%.2f_tsteps%d', ...
            N, iniON, Non(end), a0(1), a0(end), K, Son, B_all(1),...
            B_all(end), tsteps), '.', 'p');
        i = 1;
        fname = fullfile(pwd, 'data','dynamics', 'B_field', strcat('protocol', protid),...
            strcat(fname_str,'-v',int2str(i),'.mat'));
        while exist(fname, 'file') == 2
            i=i+1;
            fname = fullfile(pwd, 'data', 'dynamics', 'B_field', strcat('protocol', protid),...
                strcat(fname_str,'-v',int2str(i),'.mat'));
        end
        save(fname,'N', 'dw', 'dq', 'Non', 'I', 'H');
    end
end