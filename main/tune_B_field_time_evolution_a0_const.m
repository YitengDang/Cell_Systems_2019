% Time evolution of a system without noise, with changing B field
% Take a0 constant (saves calculating fN at each step)
% Takes only trajectories that start in autonomy

% Results to be analyzed with ...
clear variables
close all
clc
warning off
rng('shuffle');

% Parameters
nruns = 400;
gridsize = 11;
N = gridsize^2;
qsave = 1;
protocol = 6;
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

% protocols 1, 1v2(3)
%{
tsteps = 20; % complete procedure in tsteps steps
Son_0 = 5; 
K_ini = (Son_0+1)/2*(1+fN);
if protocol == 3
    load('tuneB_protocol3_K4p5to4p5_Son5to11p5_B0to5.mat');
elseif find(protocol == [1 5], 1)
    Son = linspace(5, 15, tsteps+1);
end
if rev
    Son = Son(end:-1:1);
end
K = linspace(K_ini, K_ini, tsteps+1);

%}
% protocol 4 & 6
%

K_ini = 12;
Son_ini = 2*K_ini/(1+fN) - 1;
%protocol 4
%tsteps = 10; K = linspace(K_ini,2,tsteps+1); 
%protocol 6
tsteps = 20; K = linspace(K_ini,4,tsteps+1);
Son = linspace(Son_ini, Son_ini, tsteps+1);
if rev
    K = K(end:-1:1);
end

B_all = (Son+1)/2*(1+fN) - K;
disp('B='); disp(B_all)
%}

% calculate and show all B
B_all = (Son+1)/2*(1+fN) - K;
disp('B='); disp(B_all)
%%
for tests = 1:nruns
    disp(strcat('run ', int2str(tests)));
    
    % initialize cells
    if fixedp
        cells = zeros(N, 1);
        cells(randperm(N, iniON)) = 1;
    else
        cells = randi(2, N, 1)-1; % choose completely random configurations
        iniON = sum(cells);
    end
    
    % initialize variables
    t = 0;
    I = [];
    Non = [];
    H = [];
    %dwB = zeros(1,tsteps);
    %dwR = zeros(1,tsteps);
    dw = zeros(1, tsteps);
    dq = zeros(1,tsteps);
    cells_hist = {};
    % store initial configuration
    cells_hist{end+1} = cells;
    Non(end+1) = sum(cells);
    I(end+1) = moranI(cells, a0*dist);
    [~, ~, mom] = update_cells(cells, dist, Son(1), K(1), a0(1), 0.2*a0(1));
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
            [cells_out, changed, dq(i), dw(i), H(i+1)] = update_cells_tune_B_field(...
                cells, dist, a0, Rcell, Son(i+1), K(i+1), Son(i), K(i));
            % update variables
            t = t+1;
            cells_hist{end+1} = cells_out;
            I(end+1) = moranI(cells_out, a0*dist);
            Non(end+1) = sum(cells_out);
            % calculate heat
            %dq(i) = -(2*B_all(i+1)+(Son(i+1)-1)*fN).*(Non(end) - Non(end-1))/N...
            %    -2*(Son(i+1)-1).*fN.*(Non(end)/N.*(1-Non(end)/N).*I(end)...
            %    - Non(end-1)/N.*(1-Non(end-1)/N).*I(end-1));
            cells = cells_out;
        end
    %end
    
    % check dW and dQ agree with other calculations
    dh = (H(2:end) - H(1:end-1))/N;
    %disp('dq<=0?'); disp(all(dq<=0));
    %disp('dh='); disp(dq+dw);
    %disp('dq+dw==dh?'); disp(dh == dw+dq);

    % if not in equilibrium after final step, continue until equilibrium is
    % reached
    while changed
        [cells_out, changed, mom] = update_cells(cells, dist, Son(end), K(end), a0(end), Rcell);
        t = t+1;
        cells_hist{end+1} = cells_out;
        I(end+1) = moranI(cells_out, a0*dist);
        Non(end+1) = sum(cells_out);
        H(end+1) = -mom;
        cells = cells_out;
    end
    % save file
    if qsave
        fname_str = strrep(sprintf('N%d_n%d_neq%d_a0_%.1f_K_%.fto%.f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d', ...
            N, iniON, Non(end), a0, K(1), K(end), Son(1), Son(end), B_all(1),...
            round(B_all(end),2), tsteps), '.', 'p');
        i = 1;
        fname = fullfile(pwd, 'data','dynamics', 'B_field', strcat('protocol', protid),...
            strcat(fname_str,'-v',int2str(i),'.mat'));
        while exist(fname, 'file') == 2
            i=i+1;
            fname = fullfile(pwd, 'data', 'dynamics', 'B_field', strcat('protocol', protid),...
                strcat(fname_str,'-v',int2str(i),'.mat'));
        end
        save(fname, 'N', 'dw', 'dq', 'Non', 'I', 'H');
    end
end
%%
%{
figure();
hold on
%plot(1:tsteps, dh);
plot(1:tsteps, dq+dw);
plot(1:tsteps, dw);
plot(1:tsteps, dq);
legend('dh','dw','dq');
%}