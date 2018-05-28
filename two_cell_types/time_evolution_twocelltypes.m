% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON
close all
clear all
warning off

%--------------------------------------------------------------------------
gridsize = 15;
N = gridsize^2;
frac_type1 = 0.5;
N1 = round(frac_type1*N);
N2 = N - N1;
a0 = 1.5;

% Parameters 
Con = [35 5]; %column j = type j 
K = [5 5];
p = [0.5 0.5];
Rcell = [0.2 0.2]*a0;

Mcomm = [1 0; 0 1]; % Communication matrix, type i reacts to type j iff M_ij=1
Mcomm = logical(Mcomm);

% Simulation parameters
nruns = 100;
tmax = 1000;

%--------------------------------------------------------------------------
% Distance and position
[dist, pos] = init_dist_hex(gridsize, gridsize);

% cell types
cell_type = zeros(N,1);
idx1 = randperm(N,N1);
idx2 = setdiff(1:N, idx1);
cell_type(idx2) = 1;
% extra randomization of cell positions by requiring that cell_type also
% has a Moran's I value of around 0
[cell_type, ~] = generate_I_new(cell_type, 0, 0.01, dist, a0);
idx2 = find(cell_type);
idx1 = setdiff(1:N, idx2);

% Matrix of cell reading
M = zeros(size(dist));

% same type
M(idx1, idx1) = Mcomm(1,1)*sinh(Rcell(1))./(a0*dist(idx1, idx1)).*exp(Rcell(1)-a0*dist(idx1, idx1));
M(sub2ind(size(M), idx1, idx1)) = Mcomm(1,1)*1; % 1: normal, 0: no self communication

M(idx2, idx2) = Mcomm(2,2)*sinh(Rcell(2))./(a0*dist(idx2, idx2)).*exp(Rcell(2)-a0*dist(idx2, idx2));
M(sub2ind(size(M), idx2, idx2)) = Mcomm(2,2)*1; % 1: normal, 0: no self communication

% different type
M(idx1, idx2) = Mcomm(1,2)*Rcell(2)/(Rcell(1))*sinh(Rcell(1))./(a0*dist(idx1, idx2)).*exp(Rcell(2)-a0*dist(idx1, idx2)); % conc. cell of type 1 senses due to cells of type 2
M(idx2, idx1) = Mcomm(2,1)*Rcell(1)/(Rcell(2))*sinh(Rcell(2))./(a0*dist(idx2, idx1)).*exp(Rcell(1)-a0*dist(idx2, idx1));

% self-communication
%M(sub2ind(size(M), 1:N, 1:N)) = 1; % 1: normal, 0: no self communication

% Interaction strengths
f11 = sum(sum(M(idx1, idx1)))/N;
f12 = sum(sum(M(idx1, idx2)))/N; 
f21 = sum(sum(M(idx2, idx1)))/N; %not equal to f12 if cell radii different
f22 = sum(sum(M(idx2, idx2)))/N;

%dist_vec = a0*dist(1,:);
%r = dist_vec(dist_vec>0); % exclude self influence
%fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
%%
for run=1:nruns
    disp(run);
    % initialize ON cells
    cells = zeros(N,1);
    cells(idx1(randperm(N1,round(p(1)*N1)))) = 1;
    cells(idx2(randperm(N2,round(p(2)*N2)))) = 1;

    % initialize figure
    %hin = figure(1);
    t = 0;
    I = [];
    Non1 = [];
    Non2 = [];
    h = [];
    cells_hist = {};
    %update_cell_figure(hin, pos, a0, cells, cell_type, t);
    cells_hist{end+1} = cells;
    Non1(end+1) = sum(cells(idx1));
    Non2(end+1) = sum(cells(idx2));
    I(end+1) = moranI(cells, a0*dist);
    [cells_out, changed, h(end+1)] = update_cells_twocells(cells, idx1, idx2, M, Con, K);
    while changed && t<=tmax
        t = t+1;
        %k = waitforbuttonpress;
        %update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
        cells_hist{end+1} = cells_out;
        I(end+1) = moranI(cells_out, a0*dist);
        Non1(end+1) = sum(cells_out(idx1));
        Non2(end+1) = sum(cells_out(idx2));
        cells = cells_out;
        [cells_out, changed, h(end+1)] = update_cells_twocells(cells, idx1, idx2, M, Con, K);
    end

    % Save trajectory
    qsave = 1;
    if qsave
        path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes\temp';
        fname_str = strrep(...
            sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2f_Rcell2_%.2f_t%d', ...
            N1, N2, p(1), p(2), a0, K(1), K(2), Con(1), Con(2), Rcell(1), Rcell(2), t), '.', 'p');
        i = 1;
        
        fname = fullfile(path,...
            strcat(fname_str,'-v',int2str(i),'.mat'));
        while exist(fname, 'file') == 2
            i=i+1;
          fname = fullfile(path, ...
              strcat(fname_str,'-v',int2str(i),'.mat'));
        end
        save(fname); 
        %save(fname, 'N', 'N1', 'a0', 'Con', 'K', 'Rcell', 't', 'Non1', 'Non2', 'I', 'h');
    end
end