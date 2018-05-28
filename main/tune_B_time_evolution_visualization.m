% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON
close all
clear all
warning off

% fixed lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;

% fixed circuit parameters
%Son = 21;
%K = 6;

% initial conditions
%p0 = 0.3;
%iniON = round(p0*N);

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% variable lattice parameters
Son_0 = 5; 
K_0 = (Son_0+1)/2*(1+fN);

tsteps = 10; % complete procedure in tsteps steps
%Son_all=[5 5.24528812001061 5.53124349331728 5.86746718837490 6.26660291113958 6.74553754586907 7.32716227705940 8.04299585757536 8.93715627322097 10.0724797646909 11.5401266664718];
Son_all = linspace(5, 15, tsteps+1);
K_all = linspace(K_0, K_0, tsteps+1);
B_all = (Son_all+1)/2*(1+fN) - K_all;

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% save movie?
savemovie = 1;
%{
fname = strrep(sprintf('N%d_n%d_a0_%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d_frame', ...
            N, iniON, a0, Son_all(1), Son_all(end), B_all(1),...
            round(B_all(end),2), tsteps), '.', 'p');
%}
fname = strrep(sprintf('N%d_inirand_a0_%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d_frame', ...
    N, a0, Son_all(1), Son_all(end), B_all(1),...
    round(B_all(end),2), tsteps), '.', 'p');
%}
%%
% initialize ON cells
%cells = zeros(N,1);
%cells(randperm(N,iniON)) = 1;
cells = randi(2, N, 1) - 1;
%cells(1:2:N) = 1;

% initialize figure
hin = figure(1);
t = 0;
I = [];
Non = [];
dW = [];
cells_hist = {};
update_cell_figure(hin, pos, a0, cells, cell_type, t);
cells_hist{end+1} = cells;
Non(end+1) = sum(cells);
I(end+1) = moranI(cells, a0*dist);

if savemovie
    %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    frame=1;
    update_cell_figure_save(hin, pos, a0, cells, cell_type, t, B_all(1), fname, frame);
else
    update_cell_figure(hin, pos, a0, cells, cell_type, t);
end

% protocol
 for i=1:tsteps
    pause(1);
    [cells_out, changed, dW(i)] = update_cells_tune_B_field(...
        cells, dist, a0, Rcell, fN, Son_all(i+1), K_all(i+1), fN, Son_all(i), K_all(i));
    % update variables
    t = t+1;
    cells_hist{end+1} = cells_out;
    I(end+1) = moranI(cells_out, a0*dist);
    Non(end+1) = sum(cells_out);
    cells = cells_out;
    if savemovie
        %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        frame = frame +1;
        update_cell_figure_save(hin, pos, a0, cells, cell_type, t, B_all(i+1), fname, frame);
    else
        update_cell_figure(hin, pos, a0, cells, cell_type, t);
    end
end
%
while changed
    pause(1);
    [cells_out, changed, dW(i)] = update_cells_tune_B_field(...
        cells, dist, a0, Rcell, fN, Son_all(i), K_all(i), fN, Son_all(i), K_all(i));
    % update variables
    t = t+1;
    cells_hist{end+1} = cells_out;
    I(end+1) = moranI(cells_out, a0*dist);
    Non(end+1) = sum(cells_out);
    cells = cells_out;
    if savemovie
        %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        frame = frame +1;
        update_cell_figure_save(hin, pos, a0, cells, cell_type, t, B_all(end), fname, frame);
    else
        update_cell_figure(hin, pos, a0, cells, cell_type, t);
    end
end
%}

%% Save result
%{
fname_str = sprintf('n%d_neq_%d_a0%d_K_%d_Son_%d_t_%d', ...
    iniON, Non(end), 10*a0, K, Son, t);
i = 1;
fname = fullfile(data_path, 'dynamics_nonoise', ...
    strcat(fname_str,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
  fname = fullfile(data_path,'dynamics_nonoise', ...
      strcat(fname_str,'-v',int2str(i),'.mat'));
end
    
save(fname)
%}