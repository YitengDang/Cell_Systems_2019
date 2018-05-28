% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON
% Saves movie frames for MATLAB; not used anymore now
% time_evolution_visualization_burst.m has its own movie function
close all
clear all
warning off

%data_path = 'D:\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';

% Lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;

% initial conditions
p = 0.1;
iniON = round(p*N);

% circuit parameters
Son = 8;
K = 16;

% burst parameters
B = 1000; %burst sizes
bt = 0; %burst times, assume min(bt)>0

% save movie
savemovie = 1; %1: save movie as collection of separate frames, 0: don't save
fname = sprintf('lattice_movie_Con_%d_K_%d_gz_%d_a0_%d_B_%d_frame', ...
Son, round(K), gridsize, 10*a0, B);

% Initialize parameters
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

% initialize ON cells
cells = zeros(N,1);
cells(randperm(N,iniON)) = 1;
%cells(1:2:N) = 1;

% initialize figure
hin = figure(1);
t = 0;
I = [];
Non = [];
mom = [];
cells_hist = {};
% --Initial step to evolve the system--
% burst at t=0
this_bt=find(bt==t, 1);
this_B = 0;
if ~isempty(this_bt)
    this_B = B(this_bt); 
end
disp(strcat('B=',num2str(this_B)));
%
if savemovie
    %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    frame=1;
    update_cell_figure_save(hin, pos, a0, cells, cell_type, t, B_all(1), fname, frame);
end
%
cells_hist{end+1} = cells;
Non(end+1) = sum(cells);
I(end+1) = moranI(cells, a0*dist);
[cells_out, changed, mom(end+1)] = update_cells_burst(cells, dist, Son, K, a0, Rcell, this_B);
% --Further steps to evolve the system--
while changed
    pause(1);
    % update time
    t = t+1;
    %k = waitforbuttonpress;
    % burst
    this_bt=find(bt==t, 1);
    this_B = 0;
    if ~isempty(this_bt)
        this_B = B(this_bt); 
    end
    disp(strcat('B=',num2str(this_B)));
    % update cells and save frame if instructed
    if savemovie
        %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        frame=t+1;
        update_cell_figure_save(hin, pos, a0, cells, cell_type, t, this_B, fname, frame);
    else
        update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
    end
    cells_hist{end+1} = cells_out;
    I(end+1) = moranI(cells_out, a0*dist);
    Non(end+1) = sum(cells_out);
    cells = cells_out;
    [cells_out, changed, mom(end+1)] = update_cells_burst(cells, dist, Son, K, a0, Rcell, this_B);
end

% Save workspace
%{
fname_str = sprintf('N%d_n%d_neq_%d_a0%d_K_%d_Son_%d_t_%d', ...
    N, iniON, Non(end), 10*a0, K, Son, t);
i = 1;
fname = fullfile(pwd, 'data','dynamics', 'no_noise', ...
    strcat(fname_str,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(pwd, 'data', 'dynamics', 'no_noise', ...
        strcat(fname_str,'-v',int2str(i),'.mat'));
end

save(fname)
%}