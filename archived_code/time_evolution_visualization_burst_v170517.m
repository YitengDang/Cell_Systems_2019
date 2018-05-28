% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON
close all
clear all
warning off

% Lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
% circuit parameters

K = 16;
Son = 8;
% initial conditions
p = 0.15;
iniON = round(p*N);

% burst parameters
B = [500]; %burst sizes
bt = [0]; %burst times, assume min(bt)>0
savemov = 1; %save movie? 0: no, 1: yes
subfolder = 'activation-deactivation region';

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
F = struct('cdata',[],'colormap',[]); % frames
% --Initial step to evolve the system--
% burst at t=0
this_bt=find(bt==t, 1);
this_B = 0;
if ~isempty(this_bt)
    this_B = B(this_bt); 
end
disp(strcat('B=',num2str(this_B)));
% initiate figure
update_cell_figure(hin, pos, a0, cells, cell_type, t);
F(1) = getframe();
% update cells
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
    % update figure
    update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
    F(end+1) = getframe();
    % update cells
    cells_hist{end+1} = cells_out;
    I(end+1) = moranI(cells_out, a0*dist);
    Non(end+1) = sum(cells_out);
    cells = cells_out;
    [cells_out, changed, mom(end+1)] = update_cells_burst(cells, dist, Son, K, a0, Rcell, this_B);
end

%% Create and export movie
pause(1);
close all
vname = sprintf('time_evolution_N%d_n%d_neq_%d_a0%d_K%d_Son%d_t%d_B%d', N, iniON, Non(end), 10*a0, K, Son, t, B);
% create versions to prevent duplicates
i=1;
vpath = fullfile(pwd, 'videos', 'burst', subfolder, strcat(vname, '-v',int2str(i)));
while exist(strcat(vpath,'.avi'), 'file') == 2
    i=i+1;
    vpath = fullfile(pwd, 'videos', 'burst', subfolder,  ...
        strcat(vname,'-v',int2str(i)));
end
if savemov
    myVideo = VideoWriter(vpath);
    myVideo.FrameRate = 1;  % Default 30
    myVideo.Quality = 100;    % Default 75
    open(myVideo);
    writeVideo(myVideo, F);
    close(myVideo);
    % Also save workspace in same directory (it contains all the frames for the
    % movies)
    save(vpath);
else
    fps=1; 
    movie(F,1,fps);
end