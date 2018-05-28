% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON
close all
clear all
warning off

% lattice parameters
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
% circuit parameters
Con = 15;
K = 6;
% initial conditions
p0 = 0.4;
iniON = round(p0*N);
I_min = -0.05;
I_max = I_min + 0.01;

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
[cells, I_ini, ~, test] = generate_I_new(cells, I_min, I_max, dist, a0);
%cells(1:2:N) = 1;

%%
% video for saving
frames = struct('cdata',[],'colormap',[]);

% initialize figure
hin = figure(1);
t = 0;
I = [];
Non = [];
h = [];
cells_hist = {};
%update_cell_figure(hin, pos, a0, cells, cell_type, t);
update_cell_figure_withI(hin, pos, dist, a0, cells, cell_type, t);
frames(1) = getframe(gcf);

% save vars and update cells
cells_hist{end+1} = cells;
Non(end+1) = sum(cells);
I(end+1) = moranI(cells, a0*dist);
[cells_out, changed, h(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
while changed
    pause(1);
    t = t+1;
    %k = waitforbuttonpress;
    %update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
    update_cell_figure_withI(hin, pos, dist, a0, cells, cell_type, t);
    frames(end+1) = getframe(gcf);
    cells_hist{end+1} = cells_out;
    I(end+1) = moranI(cells_out, a0*dist);
    Non(end+1) = sum(cells_out);
    cells = cells_out;
    [cells_out, changed, h(end+1)] = update_cells(cells, dist, Con, K, a0, Rcell);
end
%% Plot trajectory in p, I space
%{
figure();
hold on
plot(Non/N, I, 'r');
plot(Non(1)/N, I(1), 'ro');
plot(Non(end)/N, I(end), 'rx');
xlim([0 1]);
ylim([-0.05 1]);
%}
%% Plot h
%{
figure();
hold on
plot(0:t, h/N, 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('h = H/N');
set(gca,'FontSize', 24);
%}
%% Save result
%{
data_path = 'H:\My Documents\Multicellular automaton\videos';
fname_str = strrep(sprintf('binary_N%d_a0%.1f_Con%.2f_K%.2f_pini%.2f_Iini_%.2f_teq%d',...
    N, a0, Con, K, p0, I_ini, t), '.', 'p');
i = 1;
fname = fullfile(data_path, ...
    strcat(fname_str,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
  fname = fullfile(data_path, ...
      strcat(fname_str,'-v',int2str(i),'.mat'));
end
    
save(fname)
%}
%%
%{
figure(6);
movie(frames, 1, 2)
%}
%% Export movies
%{
fname_str = strrep(sprintf('binary_N%d_a0%.1f_Con%.2f_K%.2f_pini%.2f_Iini_%.2f_teq%d',...
    N, a0, Con, K, p0, I_ini, t), '.', 'p');
i = 1;
fname = fullfile('H:\My Documents\Multicellular automaton\videos',...
    strcat(fname_str,'-v',int2str(i),'.avi'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile('H:\My Documents\Multicellular automaton\videos',...
        strcat(fname_str,'-v',int2str(i),'.avi'));
end

myVideo = VideoWriter(fname, 'Uncompressed AVI');
myVideo.FrameRate = 5;  % Default 30
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);
%}