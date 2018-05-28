clear variables
close all

% Make the figures for spatial order examples

% Parameters of the system
gridsize = 20;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
save_fig = 1;
Non = 72;

cells = zeros(N, 1);

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = a0*dist_mat(pos,gridsize,gridsize,ex,ey);

% Stripe
cells(1:round(N/2)) = 1;

% Generate circle
% centerx = gridsize/2*a0 + gridsize/4*a0;
% centery = gridsize/4*a0*sqrt(3);
% 
% pos_mod = a0*pos;
% pos_mod(:,1) = pos_mod(:,1) - centerx;
% pos_mod(:,2) = pos_mod(:,2) - centery;
% 
% r = sqrt(sum(pos_mod.^2, 2));
% 
% cells(r<gridsize*a0/3) = 1;

% Chess board pattern
% for i = 1:gridsize
%     if mod(i,2) == 1
%         cells((i-1)*gridsize+1:2:i*gridsize) = 1;
%     else
%         cells((i-1)*gridsize+2:2:i*gridsize) = 1;
%     end
% end

% Random
% cells(randperm(N, N/2)) = 1;

I = moranI(cells, dist);
h1 = figure(1);
update_cell_figure(h1, pos, a0, cells, zeros(N,1), 0)
if save_fig > 0
    out_file = fullfile(pwd, 'figures', sprintf('moranI_N%d_I%d.pdf',N,round(100*I)));
    save_figure_pdf(h1, 10, 10, out_file);
end

