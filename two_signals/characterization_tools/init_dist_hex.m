function [dist, pos] = init_dist_hex(gridsizex, gridsizey)
% Try to load a previously calculated matrix. If not found calculate and
% save a new one initializing all cells in an hexagonal grid.

% Preiously saved matrices are saved in ./data/dist_matrix_hex with the
% filename <gridsizex><gridsizey>.mat
folder = fullfile('..', '..', 'data', 'dist_matrix_hex');
filename = fullfile(folder,...
    sprintf('%d%d.mat', gridsizex, gridsizey));

if exist(filename, 'file') == 2
    tmp = load(filename);
    dist = tmp.dist;
    pos = tmp.pos;
else
    [pos,ex,ey] = init_cellpos_hex(gridsizex,gridsizey);
    dist = dist_mat(pos,gridsizex,gridsizey,ex,ey);
    save(filename, 'pos', 'dist')
end