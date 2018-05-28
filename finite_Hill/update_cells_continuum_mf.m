function [cells_out, changed, H] = ...
    update_cells_continuum_mf(cells, Con, K, hill, fN)

% Update cells without noise in a positive feedback loop with finite hill
% coefficient: continuum of cell states
% mf: mean-field approach

Xm = mean(cells);
Yself = cells.*(Con-1)+1;
Ymf = fN*(Xm*(Con-1)+1);
Y = Yself + Ymf;

% Continuous output for finite Hill coefficient
%S = 1 + (Con-1)*Y.^hill./(K^hill+Y.^hill);% new secretion rate
cells_out = Y.^hill./(K^hill+Y.^hill);% new cell state 
cells_out = round(cells_out, 6);

% Energy of configuration (Multicellular Hamiltonian)
%H = -sum((2*cells-1).*(Y-K));
H = sum((cells_out-cells).^2);

% Configuration changed? Continuous state: round Xi to compare cells
prec = 5; %precision / #decimals
changed = ~isequal(round(cells_out, prec), round(cells, prec)); 