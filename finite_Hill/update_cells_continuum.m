function [cells_out, changed, H] = ...
    update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill, prec)
%%
% Update cells without noise in a positive feedback loop with finite hill
% coefficient: continuum of cell states

% Account for self-influence
% M: Matrix of cell reading, eq. S10 p.S5
idx = dist>0;
M = ones(size(dist)); 
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

% Concentration in each cell
C0 = 1 + (Con-1).*cells; % term in brackets of Eq. S9 p.S5

% Reading of each cell
Y = M*C0; % Eq. S9 p.S5

% Continuous output for finite Hill coefficient
%S = 1 + (Con-1)*Y.^hill./(K^hill+Y.^hill);% new secretion rate
cells_out = Y.^hill./(K^hill+Y.^hill);% new cell state 
%cells_out = round(cells_out, 10);

% Energy of configuration (Multicellular Hamiltonian)
H = -sum((2*cells-1).*(Y-K));
%H = sum((cells_out-cells).^2);

% Use h(p,I)
%[It, ~] = moranI(cells, a0*dist);
%pt = mean(cells);
%H = -0.5*(Con-1)*(1 + 4*fN.*pt.*(1-pt).*It + fN*(2*pt-1).^2) ...
%    -(2*pt-1).*(0.5*(Con+1)*(1+fN) - K);
%%
% Configuration changed? Continuous state: round Xi to compare cells
%prec = 5; %precision / #decimals; standard = 5
changed = ~isequal(round(cells_out, prec), round(cells, prec));  
end