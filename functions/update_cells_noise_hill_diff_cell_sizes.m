function [cells_out, changed] = ...
    update_cells_noise_hill_diff_cell_sizes(cells, dist, M_int, Con, K, a0, Rcell_all, noise, hill, prec)
% Update cells using noise in a positive feedback loop with finite hill
% coefficient
%%
% check input parameters
if (abs(M_int)~=1)
    error('M_int has the wrong form');
    return
end

% Calculate interaction terms
N = size(cells, 1);
M = ones(size(dist)); % Account for self-influence
for ii=1:N
    idx = setdiff(1:N, ii);
    M(ii, idx) = (Rcell_all(idx)./Rcell_all(ii)).*sinh( Rcell_all(ii) )...
        .*exp(Rcell_all(idx)-a0*dist(ii, idx)')./(a0*dist(ii, idx)');
end

% Concentration in each cell
C0 = 1 + (Con-1).*cells;

% Reading of each cell
Y = M*C0;

%dY = noise*sqrt(Y); % noise a la Berg-Purcell
%Y = Y + dY.*(2*rand(size(Y))-1);
dK = normrnd(0, noise, size(Y)); % Gaussian noise
K = K+dK;

if hill==Inf
    % binary cells
    cells_out = M_int*(Y - K) > 0;
    changed = ~isequal(cells_out, cells);
else
    % continuous cells
    if M_int==1
        cells_out = Y.^hill./(Y.^hill + K.^hill);
    elseif M_int==-1
        cells_out = K.^hill./(Y.^hill + K.^hill);
    end
    changed = ~isequal(round(cells_out, prec), round(cells, prec));
end