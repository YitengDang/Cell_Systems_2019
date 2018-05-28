function [cells_out, changed, dY_mean, dY_std] = ...
    update_cells_noise_hill(cells, dist, Con, K, a0, Rcell, noise, hill, prec)
% Update cells using noise in a positive feedback loop with finite hill
% coefficient
%%
% Account for self-influence
idx = dist>0;
M = ones(size(dist));

% Matrix of cell reading
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

% Concentration in each cell
C0 = 1 + (Con-1).*cells;

% Reading of each cell
Y = M*C0;

%dY = noise*sqrt(Y); % noise a la Berg-Purcell
%Y = Y + dY.*(2*rand(size(Y))-1);
dY = normrnd(0, noise, size(Y)); % Gaussian noise
idx = find(Y+dY<0);
%counter = 0;
while ~isempty(idx)
    %disp(numel(idx));
    dY_new = normrnd(0, noise, size(idx));
    dY(idx) = dY_new;
    idx = find(Y+dY < 0);
    %counter = counter+1;
end
Y_final = Y+dY;
dY_mean = mean(dY);
dY_std = std(dY);
%%
cells_out = Y_final.^hill./(Y_final.^hill + K.^hill);
changed = ~isequal(round(cells_out, prec), round(cells, prec));



        