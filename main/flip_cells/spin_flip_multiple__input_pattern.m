function [cells_flipped, cells_accepted] = ...
    spin_flip_multiple__input_pattern(cells_in, dist, Con, K, a0, Rcell, flip, nflips_in)
% Given a particular input pattern (cells), try flipping 'flip' number of
% spins and see if the flipped pattern returns to the initial pattern.

% Return list of flipped cells and an array (cells_accpeted) indicating
% which are the trials that gave these results.

N = size(dist,1); % number of cells
%nflips = 10; % number of random flips per lattice

% Test nflips random initial configurations and record the ones
% that lead to desired pattern
cells_flipped = cell(nflips_in, 1);
cells_accepted = zeros(nflips_in, 1);

% Limit trials for 0, 1 or 2 flips
% 1 flip: try all N possibilities
if flip<1 || ~(~mod(flip,1)) %non-positive or non-integer flip
    cells_flipped = {cells_in};
    cells_accepted = 1;
    return
elseif flip==1
    parfor idx = 1:N
        % flip cells
        cells = cells_in; 
        cells(idx) = ~cells(idx);
        cells_flipped{idx} = cells; %store initial config after flipping
        % update until equilibrium
        [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
        while changed
            [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
        end
        % compare cells
        if sum(abs(cells-cells_in))==0
            cells_accepted(idx) = 1;
        end
    end
    return 
end

% 2 flips: do not do more than nchoosek(N,2) flips
if flip==2
    nflips = min(nchoosek(N,2), nflips_in);
else
    nflips = nflips_in; 
end

parfor idx = 1:nflips
    % flip cells
    cells = cells_in; 
    idx_flip = randperm(N, flip);
    cells(idx_flip) = ~cells(idx_flip);
    cells_flipped{idx} = cells; %store initial config after flipping
    % update until equilibrium
    [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
    while changed
        [cells, changed] = update_cells(cells, dist, Con, K, a0, Rcell);
    end
    % compare cells
    if sum(abs(cells-cells_in))==0
        cells_accepted(idx) = 1;
    end
end


end
