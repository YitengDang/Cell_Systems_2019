function cond_set = trav_wave_conditions_nn(dist, Con, K, rcell, a0,...
    lambda, wave_type, nbands)
    % Calculates the conditions under which a travelling wave 
    % can propagate for a given input set, under the nearest-neighbour
    % approximation
    
    % Assume that the travelling wave consists of a single band. Diagonal
    % waves are not taken into account
    
    % input:
    % wave_type - type of travelling wave to consider, choose from "plane",
    % "inward", "outward" (latter two specify waves with one or more inward
    % and/or outward bends)
    
    % default settings:
    % z - number of nearest neighbours (default = 6) / coordination number
    gz = sqrt(size(dist, 1));
    Rcell = rcell*a0;
    z = 6;
    
    % background p depends on number of bands of the wave 
    if nargin<8
        nbands = 1; 
    end
    p = nbands.*[2/gz 2/gz]; 
    
    % ---------
    % Preliminary calculations
    fnn = zeros(1, 2);
    rnn = a0;
    fnn(1) = sinh(Rcell)*exp((Rcell-rnn)./lambda(1)).*(lambda(1)./rnn) ; % calculate signaling strength
    fnn(2) = sinh(Rcell)*exp((Rcell-rnn)./lambda(2)).*(lambda(2)./rnn) ; % calculate signaling strength
    
    idx = gz + round(gz/2); % pick cell not at corner of grid
    dist_vec = a0*dist(idx,:);
    r = dist_vec(dist_vec>0);
    fN = zeros(2,1);
    fN(1) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(1)).*(lambda(1)./r));
    fN(2) = sum(sinh(Rcell)*exp((Rcell-r)./lambda(2)).*(lambda(2)./r));
    
    Y_mf(1) = (fN(1) - z*fnn(1))*( Con(1)*p(1) + (1-p(1)) );
    Y_mf(2) = (fN(2) - z*fnn(2))*( Con(2)*p(2) + (1-p(2)) );
    
    % --------------
    % Default cell states for different wave types
    
    Xself = [0 0; 0 1; 1 1; 1 0; 0 0]; % fixed
    if strcmp(wave_type, 'plane')
        % plane wave
        Xnei = [0 2; 2 4; 4 4; 4 2; 2 0]; 
    elseif strcmp(wave_type, 'inward')
        % inward bend
        Xnei = [0 3; 3 5; 5 3; 3 1; 1 0]; 
    elseif strcmp(wave_type, 'outward')
        % outward bend
        Xnei = [0 1; 1 3; 3 5; 5 3; 3 0]; 
    end
    
    %%
    fnn_mat = repmat(fnn, 5, 1);
    Con_mat = repmat(Con, 5, 1);

    Y_self = (Con_mat-1).*Xself + 1;
    Y_nei = fnn_mat.*(Con_mat.*Xnei + (z-Xnei));
    Y_mf_mat = repmat(Y_mf, 5, 1);
    Y_all_mat = Y_self + Y_nei + Y_mf_mat;

    sgn1 = [-1; 1; 1; -1; -1];
    cond1 = (Y_all_mat(:,2) - K(1,2)).*sgn1 > 0;

    sgn2 = [-1 1; -1 1; 1 -1; 1 -1; 1 -1];
    %sgn2 = [-1 -1 1 1 1; ...
    %        1 1 -1 -1 -1];
    cond2a = (Y_all_mat(:,1) - K(2,1)).*sgn2(:, 1) > 0;
    cond2b = (Y_all_mat(:,2) - K(2,2)).*sgn2(:, 2) > 0;    
    cond2 = zeros(5, 1);
    cond2(1:2) = cond2a(1:2) & cond2b(1:2); % AND for first two conditions
    cond2(3:5) = cond2a(3:5) | cond2b(3:5); % OR for last two conditions

    cond_set = cond1 & cond2;
end