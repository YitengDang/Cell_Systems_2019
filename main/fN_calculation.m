% calculates the value for fN and gN given system parameters
function [fN, gN] = fN_calculation(gridsize, a0)
    % Parameters of the system
    %gridsize = 50;
    N = gridsize^2;
    %a0 = 0.5;
    Rcell = 0.2*a0;

    % use hexagonal lattice
    [dist, pos] = init_dist_hex(gridsize, gridsize);

    dist_vec = a0*dist(1,:);
    r = dist_vec(dist_vec>0); % exclude self influence
    fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
    gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2)); % calculate signaling strength
