function [out, map] = get_phenotype_map(a0, dist, Rcell, Son_vec, K_vec)
% Calculate the phenotype map limiting regions for a given set of
% parameters. This gives the regions: ON (R1), OFF (R3), and autonomy (R2
% if Son > K, R4 if Son < K).

dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

[K, Son] = meshgrid(K_vec, Son_vec);

% Get the signalling length
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere

% Make 4 limiting regions as boolean matrices
R1 = (K-1 < fN); % Everything ON
R2 = (Son > K & (K-1)./Son > fN); % autonomous cells for Son > K
R3 = (K./Son - 1 > fN); % Everything OFF
R4 = (Son <= K & K - Son < fN & (K-1)./Son > fN); % autonomous cells for Son < K

out = R1 + 2*R2 + 3*R3 + 4*R4; % regions do not overlap
% R1 -> red, R2 -> yellow, R3 -> magenta, R4 -> green, none -> white
map = [1, 1, 1
    1, 0, 0
    1, 1, 0
    1, 0, 1
    0, 1, 0];
map = map(1:max(max(out))+1,:);

