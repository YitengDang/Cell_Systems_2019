clear variables

% Plot a line of the entropy map of a saved simulation and compare to the
% analytical formula (Con is changed along the line)

fname = fullfile(pwd,'data','20160518_0153_changeSon_N11K10_hexagonal.mat');
load(fname)

gridsize=11;
N = gridsize^2;
% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = dist(1,:);

omega = zeros(size(Son));
for i=1:size(Son)
    omega(i) = entropy_eq_sphere(dist_vec, Son(i), K, a0, Rcell);
end

% Plot the result
figure(1)
plot(Son,log(2^N*omega), 'r-')
hold on
plot(Son2,log(2^N*omega_sim), 'ko')
hold off
xlabel('S_{ON}')
ylabel('Entropy')
title(sprintf('K = %d, a_0 = %.1f', K, a0))