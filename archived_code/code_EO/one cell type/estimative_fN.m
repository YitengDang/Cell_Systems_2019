% Compare fN with its estimative fN_est. The estimative is given using the
% integral approach and a continuous of cells
clear variables
close all

% Parameters of the system
gridsize = 21;
N = gridsize^2;
a0_vec = linspace(0.1, 4, 50);

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence

fN = zeros(size(a0_vec));
fN_est = zeros(size(a0_vec));
correction = zeros(size(a0_vec));

for i = 1:numel(a0_vec)
    a0 = a0_vec(i);
    Rcell = 0.2*a0;
    % Get the signalling length
    fN(i) = sum(Rcell*sum(exp(Rcell-a0*r)./(a0*r))); % calculate signaling strength sphere

    L0 = 1.5*a0;
    N0 = 7;
    correction(i) = 6*sinh(Rcell)/a0*exp(Rcell-a0);

    Acell = 0.5*sqrt(3)*a0^2;
    alpha = 2*pi*sinh(Rcell)*exp(Rcell)/Acell;
    L = sqrt(Acell*(N-N0)/pi + L0^2);

    fN_est(i) = correction(i) + alpha*(exp(-L0)-exp(-L));
end

h1 = figure(1);
plot(a0_vec, fN, 'k', 'LineWidth', 1.5)
hold on
plot(a0_vec, fN_est, '--r', 'LineWidth', 1.5)
hold off
legend({'Exact', 'Approx.'})

set(gca, 'Fontsize', 20)
xlabel('a_0', 'Fontsize', 24)
ylabel('f_N', 'Fontsize', 24)

outfile = fullfile(pwd,'figures', 'est_fN.pdf');
save_figure_pdf(h1, 7, 6, outfile);

h2 = figure(2);
plot(a0_vec, 100*abs(fN-fN_est)./fN, '', 'LineWidth', 1.5);

set(gca, 'Fontsize', 20)
xlabel('a_0', 'Fontsize', 24)
ylabel('Rel. Error f_N (%)', 'Fontsize', 24)

outfile = fullfile(pwd,'figures', 'est_fN_error.pdf');
save_figure_pdf(h2, 7, 6, outfile);

xlim([a0_vec(1) 1])

outfile = fullfile(pwd,'figures', 'est_fN_error_inset.pdf');
save_figure_pdf(h2, 7, 6, outfile);
