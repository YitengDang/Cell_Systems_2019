close all
clear all
warning off

% Calculate the pin_pout map for a fized conditions and plot the data

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
save_fig = 0;

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

K = 17;
Son = 7;

fname = sprintf('pin_pout_Con_%d_K_%d_gz_%d_a0_%d', ...
    Son, K, gridsize, 10*a0);

dist_vec = dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
kon = N*(K-fN-Son)/(Son-1)/fN;
koff = N*(K-fN-1)/(Son-1)/fN;

[count, t_av, nn_av] = count_eq(dist, Son, K, a0, Rcell);

p = (0:N)./N;

prob = transpose(count./repmat(sum(count,2),1,N+1));

h1 = figure(1);
im_fig = imagesc(p,p,prob);
title(sprintf('N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f', ...
    N, K, Son, a0, Rcell),'FontSize', 18)
set(gca,'Ydir','normal','FontSize', 20)
set(im_fig, 'AlphaData', count' > 0);
% colormap('copper');
% cmap = colormap;
% cmap = flipud(cmap);
% cmap(1,:) = [1 1 1];
% colormap(cmap)
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{in}', 'FontSize', 24)
ylabel('p_{out}', 'FontSize', 24)

% Organize and save
if save_fig > 0
    set(h1,'Units','Inches');
    set(h1, 'Position', [0 0 10 6 ])
    pos = get(h1,'Position');
    set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    out_file = fullfile(pwd, 'figures', strcat(fname,'_map'));
    print(h1, out_file,'-dpdf','-r0')
end

h2 = figure(2);
plot(p, t_av, 'r-o')
set(gca,'FontSize', 20)
title(sprintf('N = %d, K = %d, S_{ON} = %d, a0 = %.1f, R = %.1f', ...
    N, K, Son, a0, Rcell),'FontSize', 18)
xlabel('k_{in}', 'FontSize', 24)
ylabel('Average # steps for eq.', 'FontSize', 24)

h3 = figure(3);
% calculate the analytical formula
[~,omegak] = entropy_eq_sphere(dist_vec, Son, K, a0, Rcell);
plot((0:N)/N , log(omegak), 'LineWidth', 1.5)

h4 = figure(4);
im_fig = imagesc(p,p,nn_av');
set(gca,'ydir','normal')
set(im_fig, 'AlphaData', count' > 0);
colorbar