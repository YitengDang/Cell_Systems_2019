% Plot all parameters of a saved dynamics data

clear variables
close all
warning off

path = 'D:\eduardopavinat\Dropbox\Matlab codes\entropy\one cell type\data\dynamics_nonoise';
fname = 'N625_n63_neq_625_a05_K_14_Son_16_t_14-v1.mat';

% Save figure?
save_fig = 0;

load(fullfile(path,fname));
[~, name, ~] = fileparts(fname);

p = zeros(numel(cells_hist),1);
for i = 1:numel(cells_hist)
    p(i) = sum(cells_hist{i})/N;
end

h6 = figure(6);
hold on
piv = (0:N)/N;
Iv = -0.2:0.1:1;
E = @(p, I) -0.5*(Son-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Son+1)*(1+fN) - K);
u = @(p, I) 2*(Son-1)*fN*(2*p-1).*(1-I)+(Son+1)*(1+fN)-2*K;
v = @(p, I) 2*(Son-1)*fN*p.*(1-p);
[p_i,Imesh] = meshgrid(piv, Iv);
contour(piv, Iv, E(p_i, Imesh), 'ShowText', 'on')
[p_i,Im] = meshgrid(piv(1:4:end), Iv);
quiver(p_i, Im, u(p_i, Im), v(p_i, Im))
plot(p, I, 'k-o')
hold off

h7 = figure(7);
hold on
surf(p_i, Im, E(p_i, Im))
plot3(p, I, E(p, I'), 'Color', 'r')
hold off
