function [pos_new, dist_new, rejections] = update_cell_positions_diff_cell_sizes(...
    n, rcell_all, pos, dist, sigma_D)
% Updates the positions of the cells according to Brownian motion
% sigma_D: width of the Gaussian for Brownian motion in units of a0
% * Move the cells after cell states have been updated
% * Move all cells simultaneously 

% Idea: can add different Brownian motion for cells of different sizes
% (e.g. Stokes-Einstein viscosity).

%------pre-lim code (for running function separately)----------
%sigma_D = 0.01;

% Start with perfect hexagonal lattice
%{
n = gz;
Lx = 1;
delx = Lx/n;
dely = sqrt(3)/2*delx;
Ly = dely*n;

[xm, ym] = meshgrid(0:n-1, 0:n-1);
x = (xm+mod(ym,2)/2)*delx;
y = ym*dely;
% store initial config
%x_ini = x;
%y_ini = y;

% Calculate distance
dist = calc_dist_periodic(x, y, Lx, Ly);
dist = round(dist, 10);
%pos = [x(:) y(:)];
%}

% Calculate default fN
%{
fN = zeros(N, 1);
for i=1:N
    dist_vec = dist(i,:);
    r = dist_vec(dist_vec>0); % exclude self influence
    fN(i) = sum(sinh(R)*sum(exp(R-r)./r)); % calculate signaling strength
end
%}
%{
fN0 = mean(fN);
if ~nodisplay
    fprintf('fN0 (default fN) = %.2f \n', fN0);
end
%}
%---------main code---------------
if sigma_D==0 % if cells are not moved, return
    pos_new = pos;
    dist_new = dist;
    rejections = 0;
    return
end
Lx = 1;
N = n^2;
R_all = rcell_all*Lx/n; % convert to box unit distances
max_rejections = 10^6;

% hexagonal placement
delx = Lx/n;
dely = sqrt(3)/2*delx;
Ly = dely*n;

% set initial position
x = pos(:, 1);
y = pos(:, 2);

% Express sigma_D in terms of box units
sigma_D_box_units = sigma_D*(Lx/n);

% sum of radii of pairs of cells, for checking overlap
rcell_mat = repmat(R_all, 1, N) + repmat(R_all', N, 1); 
%%
% Key parameters
cells_updated = []; % number of cells that have been updated
cells_to_update = 1:N;
rejections = 0;
while numel(cells_updated) < N
    %%
    if mod(rejections, 10^4)==0
        fprintf('Rejections: %d \n', rejections);
        fprintf('Cells to update: %d \n', numel(cells_to_update) );
    end
    %disp(step)
    %fprintf('Cells to do = %d \n', N-numel(cells_updated) );
    %fprintf('Rejections: %d \n', rejections);
    
    cells_to_update = setdiff(1:N, cells_updated);
    
    % choose random cell to update
    idx = randperm(numel(cells_to_update), 1);
    cell_i = cells_to_update(idx);
    
    % update cell position (tentative)
    x_new = x; 
    y_new = y;
    x_new(cell_i) = mod(x_new(cell_i) + sigma_D_box_units*randn(), Lx);
    y_new(cell_i) = mod(y_new(cell_i) + sigma_D_box_units*randn(), Ly);
    
    % Calculate distance
    dist_new = calc_dist_periodic(x_new, y_new, Lx, Ly);
    %%
    cond_mat = (dist_new > rcell_mat);
    if all(cond_mat(dist_new>0))
        %all(dist_new(cell_i, setdiff(1:N, cell_i) ) >= 2*R) % check only updated cell, exclude distance to self
        x = x_new;
        y = y_new;
        dist_new = round(dist_new, 10);

        % update list of updated cells
        cells_updated(end+1) = cell_i; 
    else
        rejections = rejections + 1;
        %fprintf('Rejected! \n');
        if rejections>=max_rejections
            warning('Too many rejections! Not all cell positions updated');
            pos_new = [x(:) y(:)];
            dist_new = dist_new/(Lx/n); 
            return
        end
    end
end
pos_new = [x(:) y(:)];
% adjust distances to scale so that nearest neighbour distances are
% 1 (in case of a perfect lattice); without normalization this
% would be Lx/n.
dist_new = dist_new/(Lx/n); 

% Draw lines between old and new positions
%{
% How to erase lines?
hold on
for ii=1:N
    plot([pos(ii,1) pos_new(ii,1)], [pos(ii,2) pos(ii,2)], '--');
end
hold off
%}
%-----------end main code----------------------
%% Draw configuration
%{
hin = figure();
clf(hin,'reset');
title(sprintf('N = %d, R = %.2f', N, Rcell), ...
    'FontSize', 24);
%set(gca,'YTick',[],'XTick',[]);
set(gca,'DataAspectRatio', [1 1 1]);

x = pos_new(:,1);
y = pos_new(:,2);

axis([0 Lx 0 Ly]);
box on
hold on
for i=1:N
    position = [x(i)-R y(i)-R 2*R 2*R];
    face_clr = 'k';
    curv = [1 1];
    rectangle('Position', position, 'FaceColor', face_clr, ...
        'EdgeColor', 'k', 'Curvature', curv);
    % draw another circle for pbc
    cond = [x(i)<R 1-x(i)<R y(i)<R 1-y(i)<R];
        % sides
    if cond(1)
        rectangle('Position', [Lx+x(i)-R y(i)-R 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    elseif cond(2)
        rectangle('Position', [-Lx+(x(i)-R) y(i)-R 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    end
    if cond(3)
        rectangle('Position', [x(i)-R Ly+y(i)-R 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    elseif cond(4)
        rectangle('Position', [x(i)-R -Ly+(y(i)-R) 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    end
    % corners
    if cond(1) && cond(3)
        rectangle('Position', [Lx+x(i)-R Ly+(y(i)-R) 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    elseif cond(2) && cond(4)
        rectangle('Position', [-Lx+(x(i)-R) -Ly+(y(i)-R) 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    elseif cond(1) && cond(4)
        rectangle('Position', [Lx+(x(i)-R) -Ly+(y(i)-R) 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
    elseif cond(2) && cond(3)
        rectangle('Position', [-Lx+(x(i)-R) Ly+(y(i)-R) 2*R 2*R], 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);    
    end
end
hold off
drawnow;
%}
%% Check whether final configuration is different enough from initial config
%{
% Delta r = distance between initial and final position of particle
% Not accurate, does not include periodic boundaries
dr = sqrt( (x-pos(:,1) ).^2  + (y-pos(:,2) ).^2 );
h=figure();
histogram(dr, 'Normalization', 'pdf');
title( strcat('$$\langle \Delta r \rangle$$ = ', sprintf('%.2f', mean(mean(dr)))), 'Interpreter', 'latex');
xlabel('$\Delta r$');
ylabel('Probability');
set(gca, 'FontSize', 16);
set(h, 'Position', [500 200 840/1.5 720/1.5]);
%}
%% Test whether distribution is random enough
%{
dist2 = reshape(dist(dist>0), [N-1 N]);
nnd = min(dist2);
nnd_av = mean(nnd);

h=figure();
hold on
histogram(nnd, 'Normalization', 'pdf');
xlabel('NND');
ylabel('Probability');
title(sprintf('%s NND %s = %.2f', '\langle', '\rangle', nnd_av));
set(gca, 'FontSize', 24);
set(h, 'Position', [500 200 840 720]);

% Clark & Evans calculation
rho = @(r) 2*pi*N/L^2*r.*exp(-pi*N/L^2*r.^2);
rvals = linspace(0, max(nnd)*1.2, 1000);
plot(rvals, rho(rvals), 'LineWidth', 2);
%xlim([0 1]);

% set image properties
%h = gcf;
%set(h,'Units','px');
%set(h, 'Position', [500 300 600 600]);
%}