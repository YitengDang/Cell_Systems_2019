function [pos, dist, n_rej] = initial_cells_random_periodic(n, Lx, Ly, R)
%% Places cells randomly in a continuous space
% Parameters
N = n^2;
%N = 10;
%L = 1;
%R = 0.05;

n_rej = 0; % number of rejections
rejected = 1;
while rejected
    % Place cells
    %pos = Lx*rand(N, 2);
    x = Lx*rand(N, 1); % x-coordinates
    y = Ly*rand(N, 1); % y-coordinates

    % Calculate distance
    dist = calc_dist_periodic(x, y, Lx, Ly);

    % Check if config is allowed
    if all(dist(dist>0) >= 2*R) 
        rejected = 0;
    else
        n_rej = n_rej + 1;
        if mod(n_rej, 100)==0
            fprintf('Rejections: %d \n', n_rej);
        end
    end
end
pos = [x(:) y(:)];
% adjust distances to scale so that nearest neighbour distances are
% 1 (in case of a perfect lattice); without normalization this
% would be Lx/n.
dist = dist/(Lx/n); 
end

%% Draw configuration
%{
cells = zeros(1,N);
cells(randperm(N, round(N/2))) = 1;
hin = figure();
update_figure_periodic(hin, pos, N, L, R, cells);
%}
% set image properties
%h = gcf;
%set(h,'Units','px');
%set(h, 'Position', [500 300 600 600]);
%}