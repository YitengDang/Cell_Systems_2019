function [pos, dist, n_rej] = initial_cells_random_periodic_alt(n, Lx, Ly, R)
%% Places cells randomly in a continuous space
% Alternative approach: place cells one at a time

% Parameters
N = n^2;
%N = 10;
%L = 1;
%R = 0.05;

%initial placement
x_all = 0;
y_all = 0; 

n_rej = 0; % number of rejections
placed = 1;
while placed<N
    % Place cells
    %pos = Lx*rand(N, 2);
    rejected = 1;
    while rejected
        x = Lx*rand(); % x-coordinates
        y = Ly*rand(); % y-coordinates
        
        % Calculate distance
        dist = calc_dist_periodic([x_all; x], [y_all; y], Lx, Ly);

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
    x_all = [x_all; x];
    y_all = [y_all; y];
    placed = placed + 1;
    %fprintf('Placed cells: %d \n', placed);
end
pos = [x_all(:) y_all(:)];
% adjust distances to scale so that nearest neighbour distances are
% 1 (in case of a perfect lattice); without normalization this
% would be Lx/n.
dist = dist/(Lx/n); 

fprintf('Placed cells: %d \n', placed);
fprintf('Total rejections: %d \n', n_rej);
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