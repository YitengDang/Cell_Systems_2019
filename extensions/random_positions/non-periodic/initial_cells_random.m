function [pos, dist, rejections] = initial_cells_random(N, L, R)
% Places cells randomly in a continuous space
% Parameters
%N = 10;
%L = 1;
%a0 = L^4/N;
%R = 0.05; %0.2*a0;

rejections = 0;
rejected = 1;
while rejected
    % Place cells
    pos = L*rand(N, 2);
    x = pos(:, 1); % x-coordinates
    y = pos(:, 2); % y-coordinates

    % Calculate distance
    [x1,x2] = meshgrid(x, x);
    [y1,y2] = meshgrid(y, y);
    dist = sqrt((x1 - x2).^2 + (y1 - y2).^2);

    % Check if config is allowed
    if all(dist(dist>0) >= 2*R) 
        rejected = 0;
    else
        rejections = rejections + 1;
        if mod(rejections, 100)==0
            fprintf('Rejections: %d \n', rejections);
        end
    end
end

end

%% Draw configuration
%{
hin = figure();
clf(hin,'reset');
title(sprintf('N = %d, R = %.2f', N, R), ...
    'FontSize', 24);
set(gca,'YTick',[],'XTick',[]);
set(gca,'DataAspectRatio', [1 1 1]);
axis([0 L 0 L]);
box on
hold on
for i=1:N
    position = [x(i) y(i) 2*R 2*R];
    face_clr = 'k';
    curv = [1 1];
    rectangle('Position', position, 'FaceColor', face_clr, ...
        'EdgeColor', 'k', 'Curvature', curv);
end
hold off
drawnow;

% set image properties
%h = gcf;
%set(h,'Units','px');
%set(h, 'Position', [500 300 600 600]);
%}