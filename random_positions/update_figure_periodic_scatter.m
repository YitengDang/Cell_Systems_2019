function update_figure_periodic_scatter(pos, N, Lx, Ly, R, cells, t)
% Short version: plot using scatter
% Not guaranteed to work yet

% Convert to x, y variables
%x = pos(:, 1);
%y = pos(:, 2);
%eta = N*pi*R^2/L^2;

% Draw configuration
%clf(hin,'reset');
title(sprintf('N = %d, R = %.2f, time = %d', N, R, t), ...
    'FontSize', 24);
%set(gca,'YTick',[],'XTick',[]);
set(gca,'DataAspectRatio', [1 1 1]);
axis([0 Lx 0 Ly]);
box on
hold on

map = [0,0,0
    1,1,1];
colormap(map);
scatter(pos(:,1), pos(:,2), 50, 'MarkerEdgeColor', 'k');
scatter(pos(:,1), pos(:,2), 50, 1-cells, 'filled');
%{
for i=1:N
    position = [x(i)-R y(i)-R 2*R 2*R];
    clr = (1 - cells(i));
    face_clr = [clr clr clr];
    %face_clr = 'k';
    curv = [1 1];
    rectangle('Position', position, 'FaceColor', face_clr, ...
        'EdgeColor', 'k', 'Curvature', curv);
    % draw another circle for pbc
    cond = [x(i)<R Lx-x(i)<R y(i)<R Ly-y(i)<R];
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
    %
end
drawnow;
%}
hold off

set(gca, 'XTick', [], 'YTick', []);
set(gcf, 'Position', [500 300 600 600]);
% set image properties
%h = gcf;
%set(h,'Units','px');
%set(h, 'Position', [500 300 600 600]);

end