function update_figure(hin, pos, N, L, R, cells)

% Convert to x, y variables
x = pos(:, 1);
y = pos(:, 2);
eta = N*pi*R^2/L^2;

% Draw configuration
%hin = figure();
clf(hin,'reset');
title(sprintf('N = %d, R = %.2f, %s = %.2f', N, R, '\eta', eta), ...
    'FontSize', 24);
%set(gca,'YTick',[],'XTick',[]);
set(gca,'DataAspectRatio', [1 1 1]);
axis([-R L+R -R L+R]);
box on
hold on
for i=1:N
    position = [x(i)-R y(i)-R 2*R 2*R];
    clr = 0.8*(1 - cells(i));
    face_clr = [clr clr clr];
    curv = [1 1];
    rectangle('Position', position, 'FaceColor', face_clr, ...
        'EdgeColor', face_clr, 'Curvature', curv);
end
hold off
drawnow;

% set image properties
%h = gcf;
%set(h,'Units','px');
%set(h, 'Position', [500 300 600 600]);

end