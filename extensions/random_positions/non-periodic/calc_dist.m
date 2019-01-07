function dist = calc_dist(pos)
    x = pos(:, 1); % x-coordinates
    y = pos(:, 2); % y-coordinates

    % Calculate distance
    [x1,x2] = meshgrid(x, x);
    [y1,y2] = meshgrid(y, y);
    dist = sqrt((x1 - x2).^2 + (y1 - y2).^2);
end