function update_cell_figure_count_neighbors(hin, pos, dist, a0, cells, cell_type, t)

    % Plot the cell state in a diagram according to their type. If
    % cell_type(i) = 0, the cell is plotted as circle, if 1 it is plotted
    % as square. The count of the number of ON cells is also shown. hin is
    % the figure handle.
    
    % clear figure and hold on to keep the graphics
    clf(hin,'reset');
    title(sprintf('ON: black, N_{ON} = %d, Time: %d', sum(cells), t))
    hold on
    
    % First neighbor distance
    eps = 1e-5;
    dist_vec = get_unique_distances(dist, eps);
    dist1 = dist_vec(2);
    idx = 1*(dist < dist1+eps & dist > dist1-eps);
    nn = idx*cells;
    
    % set the diameter of the graphics (unitary)
    diameter = 0.5;
    
    for i = 1:size(pos,1)
        switch cell_type(i)
            case 0
                curv = [1 1];
            otherwise
                curv = [0 0];
        end
        
        if cells(i) == 0
            r = 0.914;
            g = 0.914;
            b = 0.914;
            face_clr = [r g b];
        else
            face_clr = 'k';
        end
        
        r = diameter/2;
        position = a0.*[pos(i,1)-r pos(i,2)-r diameter diameter];
        
        % draw the rectangle in the figure with handle hin
        rectangle('Position', position, 'FaceColor', face_clr, ...
            'EdgeColor', 'k', 'Curvature', curv);
        
        text(a0*pos(i,1)-r, a0*pos(i,2), int2str(nn(i)),'Color', 'r')
    end
    
    % set hold off and draw figure
    hold off
    drawnow;
        
                
    
    

