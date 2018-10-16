function update_cell_figure(hin, pos, a0, cells, cell_type, t)

    % Plot the cell state in a diagram according to their type. If
    % cell_type(i) = 0, the cell is plotted as circle, if 1 it is plotted
    % as square. hin is the figure handle.
    
    % clear figure and hold on to keep the graphics
    clf(hin,'reset');
    title(sprintf('N_{ON} = %d, Time: %d', sum(cells), t), ...
        'FontSize', 24);
    set(gca,'YTick',[],'XTick',[]);
    box on
    hold on
    
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
    end
    
    % set hold off and draw figure
    hold off
    drawnow;

    % set image properties
    h = gcf;
    %set(h,'Units','px');
    set(h, 'Position', [500 300 560 420]);
    
    
        
                
    
    

