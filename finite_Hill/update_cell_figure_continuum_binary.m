function update_cell_figure_continuum_binary(hin, pos, a0, cells, cell_type, t, threshold)
    % Plot the cell state in a diagram according to their type. If
    % cell_type(i) = 0, the cell is plotted as circle, if 1 it is plotted
    % as square. hin is the figure handle.
    
    % _binary: plots cells as binary (black = ON, white = OFF) depending on
    % whether the states are above a certain threshold (unstable fixed
    % point)
    
    % clear figure and hold on to keep the graphics
    clf(hin,'reset');
    title(strcat('$$\langle X_i \rangle$$',...
        sprintf('= %.2f, Time: %d', mean(cells), t)), ...
        'FontSize', 24, 'Interpreter', 'latex');
    %title(sprintf('Noff = %d, Non =  %d,  Time: %d', sum(~cells), sum(cells), t), ...
    %    'FontSize', 24, 'Interpreter', 'latex');
    set(gca,'YTick',[],'XTick',[]);
    box on
    hold on

    
    % set the diameter of the graphics (unitary)
    diameter = 0.5;
    
    face_clr_cells = (cells < threshold); %cells below threshold become white
    
    for i = 1:size(pos,1)
        switch cell_type(i)
            case 0
                curv = [1 1];
            otherwise
                curv = [0 0];
        end
        
        % greyscale, cells(i) = Xi indicates greyness
        face_clr = face_clr_cells(i)*ones(1,3);
        
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
    
    
        
                
    
    

