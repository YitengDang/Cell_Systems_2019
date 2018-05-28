function update_cell_figure_save(hin, pos, a0, cells, cell_type, t, B, fname, frame)
    % Plot the cell state in a diagram according to their type. If
    % cell_type(i) = 0, the cell is plotted as circle, if 1 it is plotted
    % as square. hin is the figure handle.
    
    % clear figure and hold on to keep the graphics
    
    % Saves figure frame; not used anymore because MATLAB has built-in
    % movie function
    clf(hin,'reset');
    title(sprintf('N_{ON} = %d, Time: %d, B = %.2f', sum(cells), t, B), ...
        'FontSize', 24)
    set(gca,'YTick',[],'XTick',[])
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
    
    % save frame
    h = gcf;
    %set(h,'Units','px');
    set(h, 'Position', [500 500 560 420]);
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);    
    i = 1;
    out_file = fullfile(pwd, 'figures', 'time_evolution', ...
        strcat(fname,int2str(frame),'-v',int2str(i)));
    %disp(out_file);
    while exist(strcat(out_file,'.fig'), 'file') == 2
        i=i+1;
        %disp(i);
        out_file = fullfile(pwd, 'figures', 'time_evolution', ...
            strcat(fname,int2str(frame),'-v',int2str(i)));
    end
    savefig(h, out_file);
    print(h, out_file,'-djpeg','-r0');
    
        
                
    
    

