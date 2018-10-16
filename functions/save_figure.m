function save_figure(h_fig, width, height, path_out, ext, qsave, colored_background)
    if nargin<7
        colored_background = 0;
    elseif nargin<6
        qsave = 1;
        colored_background = 0;
    elseif nargin<5
        ext = '.pdf';
        qsave = 1;
        colored_background = 0;
    end
    
    if ~qsave
        return
    end
    
    if width==0 && height==0 % use input figure size        
        h_fig.Units = 'Inches';
        width = h_fig.Position(3);
        height = h_fig.Position(4);
    end
    
    % Adjust the height and width of a figure and save it in PDF
    set(h_fig, 'Units', 'Inches');
    set(h_fig, 'Position', [0 0 width height])
    set(h_fig,'PaperPositionMode','Auto', 'PaperUnits', 'Inches', 'PaperSize', [width, height])
        
    if colored_background
        % save with current figure background color (instead of white)
        h_fig.InvertHardcopy = 'off'; 
    end
    
    if strcmp(ext, '.pdf')
        print(h_fig, path_out,'-dpdf','-r0');
    elseif strcmp(ext, '.eps')
        print(h_fig, path_out,'-depsc','-r0');
    elseif strcmp(ext, '.svg')
        print(h_fig, path_out,'-dsvg','-r0');
    end
    
    fprintf('Figure saved at %s \n', strcat(path_out, ext));
end