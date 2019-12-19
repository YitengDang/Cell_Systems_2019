function save_movie_gif(cells_hist, rcell, pos, pos_hist, disp_mol, fname_out,...
    frame_rate, t_ini, t_out)
    %frames = struct('cdata',[],'colormap',[]);
    
    if nargin<9
        t_ini = 0;
        t_out = length(cells_hist)-1;
    end
    % Options
    %frame_rate = 5; % frames/second
    % format = 'MPEG-4'; %'Motion JPEG AVI'; %movie format 
    % 'Motion JPEG AVI' <- default, works best
    % 'Uncompressed AVI' <- high quality(?), large file
    % 'MPEG-4' <- .mp4
    % 'Archival' <- unknown ext
    % 'Motion JPEG 2000' <- unknown ext
    
    
    % replay trajectory externally
    %tt = 0;
    h = figure;
    clf(h, 'reset');
    if isempty(pos_hist)
        % same positions every time step
        [h_cells, h_borders] = reset_cell_figure(h, pos, rcell);
        for tt=t_ini+1:t_out+1
            cells = cells_hist{tt};
            update_cell_figure_external(h_cells, h_borders, cells, tt-1, disp_mol, pos);
            
            frame = getframe(gcf);
            img =  frame2im(frame);
            [img, cmap] = rgb2ind(img,256);
            if tt == t_ini+1
                imwrite(img,cmap,strcat(fname_out, '.gif'),'gif',...
                    'LoopCount',Inf,'DelayTime', 1/frame_rate, 'Compression', 'jpeg');
            else
                imwrite(img,cmap,strcat(fname_out, '.gif'),'gif',...
                    'WriteMode','append','DelayTime', 1/frame_rate, 'Compression', 'jpeg');
            end
            
        end
    else
        % new positions every time step
        [h_cells, h_borders] = reset_cell_figure(h, pos, rcell);
        for tt=t_ini+1:t_out+1
            cells = cells_hist{tt};
            pos = pos_hist{tt};
            update_cell_figure_external(h_cells, h_borders, cells, tt-1, disp_mol, pos);
            
            frame = getframe(gcf);
            img =  frame2im(frame);
            [img, cmap] = rgb2ind(img,256);
            if tt == t_ini+1
                imwrite(img,cmap,strcat(fname_out, '.gif'),'gif',...
                    'LoopCount',Inf,'DelayTime', 1/frame_rate, 'Compression', 'jpeg');
            else
                imwrite(img,cmap,strcat(fname_out, '.gif'),'gif',...
                    'WriteMode','append','DelayTime', 1/frame_rate, 'Compression', 'jpeg');
            end
        end
    end
    
    % add final state to show equilibrium (not applicable if not in
    % equilibrium)
    %{
    cells = cells_hist{t};
    t = t+1;
    update_cell_figure_external(h, pos, cells, cell_type, t-1, disp_mol, rcell)     
    frames(t) = getframe(h);
    frame = getframe(h);
    %}
end