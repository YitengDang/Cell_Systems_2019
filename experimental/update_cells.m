function cells_new = update_cells(cells, kappa, r, T)
   cells_new = zeros(size(cells));
   l = size(cells, 1);
   [dx, dy] = meshgrid(-r:r, -r:r); % Moore neighborhood
   
    parfor i=1:l^2
        [cell_x, cell_y] = ind2sub([l l], i);
        x_nei = mod(cell_x+dx-1,l)+1; % x-indices of neighbours
        y_nei = mod(cell_y+dy-1,l)+1; % y-indices of neighbours
        cells_nei = sub2ind(size(cells),x_nei,y_nei); % linear indices of cells

        cell_state = cells(cell_x, cell_y);
        new_state = mod(cell_state+1, kappa);
        if sum(sum(cells(cells_nei)==new_state)) >= T
            cells_new(i) = new_state;
        else
            cells_new(i) = cells(cell_x, cell_y);
        end
    end
    
end