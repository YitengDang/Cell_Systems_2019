function bin = dos_Wang_Landau_binning(x, edges)
    % Determines the number of the bin the value x falls, given edges of
    % the bins. Takes a single numerical value x.
    nbins = length(edges) - 1;
    if x < edges(1) 
        bin = 1;
    elseif x > edges(end)
        bin = nbins;
    else
        for i=1:nbins
            if x <= edges(i+1)
                bin = i;
                return
            end
        end
    end
end