function cells = init_random_cells_binary_constraint(N, n, fp)
    % n: number of ON cells
    % Monte Carlo algorithm to generate random configuration with fixed N
    % and fraction of ON cells (determined by basin of attraction, determined by fixed points)
    % Works for Hill coefficient n=2 
    if numel(fp)~=3
        fprintf('System not binary / input fixed points wrong');
        return
    end
    %%
    % If p is not 0 or 1
    mcsteps = 10^4;
    %n = round(p*N);
    cells = fp(1)*ones(N, 1);
    cells( randperm(N, n) ) = fp(3);
    r = randi(N, mcsteps, 1);
    %accepted = 0;
    alpha = min(fp(3)-fp(2), fp(2)-fp(1))/3; % arbitrary choice, tuned to give satisfactory fraction of accepted statess
    for i=1:mcsteps
        %delta = (2*rand()-1)/2; % -0.5 <= delta <= 0.5
        delta = alpha*randn(); % delta ~ N(0, minimum distance to unstable fixed point) 
        cell_new = cells(r(i)) + delta;
        state_old = cells(r(i)) > fp(2); 
        state_new = cell_new > fp(2);
        if all(state_old == state_new) && all(abs(2*cell_new-1) <= 1) % if cell maintains state and is in valid range
            cells(r(i)) = cell_new;
        end
    end
    %fprintf('accepted = %d out of %d trials \n', accepted, mcsteps);
end