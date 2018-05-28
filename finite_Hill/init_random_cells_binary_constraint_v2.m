function cells = init_random_cells_binary_constraint_v2(N, n, fp)
    % v2: keeps <Xi> constant (like in Monte Carlo algorithm) by adjusting
    % two cells at a time.
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
    r = randi(N, mcsteps, 2);
    %accepted = 0;
    for i=1:mcsteps
        %delta = (2*rand()-1)/2; % -0.5 <= delta <= 0.5
        alpha = min(fp(3)-fp(2), fp(2)-fp(1))/3; % arbitrary choice, tuned to give satisfactory fraction of accepted statess
        delta = alpha*randn(); % delta ~ N(0, minimum distance to unstable fixed point) 
        cells_new = [cells(r(i, 1))+delta; cells(r(i, 2))-delta];
        states_old = cells(r(i, :)) > fp(2); 
        states_new = cells_new > fp(2);
        if all(states_old == states_new) && all(abs(2*cells_new-1) <= 1) % if both cells still in range [0, 1]
            %accepted = accepted+1;
            cells(r(i, :)) = cells_new;
        end
    end
    %fprintf('accepted = %d out of %d trials \n', accepted, mcsteps);
end