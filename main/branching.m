function [pn, t_av] = branching(pmatrix, pEn, branch_prob, n, pn, t_av, t)
% Perform the bhanching algorithm iteratively
    eps = 1e-8;
    % Chance of equilibrium at this n
    pn(n+1) = pn(n+1) + branch_prob*pEn(n+1);
    t_av = t_av + branch_prob*pEn(n+1)*t;
    
    % Now check each branch and branch if needed
    N = size(pmatrix,1)-1;
    % if the combination of all branches (1-pEn) is below the limiting value
    % do not start branching
    if branch_prob*(1-pEn(n+1)) > eps
        for nf = 0:N
            % multiplicative probability of the branch
            pbranch = branch_prob*pmatrix(nf+1, n+1);
            if pbranch > eps
                % if branch probability is above threshold
                [pn, t_av] = branching(pmatrix, pEn, pbranch, nf, pn, t_av, t+1);
            end
        end
    end
