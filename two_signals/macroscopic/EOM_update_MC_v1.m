function [n_out, Peq] = EOM_update_MC(N, M_int, Con, Coff, K, fN, gN, n_in, I_in)
    % n_in format: index [1,1] corresponds to state (0,0), [1,2] = (0,1), [2,1]=(1,0), [2,2]=(1,1)
    % in linear indices, n_in => (0,0), (1,0), (0,1), (1,1)
    if nargin<9
        I_in = [0 0];
    end
    p_in = n_in/N;
    
    % get transition matrix
    W = transition_prob_two_signals_pij(M_int, Con, Coff, K, fN, gN, p_in, I_in);
    % (0,0), (0,1), (1,0), (1,1)
    %disp(W);
    
    % obtain output states
    n_out = zeros(2);
    for k=1:4
        %disp(n_in(k));
        p = [0 cumsum(W(k, :), 2)];
        r = rand(n_in(k), 1);
        for i=1:n_in(k) % can we vectorize this?
            out_idx = find(r(i) < p, 1)-1;
            n_out(out_idx) = n_out(out_idx) + 1;
        end
    end
    
    % Calculate equilibrium probability
    p_out = n_out/N;
    W = transition_prob_two_signals_pij(M_int, Con, Coff, K, fN, gN, p_out);
    iniON = reshape(n_out', 4, 1); % n_out -> [(0,0), (0,1); (1,0), (1,1)] 
    Peq = prod(diag(W).^iniON);
end