function [p_out, Peq] = EOM_update_markov(N, M_int, Con, Coff, K, fN, gN, p_in, I_in)
% Update the system variables using a Markovian approach (calculate p^{ij}
% using transition matrix)
    if nargin<9
        I_in = [0 0];
    end
    %p_in = n_in/N;
    
    % get transition matrix
    W = transition_prob_states_two_signals_analytical_calc(M_int, Con, Coff, K, fN, gN, p_in, I_in);
    % (0,0), (0,1), (1,0), (1,1)
    %disp(W);
    
    % Calculate equilibrium probability
    p_out = W'*p_in(:);
    p_out = reshape(p_out, 2, 2);
    
    % Calculate Peq
    n_out = round(N*p_out); % approximate, might have rounding error
    W = transition_prob_states_two_signals_analytical_calc(M_int, Con, Coff, K, fN, gN, p_out);
    iniON = reshape(n_out, 4, 1); % (0,0), (1,0), (0,1), (1,1)
    Peq = prod(diag(W).^iniON);
end