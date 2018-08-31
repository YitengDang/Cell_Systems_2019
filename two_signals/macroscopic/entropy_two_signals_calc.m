function [Omega_E, S, frac_E] = entropy_two_signals_calc(N, M_int, Con, Coff, K, fN, gN)
    Omega_E_all = [];
    for n00=0:N
        %disp(n00);
        for n01=0:N-n00
            for n10=0:N-n00-n01
                n = [n00 n01; n10 N-n00-n01-n10];
                %disp(n);
                pij = n/N;

                W = transition_prob_states_two_signals_analytical_calc(M_int, Con, Coff, K, fN, gN, pij);
                iniON = reshape(round(pij.*N)',4,1); % iniON: (0,0), (0,1), (1,0), (1,1)
                Peq = prod(diag(W).^iniON); % diag(W): W(0,0), W(0,1), W(1,0), W(1,1)
                %disp(Peq)

                mult_coeff = prod([1:N]./[1:n(1) 1:n(2) 1:n(3) 1:n(4)]); 
                Omega_En = mult_coeff*Peq; 
                
                Omega_E_all(end+1) = Omega_En;
            end
        end
    end

    if all(Omega_E_all<1)
        Omega_E_log = -Inf;
    end
    Omega_E = sum(Omega_E_all);
    disp(Omega_E);
    %% Get fraction of equilibrium states
    S = log(Omega_E);
    Smax = N*log(2);
    frac_E = exp(S - Smax);
    disp(frac_E);
    
    %% Check all Omega_E
    % figure;
    % histogram(Omega_E_all);
end
