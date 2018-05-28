function [p1_out, p2_out, theta11_out, theta12_out, theta22_out, t_out] = ...
    time_evolution_twocelltypes_montecarlo_in_out(N1, N2, p_in, theta_in, f_mat, g_mat, Con, Coff, K, tmax)
    % return_output: returns only the output variables and doesn't save the
    % entire trajectory

    N = N1+N2;
    
    %% Trajectory    
    t = 0;
    [p_out, theta_out, pe, cont] = update_monte_carlo_with_theta(N1, N2, p_in, theta_in, Con, Coff, K, f_mat, g_mat);
    while cont && t<tmax
        t = t+1;
        p_in = p_out;
        theta_in = theta_out;
        [p_out, theta_out, pe, cont] = update_monte_carlo_with_theta(N1, N2, p_in, theta_in, Con, Coff, K, f_mat, g_mat);       
    end
        
    %% Values to return
    p1_out = p_out(1);
    p2_out = p_out(2);
    theta11_out = theta_out(1,1);
    theta12_out = theta_out(1,2);
    theta22_out = theta_out(2,2);
    t_out = t;
    pe_out = pe;

    