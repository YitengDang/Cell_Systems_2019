function [p1_out, p2_out, theta11_out, theta12_out, theta22_out, t_final_out]...
    = input_output_twocelltypes_montecarlo_run(p1_in, p2_in, N1, N2, f_mat, g_mat, Coff, Con, K,...
    tmax, nruns)

dim1 = numel(p1_in);
dim2 = numel(p2_in);

% variables to save
p1_out = zeros(dim1, dim2);
p2_out = zeros(dim1, dim2);
theta11_out = zeros(dim1, dim2);
theta12_out = zeros(dim1, dim2);
theta22_out = zeros(dim1, dim2);
t_final_out = zeros(dim1, dim2);
%I_out = zeros(dim1, dim2);
%h_out = zeros(dim1, dim2);
%f_mat_out = cell(nruns, 1);
%g_mat_out = cell(nruns, 1);

%---running time estimate
%dt = 1.4; % estimated time per li5
%t_est = nruns*numel(K_loop)^2*numel(Con_loop)^2*dt;
%fprintf('Estimated running time: %.f s = %.2f min = %.2f h \n', t_est, t_est/60, t_est/3600);

fprintf('K1 = %d, K2 = %d, Con1 = %d, Con2 = %d \n', K(1), K(2), Con(1), Con(2));
for li1=1:numel(p1_in)
    for li2=1:numel(p2_in)
        p_in = [p1_in(li1); p2_in(li2)];
        
        [px, py] = meshgrid(p_in, p_in);
        theta_in = f_mat.*(2*px-1).*(2*py-1);
        
        fprintf('p1 = %.2f, p2 = %.2f \n', p_in(1), p_in(2));
        this_p1 = zeros(nruns,1);
        this_p2 = zeros(nruns,1);
        this_theta11 = zeros(nruns,1);
        this_theta12 = zeros(nruns,1);
        this_theta22 = zeros(nruns,1);
        this_tf = zeros(nruns,1);
        %this_h = zeros(nruns,1);
        %this_I = zeros(nruns,1);
        parfor li=1:nruns
            %disp(li);
            try
                t = 0;
                p = p_in;
                theta = theta_in
                [p_out, theta_out, ~, cont] = update_monte_carlo_with_theta(N1, N2, p, theta, Con, Coff, K, f_mat, g_mat);
                while cont && t<tmax
                    t = t+1;
                    p = p_out;
                    theta = theta_out;
                    [p_out, theta_out, ~, cont] = update_monte_carlo_with_theta(N1, N2, p, theta, Con, Coff, K, f_mat, g_mat);       
                end
                this_p1(li) = p_out(1);
                this_p2(li) = p_out(2);
                this_theta11(li) = theta_out(1,1);
                this_theta12(li) = theta_out(1,2);
                this_theta22(li) = theta_out(2,2);
                this_tf(li) = t;
                %this_h(li) = h;
                %this_I(li) = I;
                %f_mat_out{li} = f_mat;
                %g_mat_out{li} = g_mat;
            catch 
                disp('error parfor loop');
            end
        end
        p1_out(li1, li2) = mean(this_p1);
        p2_out(li1, li2) = mean(this_p2);
        theta11_out(li1, li2) = mean(this_theta11);
        theta12_out(li1, li2) = mean(this_theta12);
        theta22_out(li1, li2) = mean(this_theta22);
        t_final_out(li1, li2) = mean(this_tf);
        %disp(this_tf);
        %h_out(li1, li2) = mean(this_h);
        %I_out(li1, li2) = mean(this_I);
    end
end