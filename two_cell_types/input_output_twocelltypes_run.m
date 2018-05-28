function [p1_out, p2_out, theta11_out, theta12_out, theta22_out, I_out, t_final_out, h_out]...
    = input_output_twocelltypes_run(p1_in, p2_in, N1, N2, dist, a0, Con, K, Rcell,...
    Mcomm, tmax, nruns)

dim1 = numel(p1_in);
dim2 = numel(p2_in);

% variables to save
p1_out = zeros(dim1, dim2);
p2_out = zeros(dim1, dim2);
theta11_out = zeros(dim1, dim2);
theta12_out = zeros(dim1, dim2);
theta22_out = zeros(dim1, dim2);
I_out = zeros(dim1, dim2);
t_final_out = zeros(dim1, dim2);
h_out = zeros(dim1, dim2);
%f_mat_out = cell(nruns, 1);
%g_mat_out = cell(nruns, 1);

%---running time estimate
%dt = 1.4; % estimated time per li5
%t_est = nruns*numel(K_loop)^2*numel(Con_loop)^2*dt;
%fprintf('Estimated running time: %.f s = %.2f min = %.2f h \n', t_est, t_est/60, t_est/3600);

fprintf('K1 = %d, K2 = %d, Con1 = %d, Con2 = %d \n', K(1), K(2), Con(1), Con(2));
for li1=1:numel(p1_in)
    for li2=1:numel(p2_in)
        p = [p1_in(li1) p2_in(li2)];
        fprintf('p1 = %.2f, p2 = %.2f \n', p(1), p(2));
        this_p1 = zeros(nruns,1);
        this_p2 = zeros(nruns,1);
        this_theta11 = zeros(nruns,1);
        this_theta12 = zeros(nruns,1);
        this_theta22 = zeros(nruns,1);
        this_h = zeros(nruns,1);
        this_tf = zeros(nruns,1);
        this_I = zeros(nruns,1);
        parfor li=1:nruns
            %disp(li);
            try
            [p1, p2, theta11, theta12, theta22, h, t_final, I] = time_evolution_twocelltypes_in_out(...
            	N1, N2, p, dist, a0, Con, K, Rcell, Mcomm, tmax);            
                this_p1(li) = p1;
                this_p2(li) = p2;
                this_theta11(li) = theta11;
                this_theta12(li) = theta12;
                this_theta22(li) = theta22;
                this_h(li) = h;
                this_tf(li) = t_final;
                this_I(li) = I;
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
        h_out(li1, li2) = mean(this_h);
        %disp(this_tf);
        t_final_out(li1, li2) = mean(this_tf);

        I_out(li1, li2) = mean(this_I);
    end
end