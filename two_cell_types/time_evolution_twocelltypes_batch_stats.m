% Batch simulations of a large collection of trajectories
% Save only the main statistics of the simulations, not detailed
% trajectories (too many data points)
close all
clear all
warning off
%--------------------------------------------------------------------------
% Set parameters of the system
gridsize = 15;
N = gridsize^2;
frac_type1 = 0.5;
N1 = round(frac_type1*N);
N2 = N - N1;
a0 = 1.5;

% Parameters 
Con = [20 20]; %column j = type j 
Coff = [1 1];
K = [2 5];
p = [0.5 0.5];
Rcell = [0.2 0.2]*a0;
%K_loop = 10; %[2 5 10]; %2:15;
%Con_loop = 10; %5:5:40; %2:1:40;

% Simulation parameters
nruns = 1000;
tmax = 1000;
save_trajectory = 0; %whether to save individual trajectories 
 
% Communication matrix
Mcomm = [0 1; 1 0]; % type i reacts to type j iff M_ij=1
Mcomm = logical(Mcomm);
Mcomm_str = sprintf('Mcomm_%d_%d_%d_%d', Mcomm(1,1), Mcomm(1,2), Mcomm(2,1), Mcomm(2,2));

% Check input
if ~all(size(Mcomm)==[2 2])
    disp('Wrong communication matrix size!');
    return
end
for idx=1:numel(Mcomm)
    if ~(Mcomm(idx)==0 || Mcomm(idx)==1)
        disp('Wrong communication matrix entries!');
        return
    end
end

% Distance and position
[dist, ~] = init_dist_hex(gridsize, gridsize);
%--------------------------------------------------------------------------
%% Run trajectories
% variables to save
p1_out = cell(nruns, 1);
p2_out = cell(nruns, 1);
theta_out = cell(nruns, 3);
I_out = cell(nruns, 1);
t_final_out = zeros(nruns, 1);
h_out = cell(nruns, 1);
f_mat_out = cell(nruns, 1);
g_mat_out = cell(nruns, 1);
savepath = fullfile(pwd, 'twocelltypes', 'data', 'time_evolution_stats_twocelltypes', Mcomm_str);

%---running time estimate
%dt = 1.4; % estimated time per li5
%t_est = nruns*numel(K_loop)^2*numel(Con_loop)^2*dt;
%fprintf('Estimated running time: %.f s = %.2f min = %.2f h \n', t_est, t_est/60, t_est/3600);

fprintf('K1 = %d, K2 = %d, Con1 = %d, Con2 = %d \n', K(1), K(2), Con(1), Con(2));
for li=1:nruns
    disp(li);
    [p1, p2, theta11, theta12, theta22, h, t_final, I, f_mat, g_mat] = time_evolution_twocelltypes_full_v2(...
        N, dist, a0, Con, K, Rcell, frac_type1, Mcomm, p, tmax, savepath, save_trajectory);
    p1_out{li} = p1;
    p2_out{li} = p2;
    theta_out{li, 1} = theta11;
    theta_out{li, 2} = theta12;
    theta_out{li, 3} = theta22;
    h_out{li} = h;
    t_final_out(li) = t_final;
    I_out{li} = I;
    f_mat_out{li} = f_mat;
    g_mat_out{li} = g_mat;
    %[p1_out{li}, p2_out{li}, theta_out{li,:}, h_out{li}, t_final_out{li}, I_out{li}] = ...
    %    time_evolution_twocelltypes_full_v2(N, dist, a0, Con, K, Rcell, ...
    %    frac_type1, Mcomm, p, tmax, savepath, save_trajectory);
end
%}
%% Calculate Peq
Peq_out = zeros(nruns, 2);
error_p_avg = zeros(nruns, 2); % i.e. column 1: avg. error in p1
for li1=1:nruns 
    %i=1;
    % define quantities
    p1 = p1_out{li1};
    p2 = p2_out{li1};
    theta_11 = theta_out{li1, 1};
    theta_12 = theta_out{li1, 2};
    theta_22 = theta_out{li1, 3};
    f_mat = f_mat_out{li1};
    g_mat = g_mat_out{li1};

    %-------------------Calculate Peq--------------------------------------
    p = [p1(end); p2(end)];
    theta = [theta_11(end) theta_12(end) theta_22(end)];
    theta_mat = [theta(1) theta(2); theta(2) theta(3)]; 

    N_12 = [N1; N2];
    xl = [N1/N; N2/N];

    Xlm_on = diag(1./(2.*p))*(theta_mat + repmat(2*p'-1, 2, 1));
    Xlm_off = diag(1./(2.*(1-p)))*(-theta_mat + repmat(2*p'-1, 2, 1));

    S_on = repmat((Con-Coff)/2, 2, 1).*Xlm_on + repmat((Con+Coff)/2, 2, 1);
    S_off = repmat((Con-Coff)/2, 2, 1).*Xlm_off + repmat((Con+Coff)/2, 2, 1);

    muon = sum(diag(1./xl)*f_mat.*S_on, 2);
    muoff = sum(diag(1./xl)*f_mat.*S_off, 2);

    S2 = ((Con'-Coff').^2.*p.*(1-p));
    sigmaon = sqrt(diag(1./xl).*g_mat*S2);
    sigmaoff = sigmaon;

    zon = (K' - Con' - muon)./sigmaon;
    zoff = (K' - 1 - muoff)./sigmaoff;
    Poffoff = normcdf(zoff);
    Ponon = 1-normcdf(zon);

    % if one of the probabilities is not defined (because there are no cells of
    % a certain type)
    if sum(isnan(Poffoff))>0 
        pe = Ponon.^N_12;
    elseif sum(isnan(Ponon))>0
        pe = Poffoff.^N_12;
    else
        pe = (Ponon.^p.*Poffoff.^(1-p)).^N_12;
    end
    Peq_out(li1, :) = pe;
    
    % ---------------Calculate average error in p1, p2---------------------
    error_p = zeros(numel(p1), 2);
    for li2=1:numel(p1)-1
        p_in = [p1(li2); p2(li2)];
        p_out = [p1(li2+1); p2(li2+1)]; %expected p_out
        theta_in = [theta_11(li2) theta_12(li2) theta_22(li2)];
        [p_MC_out] = ...
            update_monte_carlo_with_theta(N1, N2, p_in, theta_in, Con, Coff, K, f_mat, g_mat);
        error_p(li2, :) = abs(p_MC_out - p_out);
    end
    error_p_avg(li1, :) = mean(error_p, 1); 
end

%%
% save statistics
qsave = 1;
if qsave
    fname_str = strrep(sprintf('stats_out_N1_%d_N2_%d_a0_%.1f_K1_%d_K2_%d_Con1_%d_Con2_%d_nruns_%d_%s',...
        N1, N2, a0, K(1), K(2), Con(1), Con(2), nruns, Mcomm_str), '.', 'p');
    savepath = fullfile(pwd, 'data', 'time_evolution_stats_twocelltypes');
    % 'L:\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes';
    %path ='K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
    fname = fullfile(savepath, fname_str);
    %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
    clear dist
    save(fname);
end
%% Plot statistics
% Histogram of P_eq
h1=figure(1);
edges = 0:0.1:1;
histogram(Peq_out(:,1), edges);
hold on
histogram(Peq_out(:,2), edges);
xlabel('$$P_{eq}(t_f)$$');
ylabel('frequency');
legend({'1', '2'}, 'Location', 'nw');

%% Histogram of average errors
h2 = figure(2);
%edges = 0:0.1:1;
histogram(error_p_avg(:,1));
hold on
histogram(error_p_avg(:,2));
xlabel('$$\Delta_p$$');
ylabel('frequency');
legend({'1', '2'}, 'Location', 'nw');
%%
h3 = figure(3);
%edges = 0:0.1:1;
histogram(error_p_avg(:,1).*N1);
hold on
histogram(error_p_avg(:,2).*N2);
xlabel('$$N_{ON} \Delta_p$$');
ylabel('frequency');
legend({'1', '2'}, 'Location', 'nw');

%% Histogram of t_final
h4=figure(4);
%edges = 0:0.1:1;
histogram(t_final_out);
xlabel('$$t_{final}$$');
ylabel('frequency');
%legend({'1', '2'}, 'Location', 'nw');