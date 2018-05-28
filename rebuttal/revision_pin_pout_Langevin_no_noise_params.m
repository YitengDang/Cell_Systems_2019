% Creates pin-pout map from Langevin equation
% Based on code for pin-pout map for p, I Monte Carlo simulations
clear variables
close all

% Parameters
gridsize = 15;
a0 = 0.5; %1.5 1.5 1.5 0.5 0.5 0.5
Rcell = 0.2*a0;
K = 14; %6 12 20 14 15 14
Con = 5; %15 15 15 12 8 5
%delta = 0.004; % % step size in Langevin equation, suggested 0.001
%noise = 2; % 2; % Lagevin noise amplitude, suggested 2
alpha = 0; % sensing noise

% Initial state
Ii = 0;

% Initialization
N = gridsize^2;
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
M = sinh(Rcell)*exp(Rcell-a0*dist)./(a0*dist);
M(1:N+1:N^2) = 0;
aux = M*M';
gN = aux(1,1);

trials = 1000;
cntmap = zeros(N+1);
t_av = zeros(N+1, 1);

if Ii == 0
    fname = strrep(sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d', ...
            N, Con, K, 10*a0), '.','p');
else
    fname = strrep(sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d_I%d', ...
            N, Con, K, 10*a0, 10*Ii), '.','p');
end

bad_traj = zeros(N+1, 1); % record the number of 'bad' trajectories
for n = 0:N
    fprintf('n=%d \n', n);
    for i = 1:trials
        %if mod(i, 100)==0
        %    fprintf('trial %d \n', i);
        %end
        p = n/N;
        I = Ii;
        cont = true; %continue simulation?
        t = 0;
        while cont
            [p_new, I_new, Peq, p_lim] = revision_update_Langevin_no_noise_params(p, I, N, Con, K, fN, gN, alpha);

            % system out of range?
            if p_lim == 1 || p_lim == -1
                cont = false;
                t_av(n+1) = t_av(n+1) + t/trials;
                nout = round(N*(p_lim+1)/2);
                cntmap(n+1, nout + 1) = cntmap(n+1, nout + 1) + 1/trials;
            end
            if I_new < -1 || I_new > 1
                cont = false;
                disp('Escaping trajectory');
                bad_traj(n+1) = bad_traj(n+1) + 1;
            end
            % system in equilibrium?
            if rand < Peq
                cont = false;
                t_av(n+1) = t_av(n+1) + t/trials;
                nout = round(N*p);
                cntmap(n+1, nout + 1) = cntmap(n+1, nout + 1) + 1/trials;
            else
                t = t + 1;
                p = p_new;
                I = I_new;
            end
        end
    end
end

%%
qsave = 1;
if qsave
    save(fullfile(pwd, 'rebuttal', 'Langevin_pin_pout', strcat('Langevin_no_noise_params_',fname) ));
end

h1=figure(1);
p = (0:N)./N;
im_fig = imagesc(p,p,cntmap');
set(gca,'Ydir','normal','FontSize', 20)
set(im_fig, 'AlphaData', cntmap' > 0);
c = colorbar;
c.Label.String = 'Probability';
caxis([0 1]);
xlabel('p_{ini}', 'FontSize', 24)
ylabel('p_{eq}', 'FontSize', 24)

qsave = 0;
if qsave
    out_file = fullfile(pwd, 'rebuttal', 'Langevin_pin_pout', strcat('Langevin_no_noise_params_map_',fname) );
    save_figure_pdf(h1, 10, 8, out_file);
end

h2 = figure(2);
plot(p, t_av)
xlabel('p_{in}', 'FontSize', 24)
ylabel('Time', 'FontSize', 24)
set(gca,'FontSize', 20)
if qsave
    out_file = fullfile(pwd, 'rebuttal', 'Langevin_pin_pout', strcat('Langevin_no_noise_params_tav_',fname) );
    save_figure_pdf(h2, 10, 8, out_file);
end

%% Plot bad trajectories
h3=figure(3);
plot(p, bad_traj/trials, 'LineWidth', 2);
xlabel('p_{in}', 'FontSize', 24)
ylabel('Fraction escaping trajectories', 'FontSize', 24)
set(gca,'FontSize', 20)
qsave = 1;
if qsave
    out_file = fullfile(pwd, 'rebuttal', 'Langevin_pin_pout', strcat('Langevin_no_noise_params_bad_traj_',fname) );
    save_figure_pdf(h3, 10, 8, out_file);
end
