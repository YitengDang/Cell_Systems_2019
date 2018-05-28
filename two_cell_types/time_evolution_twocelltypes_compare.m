% Compares the dynamics of the CA with MC dynamics where we take into
% account theta
close all
clear all
warning off
set(0, 'defaulttextinterpreter', 'latex');
%--------------------------------------------------------------------------
% Set parameters of the system
gridsize = 15;
N = gridsize^2;
frac_type1 = 0.5;
N1 = round(frac_type1*N);
N2 = N - N1;
a0 = 1.5;
rcell = [0.2 0.2];
Rcell = rcell*a0;
Mcomm = [1 1; 1 1]; % Communication matrix, type i reacts to type j iff M_ij=1
Mcomm = logical(Mcomm);

% Circuit parameters 
Con = [25 25]; %column j = type j 
K = [10 10];
Coff = [1 1]; %default [1 1]

% Simulation parameters
p = [0.5 0.5]; %fraction of type-1 and type-2 cells which are ON
tmax = 1000;
%--------------------------------------------------------------------------
%% Load trajectory
tend = 23;
version = 1;
path = 'H:\My Documents\Multicellular automaton\two_cell_types\figures\time_evolution_examples';
fname_str = strrep(sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2fa0_Rcell2_%.2fa0_t%d_tmax%d-v%d',...
    N1, N2, p(1), p(2), a0, K(1), K(2), Con(1), Con(2), rcell(1), rcell(2), tend, tmax, version), '.', 'p');
load(fullfile(path, fname_str), 'cells_hist', 'cell_type', 'idx1', 'idx2', 'p1', 'p2', 'theta_11', 'theta_12', 'theta_22',...
    'f11', 'f12', 'f21', 'f22', 'g11', 'g12', 'g21', 'g22');

%% Dynamics with only p1, p2
p1_all = cell(1, tend+1);
p2_all = p1_all;
theta_all = cell(3, tend+1);
p_MC_out = zeros(1, tend+1);
pe = zeros(tend+1, 2);
f_mat = [f11 f12; f21 f22];
g_mat = [g11 g12; g21 g22];
cont_all = zeros(1, tend+1);
for i=1:tend+1
    p_in = [p1(i); p2(i)];
    theta_in = [theta_11(i) theta_12(i) theta_22(i)];
    [p_MC_out, theta_out,  pe(i, :), cont_all(i)] = ...
        update_monte_carlo_with_theta(N1, N2, p_in, theta_in, Con, Coff, K, f_mat, g_mat);
    p1_all{i} = [p1(i) p_MC_out(1)];
    p2_all{i} = [p2(i) p_MC_out(2)];
    theta_all{1,i} = theta_out(1);
    theta_all{2,i} = theta_out(2);
    theta_all{3,i} = theta_out(3);
end
%% Plot dynamics
h1 = figure(1);
hold on
plot(0:tend, p1, 'bo-');
plot(0:tend, p2, 'ro-');
% dashed lines --> MC evolution at given time step given the p, theta
for i=1:tend+1
    plot([i-1 i], p1_all{i}, 'bx--');
    plot([i-1 i], p2_all{i}, 'rx--');
end
legend({'p_1', 'p_2'});
xlabel('time');
ylabel('fraction ON cells');
ylim([0 1]);
set(gca, 'FontSize', 24);
set(h1, 'Units', 'Inches', 'Position', [1 1 10 8]);

% equilibrium probabilities
h2 = figure(2);
hold on
plot(0:tend, pe(:,1), 'b.-'); % termination probability
plot(0:tend, pe(:,2), 'r.-'); % termination probability
%plot(0:tend, terminate_all, 'kx--'); % termination probability
legend({'p_1', 'p_2', 'p_{eq}'}, 'Location', 'nw');
xlabel('time');
ylabel('eq. probability');
ylim([0 1]);
set(gca, 'FontSize', 24);
set(h2, 'Units', 'Inches', 'Position', [1 1 10 8]);


 %% Compare calculated average sensed concentrations
time = tend+1;
cells = cells_hist{time};
p = [p1(time); p2(time)];

% (1) from exact simulation
N = numel(cells);
Con_cells = zeros(N, 1);
K_cells = zeros(N, 1);
Con_cells(idx1) = Con(1);
Con_cells(idx2) = Con(2);
K_cells(idx1) = K(1);
K_cells(idx2) = K(2);

% Concentration in each cell
C0 = 1 + (Con_cells-1).*cells; % term in brackets of Eq. S9 p.S5

% Reading of each cell
[dist, ~] = init_dist_hex(gridsize, gridsize);
M = zeros(size(dist));
    % same type
M(idx1, idx1) = Mcomm(1,1)*sinh(Rcell(1))./(a0*dist(idx1, idx1)).*exp(Rcell(1)-a0*dist(idx1, idx1));
%M(sub2ind(size(M), idx1, idx1)) = Mcomm(1,1)*1; % 1: normal, 0: no self communication
M(idx2, idx2) = Mcomm(2,2)*sinh(Rcell(2))./(a0*dist(idx2, idx2)).*exp(Rcell(2)-a0*dist(idx2, idx2));
%M(sub2ind(size(M), idx2, idx2)) = Mcomm(2,2)*1; % 1: normal, 0: no self communication
    % different type
M(idx1, idx2) = Mcomm(1,2)*Rcell(2)/(Rcell(1))*sinh(Rcell(1))./(a0*dist(idx1, idx2)).*exp(Rcell(2)-a0*dist(idx1, idx2)); % conc. cell of type 1 senses due to cells of type 2
M(idx2, idx1) = Mcomm(2,1)*Rcell(1)/(Rcell(2))*sinh(Rcell(2))./(a0*dist(idx2, idx1)).*exp(Rcell(1)-a0*dist(idx2, idx1));
    % self-communication
M(sub2ind(size(M), 1:N, 1:N)) = 0; % 1: normal, 0: no self communication

Y_nei = M*C0; % Eq. S9 p.S5
Y_type1 = Y_nei(cell_type==0);
Y_type2 = Y_nei(cell_type==1);
%muon_exact = mean(Y(cells==1));
%muoff_exact = mean(Y(cells==0));
muon_exact = [mean(Y_type1((cells(cell_type==0))==1)); ...
    mean(Y_type2((cells(cell_type==1))==1))];
muoff_exact = [mean(Y_type1((cells(cell_type==0))==0)); ...
    mean(Y_type2((cells(cell_type==1))==0))];
fprintf('---------------------------- \n');
fprintf('muon exact, type 1 = %.3f \n', muon_exact(1));
fprintf('muon exact, type 2 = %.3f \n', muon_exact(2));
fprintf('muoff exact, type 1 = %.3f \n', muoff_exact(1));
fprintf('muoff exact, type 2 = %.3f \n', muoff_exact(2));
fprintf('---------------------------- \n');

% (2) from MC without spatial order
S = (Con'-Coff')./2.*(2.*p-1) + (Con'+Coff')/2;
%muon1 = sum(f_mat*S, 2);
muon1 = ([N/N1; N/N2].*f_mat*S);
%muoff1 = muon1;
fprintf('mu MC, type 1 = %.3f \n', muon1(1));
fprintf('mu MC, type 2 = %.3f \n', muon1(2));
fprintf('---------------------------- \n');

% (3) from MC with spatial order
N_12 = [N1; N2];
Non = round([N1; N2].*p);
theta = [theta_11(end) theta_12(end) theta_22(end)];
theta_mat = [theta(1) theta(2); theta(2) theta(3)]; 
xl = [N1/N; N2/N];
%xl = [1; 1];

Xlm_on = diag(1./(2.*p))*(theta_mat + repmat(2*p'-1, 2, 1));
Xlm_off = diag(1./(2.*(1-p)))*(-theta_mat + repmat(2*p'-1, 2, 1));

S_on = repmat((Con-Coff)/2, 2, 1).*Xlm_on + repmat((Con+Coff)/2, 2, 1);
S_off = repmat((Con-Coff)/2, 2, 1).*Xlm_off + repmat((Con+Coff)/2, 2, 1);

muon = sum(diag(1./xl)*f_mat.*S_on, 2);
muoff = sum(diag(1./xl)*f_mat.*S_off, 2);
fprintf('muon MC with theta, type 1 = %.3f \n', muon(1));
fprintf('muon MC with theta, type 2 = %.3f \n', muon(2));
fprintf('muoff MC with theta, type 1 = %.3f \n', muoff(1));
fprintf('muoff MC with theta, type 2 = %.3f \n', muoff(2));
fprintf('---------------------------- \n');

%% compare Ponon Poffoff
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

fprintf('Poffoff MC, type 1 = %.3f \n', Poffoff(1));
fprintf('Poffoff MC, type 2 = %.3f \n', Poffoff(2));
fprintf('Ponon MC, type 1 = %.3f \n', Ponon(1));
fprintf('Ponon MC, type 2 = %.3f \n', Ponon(2));
fprintf('---------------------------- \n');

Poffoff_exact = normcdf((K' - 1 - muoff_exact)./sigmaoff);
Ponon_exact = 1-normcdf( (K' - Con' - muon_exact)./sigmaon);
fprintf('Poffoff exact, type 1 = %.3f \n', Poffoff_exact(1));
fprintf('Poffoff exact, type 2 = %.3f \n', Poffoff_exact(2));
fprintf('Ponon exact, type 1 = %.3f \n', Ponon_exact(1));
fprintf('Ponon exact, type 2 = %.3f \n', Ponon_exact(2));
fprintf('---------------------------- \n');

%% Plots
%{
% hamiltonian
h2 = figure(2);
plot(0:t, h/N, 'o-');
xlabel('time');
ylabel('pseudo-energy h');
%% Plot spatial index I
h3 = figure(3);
plot(0:t, I, 'o-');
xlabel('time');
ylabel('spatial index I');
ylim([-0.2 1]);
%% Plot p1, p2
h4 = figure(4);
hold on
plot(0:t, Non1/N1, 'ro-');
plot(0:t, Non2/N2, 'bo-');
xlabel('time');
ylabel('fraction ON cells');
legend({'type 1', 'type 2'});
ylim([0 1]);
%% Plot I11, I22, I12
%{
h5 = figure(5);
hold on
plot(0:t, I_11, 'ro-');
plot(0:t, I_22, 'bo-');
plot(0:t, I_12, 'go-');
xlabel('time');
ylabel('spatial order I');
legend({'11', '22', '12'});
ylim([-0.5 1]);
%}
%% Plot theta_11, theta_22, theta_12
h6 = figure(6);
hold on
plot(0:t, theta_11, 'ro-');
plot(0:t, theta_22, 'bo-');
plot(0:t, theta_12, 'go-');
xlabel('time');
ylabel('spatial order theta');
legend({'11', '22', '12'});
ylim([-0.5 1]);
%}
%% Save trajectory
qsave = 0;
if qsave
    fname_str = strrep(...
        sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2f_Rcell2_%.2f_t%d', ...
        N1, N2, p(1), p(2), a0, K(1), K(2), Con(1), Con(2), Rcell(1), Rcell(2), t), '.', 'p');
    i = 1;
    %fname = fullfile('K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution', fname_str);
    fname = fullfile(pwd, 'dynamics', 'twocelltypes', ...
        strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
      fname = fullfile(pwd, 'dynamics', 'twocelltypes', ...
          strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    save(fname)
end