% Time evolution of a system with visualization of the dynamics without
% noise showing the count of nearest neighbors that are ON
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

% Parameters 
Con = [30 30]; %column j = type j 
Coff = [1 1];
K = [5 5];
p = [0.2 0.2];
rcell = [0.2 0.2];
Rcell = rcell*a0;
Mcomm = [1 0; 0 1]; % Communication matrix, type i reacts to type j iff M_ij=1
Mcomm = logical(Mcomm);
Mcomm_str = sprintf('Mcomm_%d_%d_%d_%d', Mcomm(1,1), Mcomm(1,2), Mcomm(2,1), Mcomm(2,2));
tmax = 50;
%% Load saved trajectory
%
i = 1;
tf = 11;
fname_str = strrep(...
    sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2fa0_Rcell2_%.2fa0_t%d_tmax%d', ...
    N1, N2, p(1), p(2), a0, K(1), K(2), Con(1), Con(2), rcell(1), rcell(2), tf, tmax), '.', 'p');
fname = fullfile(pwd, 'figures' , 'time_evolution_examples',...
    strcat(fname_str,'-v',int2str(i), '.mat'));
load(fname);

%
frames = struct('cdata',[],'colormap',[]);
hin = figure();
for tau=1:t+1
    cells = cells_hist{tau};
    update_cell_figure(hin, pos, a0, cells, cell_type, tau);
    frames(tau) = getframe(gcf);
end

% Save video
fname_str = strrep(...
    sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2fa0_Rcell2_%.2fa0_t%d_tmax%d', ...
    N1, N2, p(1), p(2), a0, K(1), K(2), Con(1), Con(2), rcell(1), rcell(2), tf, tmax), '.', 'p');
fname = fullfile(pwd, 'figures' , 'time_evolution_examples', ...
    strcat(fname_str,'-v',int2str(i)));
myVideo = VideoWriter(fname, 'Uncompressed AVI');
myVideo.FrameRate = 1;  % Default 30
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);
%}
%%
%--------------------------------------------------------------------------
% Distance and position
[dist, pos] = init_dist_hex(gridsize, gridsize);

% cell types
cell_type = zeros(N,1);
idx1 = randperm(N,N1);
idx2 = setdiff(1:N, idx1);
cell_type(idx2) = 1;
% extra randomization of cell positions by requiring that cell_type also
% has a Moran's I value of around 0
[cell_type, ~] = generate_I_new(cell_type, 0, 0.01, dist, a0);
idx2 = find(cell_type);
idx1 = setdiff(1:N, idx2);

% Matrix of cell reading
M = zeros(size(dist));
    % same type
M(idx1, idx1) = Mcomm(1,1)*sinh(Rcell(1))./(a0*dist(idx1, idx1)).*exp(Rcell(1)-a0*dist(idx1, idx1));
M(sub2ind(size(M), idx1, idx1)) = Mcomm(1,1)*1; % 1: normal, 0: no self communication

M(idx2, idx2) = Mcomm(2,2)*sinh(Rcell(2))./(a0*dist(idx2, idx2)).*exp(Rcell(2)-a0*dist(idx2, idx2));
M(sub2ind(size(M), idx2, idx2)) = Mcomm(2,2)*1; % 1: normal, 0: no self communication
    % different type
M(idx1, idx2) = Mcomm(1,2)*Rcell(2)/(Rcell(1))*sinh(Rcell(1))./(a0*dist(idx1, idx2)).*exp(Rcell(2)-a0*dist(idx1, idx2)); % conc. cell of type 1 senses due to cells of type 2
M(idx2, idx1) = Mcomm(2,1)*Rcell(1)/(Rcell(2))*sinh(Rcell(2))./(a0*dist(idx2, idx1)).*exp(Rcell(1)-a0*dist(idx2, idx1));
    % self-communication
%M(sub2ind(size(M), 1:N, 1:N)) = 1; % 1: normal, 0: no self communication

% Exclude self-influence (only for calculation, not in simulations)
M(sub2ind(size(M), 1:N, 1:N)) = 0;

% Interaction strengths
f11 = sum(sum(M(idx1, idx1)))/N;
f12 = sum(sum(M(idx1, idx2)))/N; 
f21 = sum(sum(M(idx2, idx1)))/N; %not equal to f12 if cell radii different
f22 = sum(sum(M(idx2, idx2)))/N;

g11 = sum(sum(M(idx1, idx1)))/N;
g12 = sum(sum(M(idx1, idx2)))/N; 
g21 = sum(sum(M(idx2, idx1)))/N; %not equal to f12 if cell radii different
g22 = sum(sum(M(idx2, idx2)))/N;

% Restore self-influence
M(sub2ind(size(M), 1:N, 1:N)) = 1;

%%
% initialize ON cells
cells = zeros(N,1);
cells(idx1(randperm(N1,round(p(1)*N1)))) = 1;
cells(idx2(randperm(N2,round(p(2)*N2)))) = 1;

% initialize figure
hin = figure(1);
t = 0;
I = [];
I_11 = [];
I_22 = [];
I_12 = [];
theta_11 = [];
theta_22 = [];
theta_12 = [];
Non1 = [];
Non2 = [];
h = [];
cells_hist = {};

update_cell_figure(hin, pos, a0, cells, cell_type, t);
cells_hist{end+1} = cells;
Non1(end+1) = sum(cells(idx1));
Non2(end+1) = sum(cells(idx2));
I(end+1) = moranI(cells, a0*dist);
[theta_11(end+1), theta_12(end+1), theta_22(end+1)]...
    = moranI_twocelltypes(cells, idx1, idx2, a0*dist);
[cells_out, changed, h(end+1)] = update_cells_twocells(cells, idx1, idx2, M, Con, K);
while changed && t<tmax
    t = t+1;
    %k = waitforbuttonpress;
    update_cell_figure(hin, pos, a0, cells_out, cell_type, t);
    cells_hist{end+1} = cells_out;
    I(end+1) = moranI(cells_out, a0*dist);
    [theta_11(end+1), theta_12(end+1), theta_22(end+1)]...
        = moranI_twocelltypes(cells_out, idx1, idx2, a0*dist);
    Non1(end+1) = sum(cells_out(idx1));
    Non2(end+1) = sum(cells_out(idx2));
    cells = cells_out;
    [cells_out, changed, h(end+1)] = update_cells_twocells(cells, idx1, idx2, M, Con, K);
end
% Alternative definitions
p1 = Non1/N1;
p2 = Non2/N2;
I_12_v2 = (theta_12 - f12*(2*p1 - 1).*(2*p2 - 1))./(4*f12*sqrt(p1.*(1-p1).*p2.*(1-p2)) );
%}

%% Plot hamiltonian
%{
h2 = figure(2);
plot(0:t, h/N, 'o-');
xlabel('time');
ylabel('pseudo-energy h');
%}
%% Plot spatial index I
%{
h3 = figure(3);
plot(0:t, I, 'o-');
xlabel('time');
ylabel('spatial index I');
ylim([-0.2 1]);
%}
%% Plot p1, p2
h4 = figure(4);
hold on
plot(0:t, p1, 'ro-', 'LineWidth', 1.5);
plot(0:t, p2, 'bo-', 'LineWidth', 1.5);
xlabel('time');
ylabel('fraction ON cells');
legend({'p_1', 'p_2'}, 'Location', 'se');
ylim([0 1]);
set(gca, 'FontSize', 24);

qsave = 1;
if qsave 
    fname_str = strrep(...
        sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2fa0_Rcell2_%.2fa0_t%d_tmax%d_trajectory',...
        N1, N2, p(1), p(2), a0, K(1), K(2), Con(1), Con(2), rcell(1), rcell(2), t, tmax), '.', 'p');
    i = 1;
    fname = fullfile(pwd, 'figures' , 'time_evolution_examples', ...
        strcat(fname_str,'-v',int2str(i), '.pdf'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd, 'figures' , 'time_evolution_examples', ...
          strcat(fname_str,'-v',int2str(i), '.pdf'));
    end
    save_figure_pdf(h4, 10, 8, fname);
end
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
%
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
close all
qsave = 1;
if qsave
    fname_str = strrep(...
        sprintf('N1_%d_N2_%d_p1_%.2f_p2_%.2f_a0_%.2f_K1_%d_K2_%d_Con1_%d_Con2_%d_Rcell1_%.2fa0_Rcell2_%.2fa0_t%d_tmax%d', ...
        N1, N2, p(1), p(2), a0, K(1), K(2), Con(1), Con(2), rcell(1), rcell(2), t, tmax), '.', 'p');
    i = 1;
    %fname = fullfile('K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution', fname_str);
    %fname = fullfile(pwd, 'dynamics', 'twocelltypes', ...
    %    strcat(fname_str,'-v',int2str(i),'.mat'));
    fname = fullfile(pwd, 'figures' , 'time_evolution_examples', ...
        strcat(fname_str,'-v',int2str(i), '.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd, 'figures' , 'time_evolution_examples', ...
            strcat(fname_str,'-v',int2str(i), '.mat')); 
        %fname = fullfile(pwd, 'dynamics', 'twocelltypes', ...
            %strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    save(fname)
end