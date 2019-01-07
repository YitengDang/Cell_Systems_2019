%% Tester code for loading a simulation file and using it 
clear all
close all
%%
folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\Fig5-sample complex trajectory';
fname_str = 'sample_complex_trajectory';
load(fullfile(folder, fname_str))

%%
M_int = save_consts_struct.M_int;
a0 = save_consts_struct.a0;
Rcell = save_consts_struct.rcell*a0;
Con = save_consts_struct.Con;
Coff = [1 1]; %save_consts_struct.Coff;
K = save_consts_struct.K;
lambda = [1 save_consts_struct.lambda12];
hill = save_consts_struct.hill;
%hill = 2;
noise = save_consts_struct.noise;

%%
t0 = 0;
option = 1;
fig_pos = [1 1 8 4];
msg = plot_I_vs_t(cells_hist(1:51), t0, a0, distances, option, fig_pos);
box on

h = gcf;
folder = 'H:\My Documents\Multicellular automaton\paper_2_draft\figures\originals\Fig5-sample complex trajectory';
fname_str = 'sample_complex_trajectory_I_vs_t_tmax_20_v3';
fname = fullfile(folder, fname_str);
qsave = 1;
save_figure(h, 8, 4, fname, '.pdf', qsave)

%%
dist = distances;
cells = cells_hist{25};
%%
N = size(cells, 1);

% Account for self-influence
idx = dist>0;
% TO DO: vectorize / combine
M1 = ones(size(dist)); 
M1(idx) = sinh(Rcell)./(a0*dist(idx)/lambda(1))...
    .*exp((Rcell-a0*dist(idx))/lambda(1));
M2 = ones(size(dist)); 
M2(idx) = sinh(Rcell)./(a0*dist(idx)/lambda(2))...
    .*exp((Rcell-a0*dist(idx))/lambda(2));

% Concentration in each cell
C0 = Coff + (Con-Coff).*cells; 

% Reading of each cell
Y1 = M1*C0(:, 1); 
Y2 = M2*C0(:, 2);
Y = [Y1 Y2];

% Add noise to K
K_cells = K.*ones(2, 2, N);
dK = normrnd(0, noise, 2, 2, N);
K_cells = max(K_cells + dK, 1); % do not allow K < 1 = Coff

K_cells_1 = squeeze(K_cells(1,:,:))';
K_cells_2 = squeeze(K_cells(2,:,:))';

fX1 = 1./( 1 + ...
    ((Y./K_cells_1).*(1-M_int(1,:))/2).^hill + ... % case repression 
    ((K_cells_1./Y).*(1+M_int(1,:))/2).^hill ... % case activation
    ).*abs(M_int(1,:)) + ... 
    (1-abs(M_int(1,:))).*ones(N,2); % case no interaction
fX2 = 1./( 1 + ...
    ((Y./K_cells_2).*(1-M_int(2,:))/2).^hill + ... % case repression 
    ((K_cells_2./Y).*(1+M_int(2,:))/2).^hill ... % case activation
    ).*abs(M_int(2,:)) + ... 
    (1-abs(M_int(2,:))).*ones(N,2); % case no interaction

fX1_old = (Y.^hill.*(1+M_int(1,:))/2 + (squeeze(K_cells(1,:,:))'.^hill.*(1-M_int(1,:))/2))...
    ./(squeeze(K_cells(1,:,:))'.^hill+Y.^hill).*abs(M_int(1,:)) + (1-abs(M_int(1,:))).*ones(N,2);
fX2_old = (Y.^hill.*(1+M_int(2,:))/2 + (squeeze(K_cells(2,:,:))'.^hill.*(1-M_int(2,:))/2))...
    ./(squeeze(K_cells(2,:,:))'.^hill+Y.^hill).*abs(M_int(2,:)) + (1-abs(M_int(2,:))).*ones(N,2);

%% 
%all(fX2(:) >= 0 &  fX2(:) <= 1)
max( abs(fX1(:)-fX1_old(:)) )
max( abs(fX2(:)-fX2_old(:)) )