%% Tester code for loading a simulation file and using it 
clear all
close all
%%
folder = 'N:\tnw\BN\HY\Shared\Yiteng\two_signals\parameter_gradient_sinusoidal\sine_wave_gradient_K21_Ax_0p0_nx_1_Ay_0p0_ny_1';
fname_str = 'Parameter_gradient_K_2_1_square_wave_t_out_29_period_Inf-v1';
load(fullfile(folder, fname_str))

%%
M_int = save_consts_struct.M_int;
a0 = save_consts_struct.a0;
Rcell = save_consts_struct.rcell*a0;
Con = save_consts_struct.Con;
Coff = [1 1]; %save_consts_struct.Coff;
K = save_consts_struct.K;
lambda = [1 save_consts_struct.lambda12];
%hill = save_consts_struct.hill;
hill = 2;
noise = save_consts_struct.noise;

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