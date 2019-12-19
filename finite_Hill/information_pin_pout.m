%% Calculates the mutual information between pin and pout for given pin-pout map
% Plots single map with mutual information as title
clear all
close all
clc
set(0, 'defaulttextinterpreter', 'tex');
%% Parameters
gz = 15;
N=gz^2;
Con = 8;
K = 15;
a0 = 0.5; %5.6;
hill = Inf;
noise = 3.0; % 10^(-0.75);

%% Load pin-pout data
%data_path = 'H:\My Documents\Multicellular automaton\data\pin_pout\noise';
%data_path = 'H:\My Documents\Multicellular automaton\finite_Hill\figures\pin_pout\data';
%[fname, path, ~] = uigetfile(data_path);
%load(fullfile(path,fname));

%{
%fname_str = strrep(sprintf('pin_pout_N%d_Con%d_K%d_a0_%d_hill%.2f', N, Con, K, a0*10, hill),...
%    '.','p');
fname_str = strrep(sprintf('pin_pout_noise_%.4f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_montecarlo',...
     noise, N, Con, K, a0, hill), '.','p');
fname = fullfile(pwd, 'figures', 'pin_pout', 'data', fname_str);
load(fname, 'prob');
%}

%--- infinite hill, vs K, Con, from M drive --- 
%{
folder = 'M:\tnw\bn\hy\Shared\Yiteng\Multicellularity\data_archived\pin_pout\pin-pout_analytical_N225_a0_1p5_K1to30_Con1to30';
fname_str = sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d',N, Con, K, a0*10);
load(fullfile(folder, fname_str));
%}
folder = 'H:\My Documents\Multicellular automaton\data\main\pin_pout\noise';
fname_str = sprintf('pin_pout_N%d_Con_%d_K_%d_gz_%d_a0_%d_noise_%d',...
    N, Con, K, gz, a0*10, noise*10);
load(fullfile(folder, fname_str));
nsmpl = 1000;
prob = count'/nsmpl;

%--- finite Hill, vs a0 and noise -------------
%{
nsmpl = 100;
folder = 'H:\My Documents\Multicellular automaton\figures\finite_Hill\pin_pout\data_binary\vs_noise';
%{
fname_str = strrep(sprintf(...
    'pin_pout_N%d_Con%d_K%d_a0_%.2f_hill%.2f_prec8_in_out_binary_nsmpl%d',...
    N, Con, K, a0, hill, nsmpl), '.', 'p');
%}
fname_str = strrep(sprintf(...
    'pin_pout_N%d_Con%d_K%d_a0_%.2f_hill%.2f_noise%.2f_prec8_in_out_binary_nsmpl%d_run1',...
    N, Con, K, a0, hill, noise, nsmpl), '.', 'p');
load(fullfile(folder, fname_str));
%}
%% Compute mutual information
% Extreme cases
%prob = diag(ones(N+1,1)); % autonomy
%prob = zeros(N+1); prob(1,:)=1; % all OFF 

% Compute S[ P(p_out) ]
prob_pout = sum(prob, 2)/(N+1); % P(p_out)
idx = prob_pout>0;
S_pout = - sum(prob_pout(idx).*log2(prob_pout(idx)));

% Compute S[ P(p_out | p_in) ]
P_pin = 1/(N+1);
S_pout_pin = zeros(N+1, 1);
for i=1:N
    P_pout_pin = prob(:, i);
    idx2 = P_pout_pin > 0;
    S_pout_pin(i) = - sum(P_pout_pin(idx2).*log2(P_pout_pin(idx2)));
end
 
% Compute I(pin, pout)
I = S_pout - sum(P_pin.*S_pout_pin);
I_scaled = I/log2(N+1); % I rescaled to lie between 0 and 1
%figure();
%plot(prob_pout);

%% Plot pin-pout map with I
p = (0:N)/N;
h1=figure(1);
im_fig = imagesc(p, p, prob);
c = colorbar;
colormap('winter');
caxis([0 1]);
ylim([-0.01 1.01]);
set(gca, 'YDir', 'Normal','FontSize', 32)
xlabel('p_{in}');
ylabel('p_{out}');
%xlabel('f_{in}');
%ylabel('f_{out}');
%title(sprintf('I=%.3f = %.3f I_{max}', I, I_scaled));
title(sprintf('I = %.3f', I));
ylabel(c, 'P(p_{out}|p_{in})', 'Interpreter', 'tex');
% set invisible parts where count is zero
set(im_fig, 'AlphaData', prob > 0);

% save figure
qsave = 1;
if qsave    
    %fname_str = strrep(sprintf('pin-pout_withI_N%d_a0_%.1f_K_%d_Con_%d',...
    %    N, a0, K, Con), '.', 'p');
    fname_str = strrep(sprintf('pin-pout_withI_N%d_a0_%.1f_K_%d_Con_%d_hill_%.2f_noise_%.2f',...
        N, a0, K, Con, hill, noise), '.', 'p');
    folder = 'H:\My Documents\Thesis\Information theory\Fig2';
    fname = fullfile(folder, fname_str);
    %fname = fullfile(pwd, 'figures', 'information', fname_str); %filename
    
    save_figure_pdf(h1, 10, 8, fname);
    %save_figure_eps(h1, 10, 8, fname);
end