%% Calculates the mutual information between pin and pout for given pin-pout map
clear all
close all
clc

%% Parameters
gridsize = 11;
N=gridsize^2;
Con = 5;
K = 10;
a0 = 0.5;

%% Load pin-pout data
fname_str = sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d', N, Con, K, a0*10);
fname = fullfile(pwd, 'data', 'pin_pout', 'analytical', fname_str);
load(fname, 'prob');

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
set(0, 'defaulttextinterpreter', 'latex');
p = (0:N)/N;
h1=figure(1);
im_fig = imagesc(p, p, prob);
c = colorbar;
set(gca, 'YDir', 'Normal','FontSize', 24)
xlabel('$$p_{in}$$');
ylabel('$$p_{out}$$');
title(sprintf('$$I=%.3f = %.3f I_{max}$$', I, I_scaled));
ylabel(c, '$$P(p_{out}|p_{in})$$', 'Interpreter', 'latex');
% set invisible parts where count is zero
set(im_fig, 'AlphaData', prob > 0);
%
% save figure
qsave = 1;
if qsave
    fname_str = sprintf('pin-pout_withI_N%d_a0%d_K_%d_Con_%d',...
        N, 10*a0, K, Con);
    fname = fullfile(pwd, 'figures', 'information', fname_str); %filename
    save_figure_pdf(h1, 10, 8, fname);
    save_figure_eps(h1, 10, 8, fname);
end