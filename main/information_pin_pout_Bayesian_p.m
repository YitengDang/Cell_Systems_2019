%% Calculates P(pin|pout) given P(pout|pin) using Bayes' rule
clear all
close all
clc

%% Parameters
gridsize = 15;
N=gridsize^2;
Con = 14;
K = 17;
a0 = 1.5;

%% Load pin-pout data
fname_str = sprintf('pin_pout_N%d_Con_%d_K_%d_a0_%d', N, Con, K, a0*10);
fname = fullfile(pwd, 'data', 'pin_pout', 'analytical', fname_str);
load(fname, 'prob');

%% Plot pin-pout map
set(0, 'defaulttextinterpreter', 'latex');
p = (0:N)/N;
h1=figure(1);
im_fig = imagesc(p, p, prob);
c = colorbar;
set(gca, 'YDir', 'Normal','FontSize', 24)
xlabel('$$p_{in}$$');
ylabel('$$p_{out}$$');
ylabel(c, '$$P(p_{out}|p_{in})$$', 'Interpreter', 'latex');
% set invisible parts where count is zero
set(im_fig, 'AlphaData', prob > 0);
xlim([0 1]);

%%
w_tot = sum(sum(prob));
P_pin = sum(prob, 1)/w_tot;
P_pout = sum(prob, 2)/w_tot;

%figure();
%plot(p, P_pout);

prob_Bayes = zeros(N+1,N+1);
for i=1:N+1
    for j=1:N+1
        prob_Bayes(i,j) = P_pin(i)/P_pout(j)*prob(j,i);        
    end
end
%%
h2=figure(2);
hold on;
im_fig = imagesc(p, p, prob_Bayes);
c = colorbar;
set(gca, 'YDir', 'Normal','FontSize', 24)
xlabel('$$p_{out}$$');
ylabel('$$p_{in}$$');
ylabel(c, '$$P(p_{in}|p_{out})$$', 'Interpreter', 'latex');
% set invisible parts where count is zero
set(im_fig, 'AlphaData', prob_Bayes > 0);
xlim([0 1]);
ylim([0 1]);

qsave = 1;
if qsave
    fname_str = sprintf('pout-pin_Bayes_N%d_a0_%d_K_%d_Con_%d',...
        N, 10*a0, K, Con);
    fname = fullfile(pwd, 'figures', 'information', fname_str); %filename
    save_figure_pdf(h2, 10, 8, fname);
    save_figure_eps(h2, 10, 8, fname);
end
%%
h3=figure(3);
prob_mean = p*prob_Bayes;
prob_n2 = (p.^2)*prob_Bayes;
prob_var = prob_n2 - prob_mean.^2; 
%plot(p, prob_mean, 'x-');
errorbar(p, prob_mean, sqrt(prob_var), 'rx-');
set(gca, 'FontSize', 24)
xlabel('$$p_{out}$$');
ylabel('Expected $$p_{in}$$');

qsave = 1;
if qsave
    fname_str = sprintf('pout-pin_Bayes_N%d_a0_%d_K_%d_Con_%d_first_moments',...
        N, 10*a0, K, Con);
    fname = fullfile(pwd, 'figures', 'information', fname_str); %filename
    save_figure_pdf(h3, 10, 8, fname);
    save_figure_eps(h3, 10, 8, fname);
end

%% Plot P(pin|pout) for fixed pout
h4=figure(4);
nout = 60;
plot(p, prob_Bayes(:,nout+1))
xlabel('$$p_{in}$$');
ylabel('$$P(p_{in}|p_{out})$$');
