%% Calculates the mutual information between pin and pout for given pin-pout map
% Plots mutual information against multiple chosen parameters
clear all
close all
clc

%% Parameters
gridsize = 15;
N=gridsize^2;
Con = 15;
K = 6;
a0 = 1.5;
%hill = 2;
%fileID = 'montecarlo';

% var1 = [a0 K Con]
% var2 = noise
var1list = [0.5 15 8; 1.5 6 15; 1.5 20 15];
var2list = {(0:0.1:0.3)*15, [0 0.05 (0.1:0.1:0.3)]*6, [0 1 1.8 2]};

Ivals = cell(size(var1list, 1), 1);
%I_scaled = varlist;
for i=1:size(var1list, 1)
    a0 = var1list(i, 1);
    K = var1list(i, 2);
    Con = var1list(i, 3);
    var2 = var2list{i};
    Ivals_temp = [];
    for j=1:numel(var2)
        noise=var2(j);

        % Load pin-pout data
        fname_str = sprintf('pin_pout_N%d_Con_%d_K_%d_gz_%d_a0_%d_noise_%d',...
            N, Con, K, gridsize, a0*10, round(noise*10));
        fname = fullfile(pwd, 'data', 'pin_pout', 'noise', fname_str);
        load(fname, 'count');
        prob = transpose(count./repmat(sum(count,2),1,N+1));

        % Compute mutual information 
        % Compute S[ P(p_out) ]
        prob_pout = sum(prob, 2)/(N+1); % P(p_out)
        idx = prob_pout>0;
        S_pout = - sum(prob_pout(idx).*log2(prob_pout(idx)));

        % Compute S[ P(p_out | p_in) ]
        P_pin = 1/(N+1);
        S_pout_pin = zeros(N+1, 1);
        for j=1:N
            P_pout_pin = prob(:, j);
            idx2 = P_pout_pin > 0;
            S_pout_pin(j) = - sum(P_pout_pin(idx2).*log2(P_pout_pin(idx2)));
        end

        % Compute I(pin, pout)
        thisI = S_pout - sum(P_pin.*S_pout_pin);
        Ivals_temp(end+1) = thisI;
        %I_scaled(i) = thisI/log2(N+1); % I rescaled to lie between 0 and 1
        %figure();
        %plot(prob_pout);
        
    end
    Ivals{i} = Ivals_temp;
end


%% Plot I against a0
set(0, 'defaulttextinterpreter', 'latex');
h1=figure(1);
hold on
for i=1:size(var1list,1)
    var2_temp = var2list{i}/var1list(i, 2);
    plot(var2_temp, Ivals{i}, '-o', 'LineWidth', 2);
end
xlabel('$$\alpha/K$$');
ylabel('$$I(p_{in}, p_{out})$$');
set(gca, 'FontSize', 24);
hold off
%fundamental limit
%plot([varlist(1) varlist(end)], [1 1], 'r--');
%hold off
%{
h2=figure(2);
plot(varlist, I_scaled, '-o', 'LineWidth', 2);
xlabel('Noise strength $$\alpha$$');
ylabel('$$I(p_{in}, p_{out})/I_{max}$$');
set(gca, 'FontSize', 24);
%}
% save figures
qsave = 1;
if qsave
    fname_str = strrep(sprintf('Information_pin-pout_noise_N%d_multiple_set1',...
        N), '.', 'p');
    fname = fullfile(pwd, 'figures', 'Information', fname_str); %filename
    save_figure_pdf(h1, 10, 8, fname);
    save_figure_eps(h1, 10, 8, fname);
end
%{
qsave = 0;
if qsave
    fname_str = strrep(sprintf('Information_pin-pout_noise_Imax_N%d_a0%d_K_%d_Con_%d_noise_%d',...
        N, 10*a0, K, Con, varlist(1), varlist(end)), '.', 'p');
    fname = fullfile(pwd, 'figures', 'pin_pout', 'noise', fname_str); %filename
    save_figure_pdf(h2, 10, 8, fname);
    save_figure_eps(h2, 10, 8, fname);
end
%}