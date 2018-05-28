%% Calculates the mutual information between pin and pout for given pin-pout map
% Plots mutual information against chosen parameter
clear variables
close all
clc
%% Parameters
gridsize = 15;
N=gridsize^2;
Con = 16;
K = 8;
a0 = 6;
hill = 2;
fileID = 'in_out_binary';
%noise = 0; %10^(-0.2); % all noise levels 10.^([-1 -0.8 -0.5 -0.2 0])
%tmax = 1000;
nsmpl = 100;
prec = 8;

outerVarlist = [0]; %10.^([-1 -0.8 -0.5 0]); % noise levels
%outerVarlist = [4 5 5.5 6.0 7 8]; %a0
%varlist = [4 5:0.1:6.0 7:10 15 20]; %a0
%varlist = [0.5 3 4 5 5.2 5.4 5.6 5.8 6 6.5 7 10]; %a0
varlist = 10.^[-1:0.1:0.5]; %10.^(-2:0.25:0.5); %noise
%varlist = 10.^(sort([-0.5:0.5:0.5 -0.25 0.25 -2 -1]));
Ivals = zeros(numel(outerVarlist), numel(varlist));
I_scaled = zeros(numel(outerVarlist), numel(varlist));
%%
for i=1:numel(outerVarlist)
   %noise = outerVarlist(i);   
   for irun=1:numel(varlist)
        %a0=varlist(irun);
        noise = varlist(irun);
        % correct file naming
        if (a0 >= 5 && a0 <= 6)
            appendix = sprintf('nsmpl100_run%d', i);
        elseif (a0 > 10)
            appendix = 'nsmpl100';
        else
            appendix = 'nsmpl1000';
        end
        %noise = varlist(irun);

        % Load pin-pout data
        % without noise
        %path = fullfile(pwd, 'figures', 'pin_pout', 'data_binary');
        %path = 'H:\My Documents\Multicellular automaton\finite_Hill\figures\pin_pout\data_binary';
        %fname_str = strrep(sprintf('pin_pout_N%d_Con%d_K%d_a0_%.2f_hill%.2f_prec%d_%s_%s',...
        %    N, Con, K, a0, hill, prec, fileID, appendix ), '.','p');
        %fname = fullfile(path, fname_str);
        %load(fname, 'prob', 'I_av');

        % with noise
        %path = 'H:\My Documents\Multicellular automaton\finite_Hill\figures\pin_pout\data';
        %path = 'H:\My Documents\Multicellular automaton\finite_Hill\data\pin_pout\noise';
        path = 'H:\My Documents\Multicellular automaton\finite_Hill\figures\pin_pout\data_binary';
        %fname_str = strrep(sprintf('pin_pout_noise_%.3f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_tmax%d_nsmpl%d_%s',...
        %    noise, N, Con, K, a0, hill, tmax, nsmpl, fileID), '.', 'p'); % var = a0
        %fname_str = strrep(sprintf('pin_pout_noise_%.4f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_%s', ...
        %    noise, N, Con, K, a0, hill, fileID), '.', 'p'); % var = a0
        %fname_str = strrep(sprintf('pin_pout_noise_%.3f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_%s', ...
        %    noise, N, Con, K, a0, hill, fileID), '.', 'p'); % var = noise
        fname_str = strrep(sprintf('pin_pout_N%d_Con%d_K%d_a0_%.2f_hill%.2f_noise%.2f_prec%d_%s_nsmpl%d_run1', ...
            N, Con, K, a0, hill, noise, prec, fileID, nsmpl), '.', 'p'); % var = noise
        %[fname, path, ~] = uigetfile(path);
        %load(fullfile(path, strcat(fname, '.mat')));
        load(fullfile(path, fname_str));
        %close all;

        % Compute mutual information
        % Compute S[ P(p_out) ]
        prob_pout = sum(prob, 2)/(N+1); % P(p_out)
        idx = prob_pout>0;
        S_pout = - sum(prob_pout(idx).*log2(prob_pout(idx)));

        % Compute S[ P(p_out | p_in) ]
        P_pin = 1/(N+1);
        S_pout_pin = zeros(N+1, 1);
        for k=1:N
            P_pout_pin = prob(:, k);
            idx2 = P_pout_pin > 0;
            S_pout_pin(k) = - sum(P_pout_pin(idx2).*log2(P_pout_pin(idx2)));
        end

        % Compute I(pin, pout)
        thisI = S_pout - sum(P_pin.*S_pout_pin);
        disp(irun); disp(thisI);

        Ivals(i, irun) = thisI;
        I_scaled(i, irun) = thisI/log2(N+1); % I rescaled to lie between 0 and 1
        %figure();
        %plot(prob_pout);
   end
end

%% Plot I against a0
set(0, 'defaulttextinterpreter', 'latex');
h1=figure(1);
%plot(varlist, Ivals, '-o', 'LineWidth', 2);
semilogx(repmat(varlist, numel(outerVarlist), 1)', Ivals', 'o--', 'LineWidth', 2);
hold on
xlabel('noise strength $$\alpha$$');
%xlabel('$$a_0$$');
ylabel('$$I(f_{in}, f_{out})$$');
set(gca, 'FontSize', 24);
%ylim([0 4]);

%{
h2=figure(2);
plot(varlist, I_scaled, '-o', 'LineWidth', 2);
%xlabel('$$a_0$$');
%semilogx(varlist, I_scaled, '-o', 'LineWidth', 2);
xlabel('noise strength $$\alpha$$');
ylabel('$$I(p_{in}, p_{out})/I_{max}$$');
set(gca, 'FontSize', 24);
ylim([0 1]);
%}
%%
% averaged
%{
h3=figure(3);
Ivals_m1 = mean(Ivals, 1);
Ivals_m2 = std(Ivals, 1)/sqrt(size(Ivals, 1));
%plot(varlist, Ivals_m1, '--o', 'LineWidth', 2);
errorbar(varlist, Ivals_m1, Ivals_m2, '--', 'LineWidth', 2);
%xlabel('noise strength $$\alpha$$');
xlabel('$$a_0$$');
ylabel('$$I(f_{in}, f_{out})$$');
set(gca, 'FontSize', 24);
%ylim([0 4]);
%}
%%
%{
h1=figure(1);
semilogx(varlist, Ivals, '-o', 'LineWidth', 2);
xlabel('noise strength $$\alpha$$');
ylabel('$$I(p_{in}, p_{out})$$');
set(gca, 'FontSize', 24);
%}
%% Temp plot (to add another noise level)
%{
Ivals = zeros(numel(varlist)-2, 1);
for irun=1:numel(varlist)-2
        a0=varlist(irun+2);
        %noise = varlist(irun);

        % Load pin-pout data
        % without noise
        %fname_str = strrep(sprintf('pin_pout_N%d_Con%d_K%d_a0_%d_hill%.2f_%s',...
        %    N, Con, K, a0*10, hill, fileID), '.','p');
        %fname = fullfile(pwd, 'figures', 'pin_pout', 'data', fname_str);
        %load(fname, 'prob');

        % with noise
        %path = 'H:\My Documents\Multicellular automaton\finite_Hill\figures\pin_pout\data';
        path = 'H:\My Documents\Multicellular automaton\finite_Hill\data\pin_pout\noise';
        fname_str = strrep(sprintf('pin_pout_noise_%.3f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_tmax%d_nsmpl%d_%s',...
            noise, N, Con, K, a0, hill, tmax, nsmpl, fileID), '.', 'p'); % var = a0
        %fname_str = strrep(sprintf('pin_pout_noise_%.4f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_%s', ...
        %    noise, N, Con, K, a0, hill, fileID), '.', 'p'); % var = a0
        %fname_str = strrep(sprintf('pin_pout_noise_%.3f_N%d_Con%d_K%d_a0_%.2f_hill%.2f_%s', ...
        %    noise, N, Con, K, a0, hill, fileID), '.', 'p'); % var = noise
        %[fname, path, ~] = uigetfile(path);
        %load(fullfile(path, strcat(fname, '.mat')));
        load(fullfile(path, fname_str));
        %close all;

        % Compute mutual information
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
        thisI = S_pout - sum(P_pin.*S_pout_pin);
        disp(irun); disp(thisI);

        Ivals(irun) = thisI;
        I_scaled(irun) = thisI/log2(N+1); % I rescaled to lie between 0 and 1
        %figure();
        %plot(prob_pout);
end
 %% Plot I against a0
h1=figure(1);
hold on
plot(varlist(3:end), Ivals, '-o', 'LineWidth', 2);
%xlabel('$$a_0$$');
%semilogx(varlist, Ivals, '-o', 'LineWidth', 2);
%xlabel('noise strength $$\alpha$$');
xlabel('$$a_0$$');
ylabel('$$I(p_{in}, p_{out})$$');
set(gca, 'FontSize', 24);
%ylim([0 4]);
%}    
%%
%{
varlist2 = [4 5:0.2:6 6.5 7];
Ivals2 = varlist2;
for irun=1:numel(varlist2)
        a0=varlist2(irun);
        
        % without noise
        fname_str = strrep(sprintf('pin_pout_N%d_Con%d_K%d_a0_%d_hill%.2f_%s',...
            N, Con, K, a0*10, hill, fileID), '.','p');
        fname = fullfile(pwd, 'figures', 'pin_pout', 'data', fname_str);
        load(fname, 'prob');
        % Compute mutual information
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
        thisI = S_pout - sum(P_pin.*S_pout_pin);
        disp(irun); disp(thisI);

        Ivals2(irun) = thisI;
        I_scaled(irun) = thisI/log2(N+1); % I rescaled to lie between 0 and 1
        %figure();
        %plot(prob_pout);
end
 %% Plot I against a0
h1=figure(1);
hold on
plot(varlist2, Ivals2, '-ok', 'LineWidth', 2);
%xlabel('$$a_0$$');
%semilogx(varlist, Ivals, '-o', 'LineWidth', 2);
%xlabel('noise strength $$\alpha$$');
xlabel('$$a_0$$');
ylabel('$$I(p_{in}, p_{out})$$');
set(gca, 'FontSize', 24);
%ylim([0 4]);
%}   
%% save figures
%
outdir = fullfile('H:\My Documents\Multicellular automaton\temp');
qsave = 1;
if qsave
    %fname_str = strrep(sprintf('Information_pin-pout_N%d_a0%d_K_%d_Con_%d_hill%.2f_nonoise',...
    %    N, 10*a0, K, Con, hill), '.', 'p');
    %fname_str = strrep(sprintf('Information_pin-pout_N%d_K_%d_Con_%d_hill%.2f_noise%.2f_vs_a0_%.1fto%.1f_all',...
    %    N, K, Con, hill, noise, varlist(1), varlist(end)), '.', 'p');
    fname_str = strrep(sprintf('Information_pin-pout_N%d_a0_%.2f_K_%d_Con_%d_hill%.2f_vs_noise_%s',...
        N, a0, K, Con, hill, fileID), '.', 'p');
    %fname_str = strrep(sprintf('Information_pin-pout_N%d_a0_%.2f_K_%d_Con_%d_hill%.2f_vs_noise_tmax%d_nsmpl%d',...
    %    N, a0, K, Con, hill, tmax, nsmpl), '.', 'p');
    fname = fullfile(outdir, fname_str); %filename
    save_figure_pdf(h1, 10, 8, fname);
    save_figure_eps(h1, 10, 8, fname);
end
%}
%{
qsave = 0;
if qsave
    %fname_str = strrep(sprintf('Information_rescaled_pin-pout_N%d_a0%d_K_%d_Con_%d_hilll%.2f',...
    %    N, 10*a0, K, Con, hill), '.', 'p');
    fname_str = strrep(sprintf('Information_pin-pout_N%d_K_%d_Con_%d_hill%.2f_noise%.2f_vs_a0_%.1fto%.1f_Imax',...
        N, K, Con, hill, noise, varlist(1), varlist(end)), '.', 'p');
    %fname_str = strrep(sprintf('Information_pin-pout_N%d_a0_%.2f_K_%d_Con_%d_hill%.2f_vs_noise_Imax',...
    %    N, a0, K, Con, hill), '.', 'p');
    %fname_str = strrep(sprintf('Information_pin-pout_N%d_a0_%.2f_K_%d_Con_%d_hill%.2f_vs_noise_tmax%d_nsmpl%d_Imax',...
    %    N, a0, K, Con, hill, tmax, nsmpl), '.', 'p');
    fname = fullfile(outdir, fname_str); %filename
    save_figure_pdf(h2, 10, 8, fname);
    save_figure_eps(h2, 10, 8, fname);
end
%}