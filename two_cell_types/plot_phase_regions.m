%% Plots the different phases as regions in projected space at fixed values of other parameters
close all
clear all
warning off
set(0, 'defaulttextinterpreter', 'latex');

%--------------------------------------------------------------------------
gridsize = 15;
N = gridsize^2;
frac_type1 = 0.5;
N1 = round(frac_type1*N);
N2 = N - N1;

a0 = 1.5;
rcell = [0.2 0.2];
Rcell = rcell*a0;
%K = [2 5];
Mcomm = [1 0; 1 0];
Mcomm = logical(Mcomm);
subfolder = sprintf('Mcomm_%d_%d_%d_%d', Mcomm(1,1), Mcomm(1,2), Mcomm(2,1), Mcomm(2,2));
Klist = [2 5 10];

%--------------------------------------------------------------------------
for i1=1:numel(Klist)
    for i2=1:numel(Klist)
        K = [Klist(i1) Klist(i2)];
        %% Directly calculate phases
        % Load interaction strength
        ntrials = 10^3;
        version = 1;
        fname_str = strrep(...
            sprintf('N1_%d_N2_%d_frac1_%.2f_a0_%.2f_Rcell1_%.2fa0_Rcell2_%.2fa0_ntrials_%d',...
            N1, N2, frac_type1, a0, rcell(1), rcell(2), ntrials), '.', 'p');
        fname = fullfile(pwd, 'data', 'interaction_strengths',...
            strcat(fname_str,'-v', num2str(version)));
        load(fname);
        f11 = N/N1*Mcomm(1,1)*mean(f11); %normalization issue (follow notation of iScience Supplement)
        f12 = N/N1*Mcomm(1,2)*mean(f12);
        f21 = N/N2*Mcomm(2,1)*mean(f21);
        f22 = N/N2*Mcomm(2,2)*mean(f22);

        %% Check manually
        %{
        Con = [30 5];
        clc;

        % Check for autonomy
        % autonomy for type 1
        cond_aut1 = [Con(1)*Mcomm(1,1) + f11 + f12 > K(1) ; ...
            1*Mcomm(1,1) + f11*Con(1)+f12*Con(2) < K(1)];
        fprintf('Type 1 autonomous? %d \n', all(cond_aut1));
        fprintf('[ON remains ON | OFF remains OFF] = [%d | %d] \n', cond_aut1(1), cond_aut1(2));

        % autonomy for type 2
        cond_aut2 = [Con(2)*Mcomm(2,2) + f21 + f22 > K(2) ; ...
            1*Mcomm(2,2) + f21*Con(1)+f22*Con(2) < K(2)];
        fprintf('Type 2 autonomous? %d \n', all(cond_aut2));
        fprintf('[ON remains ON | OFF remains OFF] = [%d | %d] \n', cond_aut2(1), cond_aut2(2));

        % Check for all ON
        cond_allON = [Mcomm(1,1) + f11 + f12 > K(1); Mcomm(2,2) + f21 + f22 > K(2)];
        fprintf('Type 1 all ON? %d \n', cond_allON(1) );
        fprintf('Type 2 all ON? %d \n', cond_allON(2) );

        % Check for all OFF
        cond_allOFF = [Con(1)*(Mcomm(1,1) + f11) + Con(2)*f12 < K(1); Con(2)*(Mcomm(2,2) + f22) + f21*Con(1) < K(2)];
        fprintf('Type 1 all OFF? %d \n', cond_allOFF(1) );
        fprintf('Type 2 all OFF? %d \n', cond_allOFF(2) );
        %}
        %% Plot across phase range
        Con_range = 1:0.5:40;
        states_type1 = zeros(numel(Con_range), numel(Con_range));
        states_type2 = states_type1;

        for i=1:numel(Con_range)
            for j=1:numel(Con_range)
                Con = [Con_range(i) Con_range(j)];
                cond_aut1 = [Con(1)*Mcomm(1,1) + f11 + f12 > K(1); ...
                    1*Mcomm(1,1) + f11*Con(1)+f12*Con(2) < K(1)];
                cond_aut2 = [Con(2)*Mcomm(2,2) + f21 + f22 > K(2) ; ...
                    1 + f21*Con(1)+f22*Con(2) < K(2)];
                cond_allON = [1*Mcomm(1,1) + f11 + f12 > K(1); 1*Mcomm(2,2) + f21 + f22 > K(2)];
                cond_allOFF = [Con(1)*(1*Mcomm(1,1) + f11) + Con(2)*f12 < K(1); Con(2)*(1*Mcomm(2,2) + f22) + f21*Con(1) < K(2)];

                if (cond_aut1(1) && cond_aut1(2))
                    states_type1(i,j) = 3; % autonomy
                elseif (cond_aut1(2)+cond_allOFF(1)+cond_aut1(1)+cond_allON(1))==0
                    states_type1(i,j) = 4; % activation / deactivation
                else
                    states_type1(i,j) = 2*cond_aut1(2)-cond_allOFF(1) + 5*cond_aut1(1)+cond_allON(1);
                end

                if (cond_aut2(1) && cond_aut2(2))
                    states_type2(i,j) = 3; %autonomy
                elseif (cond_aut2(2)+cond_allOFF(2)+cond_aut2(1)+cond_allON(2))==0
                    states_type2(i,j) = 4; % activation / deactivation
                else
                    states_type2(i,j) = 2*cond_aut2(2)-cond_allOFF(2) + 5*cond_aut2(1)+cond_allON(2);
                end
            end
        end

        h1=figure(1);
        imagesc(Con_range, Con_range, states_type1');
        title(sprintf('Cell type 1, $$K_1=%d, K_2=%d$$', K(1), K(2)));
        xlabel('$$C_{ON, 1}$$');
        ylabel('$$C_{ON, 2}$$');
        set(gca, 'YDir', 'normal', 'FontSize', 24);
        c = colorbar('Ticks', 1:6,...
            'TickLabels', {'all OFF', 'OFF->OFF', 'Autonomy', 'Act.-deact.', 'ON->ON', 'all ON'});
        colormap('jet');
        caxis([1 6]);

        h2=figure(2);
        imagesc(Con_range, Con_range, states_type2');
        title(sprintf('Cell type 2, $$K_1=%d, K_2=%d$$', K(1), K(2)));
        xlabel('$$C_{ON, 1}$$');
        ylabel('$$C_{ON, 2}$$');
        set(gca, 'YDir', 'normal', 'FontSize', 24);
        c = colorbar('Ticks', 1:6,...
                 'TickLabels', {'all OFF', 'OFF->OFF', 'Autonomy', 'Act.-deact.', 'ON->ON', 'all ON'});
        colormap('jet');
        caxis([1 6]);

        %% Save figures
        qsave = 1;
        if qsave 
            fname_str = strrep(sprintf('N%d_a0_%.1f_N1_%d_K1_%d_K2_%d_%s',...
                N, a0, N1,  K(1), K(2), subfolder), '.', 'p');
            path = fullfile('H:\My Documents\Multicellular automaton\two_cell_types\figures\phase_maps', subfolder); 
            %path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\figures\time_evolution_twocelltypes\Mcomm_1_1_1_1';
            %path = 'L:\BN\HY\Shared\Yiteng\Multicellularity\data\time_evolution_twocelltypes';
            %path ='K:\bn\hy\Shared\Yiteng\Multicellularity\data_twocelltypes\time_evolution';
            fname1 = fullfile(path, strcat('CellType1_', fname_str));
            fname2 = fullfile(path, strcat('CellType2_', fname_str));
            %fname = fullfile(pwd, 'data', 'time_evolution', fname_str);
            save_figure_pdf(h1, 10*1.4, 8*1.4, fname1);    
            save_figure_pdf(h2, 10*1.4, 8*1.4, fname2);
        end
        close all;
    end
end
%% Plot phase boundaries
%{
% Load interaction strength
fname_str = strrep(...
    sprintf('N1_%d_N2_%d_frac1_%.2f_a0_%.2fto%.2f_Rcell1_%.2fa0_Rcell2_%.2fa0_ntrials_%d',...
    N1, N2, frac_type1, a0_range(1), a0_range(end), rcell(1), rcell(2), ntrials), '.', 'p');
fname = fullfile(pwd, 'data', 'interaction_strengths',...
    strcat(fname_str,'-v1'));
load(fname);


idx = 20;
a0 = a0_range(idx);
f11 = f11_avg(idx);
f12 = f12_avg(idx);
f21 = f21_avg(idx);
f22 = f22_avg(idx);

K1_range = 1:40;
Con2_range = [2 5 10 20];
for i=1:numel(Con2_range)
    Con2 = Con2_range(i);
    % type 1
    %Con2 = 10; % fix Con
    figure(i);
    plot([1+f11+f12 1+f11+f12], [1 40], 'g-'); % all ON
    hold on
    plot(K1_range, (K1_range - f12*Con2)./(1+f12), 'r-'); % all OFF
    plot(K1_range, K1_range-f11-f12, 'k--'); % ON cannot turn OFF
    plot(K1_range, (K1_range-1-f12*Con2)/f11, 'k-.'); % Off cannot turn ON
    xlabel('$$K_1$$');
    ylabel('$$C_{ON,1}$$'); 
    set(gca, 'FontSize', 24);
    title(sprintf('a0 = %.2f, Con2 = %.1f', a0, Con2));
    xlim([1 40]);
    ylim([1 40]);
end
%}
%% Calculate phase for specific set of parameters
%{
idx = 15;

gridsize = 15;
N = gridsize^2;
a0 = a0_range(idx);
rcell = [0.2 0.2];
Rcell = rcell*a0;
Con = [20 20];
K = [10 10];

% plot data for specific set of parameters

Con2 = Con(2);
figure(5);
plot([1+f11+f12 1+f11+f12], [1 40], 'g-'); % all ON
hold on
plot(K1_range, (K1_range - f12*Con2)./(1+f12), 'r-'); % all OFF
plot(K1_range, K1_range-f11-f12, 'k--'); % ON cannot turn OFF above line
plot(K1_range, (K1_range-1-f12*Con2)/f11, 'k-.'); % Off cannot turn ON below line
plot(K(1), Con(1), 'rx');
xlabel('$$K_1$$');
ylabel('$$C_{ON,1}$$'); 
set(gca, 'FontSize', 24);
title(sprintf('a0 = %.2f, Con2 = %.1f', a0, Con2));
xlim([1 40]);
ylim([1 40]);
%}