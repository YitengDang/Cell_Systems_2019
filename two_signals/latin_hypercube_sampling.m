%% Implements Latin hypercube sampling for a region of K, Con we wish to examine
clear all 
close all
% take square in K, Con space with lower left corner at (K_min, Con_min)
% square size = d(1) * d(2)

nvals = 30;
nvar = 4;
x = lhsdesign(nvals, nvar);

% K12 K21 Con1 Con2
K12r = [2 40]; 
K21r = [2 40];
Con1r = [2 40];
Con2r = [2 40];

K12_all = (K12r(2)-K12r(1))*x(:,1) + K12r(1);
K21_all = (K12r(2)-K12r(1))*x(:,2) + K21r(1);
Con1_all = (Con1r(2)-Con1r(1))*x(:,3) + Con1r(1);
Con2_all = (Con2r(2)-Con2r(1))*x(:,4) + Con2r(1);

% plot projections
figure();
hold on
%plot(K12_all, K21_all, 'x');
plot(Con1_all, Con2_all, 'x');
grid on;
dK12 = (K12r(2)-K12r(1))/nvals;
dK21 = (K21r(2)-K21r(1))/nvals;
dCon1 = (Con1r(2)-Con1r(1))/nvals;
dCon2 = (Con2r(2)-Con2r(1))/nvals;
%set(gca, 'XTick', K12r(1):dK12:K12r(2), 'YTick',  K21r(1):dK21:K21r(2));
set(gca, 'XTick', Con1r(1):dCon1:Con1r(2), 'YTick',  Con2r(1):dCon2:Con2r(2));
xlim([0 40]);
ylim([0 40]);
close all
%% Export values
qsave = 1;
if qsave
    path = fullfile('H:\My Documents\Multicellular automaton\data\two_signals\LHS_sample');
    fname_str = sprintf('LHS_values_K12_K21_Con1_Con2_nvals%d',...
        nvals);
    i = 1;
    fname = fullfile(path, ...
        strcat(fname_str,'-v',int2str(i),'.mat'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(path, ...
            strcat(fname_str,'-v',int2str(i),'.mat'));
    end
    save(fname);
end