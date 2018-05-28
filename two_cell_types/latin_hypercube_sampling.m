%% Implements Latin hypercube sampling for a region of K, Con we wish to examine

% take square in K, Con space with lower left corner at (K_min, Con_min)
% square size = d * d

K_min = 2;
Con_min = 5;
d = 15; 
nvals = d;
nvar = 2;
x = lhsdesign(nvals,nvar);
Kvals = d*x(:,1) + K_min;
Convals = d*x(:,2) + Con_min;

figure();
hold on
%plot(x(:,1), x(:,2), 'x');
plot(Kvals, Convals, 'x');
grid on;
%set(gca, 'XTick', 0:1/nvals:1, 'YTick', 0:1/nvals:1);
set(gca, 'XTick', K_min:1:K_min+d, 'YTick', Con_min:1:Con_min+d);
xlim([0 30]);
ylim([0 30]);

%% Export values
qsave = 1;
if qsave
path = fullfile(pwd, 'data', 'LHS');
fname_str = sprintf('LHS_values_K%d_to_%d_Con%d_to_%d',...
    K_min, K_min+d, Con_min, Con_min+d);
i = 1;
fname = fullfile(path, ...
    strcat(fname_str,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(path, ...
        strcat(fname_str,'-v',int2str(i),'.mat'));
end
save(fname, 'Kvals', 'Convals');
end