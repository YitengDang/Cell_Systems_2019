% This script calculates the map of p_in p_eq using exact simulation.
close all
clear all
%warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
Con = 5;
K = 11;
%7.3 8.5 9.9 11.0 12.5 

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);

% Calculate the signaling strength
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(Rcell*sum(exp(Rcell-r)./r)); % calculate signaling strength
kon = N*(K-fN-Con)/(Con-1)/fN;
koff = N*(K-fN-1)/(Con-1)/fN;

%% calculate the map
[count, t_av] = count_eq_parallel(dist, Con, K, a0, Rcell);
prob = transpose(count./repmat(sum(count,2),1,N+1));

% save mat file
%save(fullfile(pwd,'figures','pin_pout',strcat(fname,'.mat')))
%% Calculate weighted hamiltonian
hamiltonian = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
pv = (0:N)./N;
Iv = -0.05:0.001:0.05;
probIp = ones(1, numel(Iv))/numel(Iv);
%probIp = normpdf(Iv, -0.01, 0.1)/sum(normpdf(Iv, -0.01, 0.03));
%figure();
%plot(Iv, probIp);

% average Hamiltonian over the Iv range
[pmesh, Imesh] = meshgrid(pv, Iv);
hmesh = hamiltonian(pmesh, Imesh);
h_weighted = probIp*hmesh;

figure();
hold on
plot(pv, hamiltonian(pv,0), 'b');
plot(pv, h_weighted, 'g');

% find max and plot
idx0 = find(max(hamiltonian(pv, 0)) == hamiltonian(pv,0)); % find index corresponding to p*
idx = find(max(h_weighted) == h_weighted); % find index corresponding to p*
plot([ pv(idx0)  pv(idx0)], [-10 1], 'g--', 'LineWidth', 2); %
plot([ pv(idx)  pv(idx)], [-10 1], 'b--', 'LineWidth', 2);
%%
% plot the map
h1 = figure(1);
hold on

% Plot line p = 1/2 - 4B/fN
B = (Con+1)/2*(1+fN) - K;
pline = 1/2 - B/(4*fN);
plot([pline pline], [0 1], 'r--', 'LineWidth', 2);

% Plot line p* defined by Hmax = H(p*)
%idx0 = find(max(hamiltonian(p, 0)) == hamiltonian(p,0)); % find index corresponding to p*
idx = find(max(h_weighted) == h_weighted); % find index corresponding to p*
ps = pv(idx); %p*
plot([ps ps], [0 1], 'g--', 'LineWidth', 2);

% Plot only H(p,0)
yyaxis right
plot(pv, hamiltonian(pv,0), 'b-', 'LineWidth', 1.5);
ylabel('h(p, 0)');

% legend
%legend({'p=1/2-B/(4fN)', 'p*', 'h(p,0)'});

% plot pin-pout map
yyaxis left
im_fig = imagesc(pv,pv,prob);
% set title and font
title(sprintf('N = %d, K = %.1f, C_{ON} = %.1f, a_0 = %.1f, R = %.1f', ...
    N, K, Con, a0, Rcell),'FontSize', 18)
set(gca,'Ydir','normal','FontSize', 20)
% set invisible parts where count is zero
set(im_fig, 'AlphaData', count' > 0);
% set colorbar and labels
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{in}', 'FontSize', 24)
ylabel('p_{out}', 'FontSize', 24)
xlim([pv(1) pv(end)]);
ylim([0 1]);

% Organize and save
save_fig = 0; % save figure? 0:no, 1: yes
if save_fig > 0
    fname = strrep(sprintf('pin_pout_N%d_Con_%.1f_K_%.1f_a0_%.1f_I0_%.1f',...
        N, Con, K, a0, I0), '.', 'p');
    set(h1,'Units','Inches');
    set(h1, 'Position', [0 0 10 6 ]);
    pos = get(h1, 'Position');
    set(h1, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
    out_file = fullfile(pwd, 'rebuttal', 'hmax_pin_pout', strcat(fname,'_map'));
    print(h1, out_file,'-dpdf','-r0')
end
%% Plot Hamiltonian map H(p,I)
h2 = figure(2);
hold on
Iv = -0.05:0.05:1;
[pmesh,Imesh] = meshgrid(pv, Iv);
contourf(pv, Iv, hamiltonian(pmesh, Imesh),'LineStyle', 'none')
colormap('summer')
c = colorbar;

% Plot Hmax at (p,I) = (p*, 0)
plot(ps, 0, 'rx');
%%
% Plot the average number of steps it takes to reach equilibrium
h3 = figure(3);
plot(pv, t_av, 'r-o')
set(gca,'FontSize', 20)
title(sprintf('N = %d, K = %.1f, S_{ON} = %.1f, a0 = %.1f, R = %.1f', ...
    N, K, Con, a0, Rcell),'FontSize', 18)
xlabel('k_{in}', 'FontSize', 24)
ylabel('Average # steps for eq.', 'FontSize', 24)
%%
h4 = figure(4);
% calculate the analytical formula and plot
[~,omegak] = entropy_eq_sphere(dist_vec, Con, K, a0, Rcell);
plot((0:N)/N , log(omegak), 'LineWidth', 1.5)
set(gca,'FontSize', 20)
xlabel('p');
ylabel('$$\Omega_k$$', 'Interpreter', 'latex');
