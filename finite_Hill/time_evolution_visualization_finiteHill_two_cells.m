% Time evolution of a system without noise and with visualization
% Finite Hill coefficient
% with only 2 cells
close all
clear all
warning off

% Lattice parameters
N=2;
a0 = 5.6;
Rcell = 0.2*a0;
% circuit parameters
K = 8;
Con = 16;
hill = 2; % Hill coefficient
% Initial configuration
% If fixedp, assume initial states ~ Normal(p,sigma)
initialID = 'uniform'; %'uniform' 'singleOFF' 'fixedp' 'fixedp_unif'
p = 0.27;
sigma = 0.1;

% simulation parameters
qsave=0;
tmax = 1000;

% use hexagonal lattice
gridsize = 2;
dist = [0 1; 1 0];
pos = [0 0; 1 0];

dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% generate cell_type (0 case type 1, 1 case type 2)
cell_type = zeros(N,1); % all the same here

%% Plot phase diagram
set(0, 'defaulttextinterpreter', 'latex');
%
% Bistability example: K=8, Con=16, a0=1.5, hill=2
set(0,'defaulttextinterpreter', 'latex');
h=figure();
hold on
S = linspace(1,Con,100);
fS = 1 + (Con-1)*S.^hill./(K^hill+S.^hill);
plot((S-1)/(Con-1), (S-1)/(Con-1), 'b-');
plot((S-1)/(Con-1), (fS-1)/(Con-1), 'r-');
xlabel('Sensed concentration $$Y_i$$', 'FontSize', 20);
ylabel('Secreted concentration $$f(Y_i)$$', 'FontSize', 20);
legend({'Y_i', 'f(Y_i)'}, 'FontSize', 20, 'Location', 'nw');
%Save fig
if qsave
    fname = strrep(sprintf('single_cell_dynamics_Con%d_K%d_hill%.2f', Con, K, hill), '.','p');
    out_file = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'figures', fname);
    save_figure_pdf(h, 10,8, out_file)
end

% Save figure
if qsave
    fname = strrep(sprintf('stability_diagram_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s',...
        N, a0, Con, K, hill, initialID), '.','p');
    out_file = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'fig_bistable_sim', fname);
    save_figure_pdf(h, 10, 8, out_file)
end
%}
%% initialize cells
cells_hist = {};
if strcmp(initialID, 'uniform')
    cells_hist{1} = rand(N, 1); %uniformly distributed
elseif strcmp(initialID, 'singleOFF')
    cells_hist{1} = ones(N, 1); 
    x = randi(N);
    cells_hist{1}([x x+1 x+gridsize x+gridsize+1]) = 0;
elseif strcmp(initialID, 'fixedp')
    cells_hist{1} = p + sigma*randn(N, 1);
    cells_hist{1}(cells_hist{1} < 0) = 0;
    cells_hist{1}(cells_hist{1} > 1) = 1;
elseif strcmp(initialID, 'fixedp_unif')
    cells_hist{1} = p*ones(N, 1); 
end
%cells = [0.62*ones(floor(N/2),1); zeros(ceil(N/2),1)]; % extreme case: separation of ON/OFF cells

% Load previous simulation
%fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
%    'N225_a01p5_Con18p00_K14p00_hill2p00_fixedp_pini0p27_sigma0p10-v1.mat');
%load(fname, 'cells_hist');
%% Run simulation
cells = cells_hist{1};
% video for saving
frames = struct('cdata',[],'colormap',[]);

% initialize figure
hin = figure(1);
t = 0;
H = [];
xmean = [];
xstd = [];
I = [];
%cells_hist = {};

update_cell_figure_continuum(hin, pos, a0, cells, cell_type, t);
frames(1) = getframe(gcf);
%cells_hist{end+1} = cells;
xmean(end+1) = mean(cells);
xstd(end+1) = std(cells);
I(end+1) = moranI(cells, a0*dist);
[cells_out, changed, H(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill);

while changed && t<tmax
    t = t+1;
    %k = waitforbuttonpress;
    %pause(0.1);
    update_cell_figure_continuum(hin, pos, a0, cells_out, cell_type, t);
    frames(end+1) = getframe(gcf);
    %cells_hist{end+1} = cells_out;
    xmean(end+1) = mean(cells_out);
    xstd(end+1) = std(cells_out);
    I(end+1) = moranI(cells, a0*dist);
    cells = cells_out;
    [cells_out, changed, H(end+1)] = update_cells_continuum(cells, dist, Con, K, a0, Rcell, hill);
end
%% Plot muX(t)=<Xi>(t) and sigmaX(t)
h=figure(3);
hold on
%errorbar(0:t, xmean, xstd, 'bo-');
plot(0:t, xmean, 'bo-');
plot(0:t, xstd, 'ro-');

set(gca, 'FontSize', 20);
xlabel('time');
ylabel('$$\langle X_i \rangle$$'); % correct label
%xlim([0 40]);
%legend({'\langle X_i \rangle', '\sigma_X'});
%legend({'\mu_X (Standard)', '\sigma_X (Standard)', '\mu_X (Mean field)', '\sigma_X (Mean field)'});

% Calculate trajectory with all cells Xm
Xm_ini = mean(cells_hist{1}); % initial Xm
tsteps = t+1;
Xm = zeros(1, tsteps+1);
Xm(1) = Xm_ini;
for i=1:tsteps
    Y = (1+fN)*(Xm(i)*(Con-1) + 1);
    Xm(i+1) = Y^hill/(K^hill + Y^hill);
end

h;
hold on
plot(0:tsteps-1, Xm(1:end-1), 'g-', 'LineWidth', 1.5);
%legend({'simulation', 'calculation'});
legend({'\mu_X', '\sigma_X',...
    '\mu_X (Unif)'}, 'Location', 'e');

%Save fig
if qsave
    fname_str = strrep(sprintf('trajectory_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
        N, a0, Con, K, hill, initialID, p, sigma), '.','p');
    i = 1;
    fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'fig_bistable_sim',...
        strcat(fname_str,'-v',int2str(i),'.pdf'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'fig_bistable_sim',...
            strcat(fname_str,'-v',int2str(i),'.pdf'));
    end
    save_figure_pdf(h, 10, 8, fname)
end
%% Plot hamiltonian
%{
% calculate hamiltonian for uniform lattice
%hm = -Xm.*((1+fN)*(Xm*(Con-1)+1) - K);
hm = (Xm(2:end) - Xm(1:end-1)).^2;

h=figure(4);
hold on
tmax = 10;
plot(0:tmax, H(1:tmax+1)/N, 'b','LineWidth', 1.5); %h = H/N
plot(0:tmax, HMF(1:tmax+1)/N, 'r','LineWidth', 1.5); %h = H/N
plot(0:tmax, hm(1:tmax+1), 'g','LineWidth', 1.5); %h = H/N
set(gca, 'FontSize', 20);
xlabel('time');
ylabel('$$h(t)$$');
legend({'Standard', 'Mean field', 'Uniform lattice'});

%Save fig
if qsave
    fname_str = strrep(sprintf('hamiltonian_new_N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
        N, a0, Con, K, hill, initialID, p, sigma), '.','p');
    i = 1;
    fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'fig_bistable_sim',...
        strcat(fname_str,'-v',int2str(i),'.pdf'));
    while exist(fname, 'file') == 2
        i=i+1;
        fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'fig_bistable_sim',...
            strcat(fname_str,'-v',int2str(i),'.pdf'));
    end
    save_figure_pdf(h, 10, 8, fname)
end
%% Plot (p,I) trajectory
figure(5);
hold on

% plot hamiltonian map
E = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
                -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);
piv = (0:N)/N;
Iv = -0.05:0.05:1;
[p_i,Imesh] = meshgrid(piv, Iv);
contourf(piv, Iv, E(p_i, Imesh),'LineStyle', 'none');
colormap('summer')
c = colorbar;

% plot trajectories
p1=plot(xmean, I, 'b');
%plot endpoints
plot(xmean(end), I(end),'bx');

% add legend
%legend([p1 p2], {'normal', 'mean field'});
%}
%% Save whole trajectory if interesting
%{
fname_str = strrep(sprintf('N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
    N, a0, Con, K, hill, initialID, p, sigma), '.','p');
i = 1;
fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
    strcat(fname_str,'-v',int2str(i),'.mat'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
        strcat(fname_str,'-v',int2str(i),'.mat'));
end
save(fname)
%}
%%
% Check movies
%figure(6);
%movie(framesMF, 1, 2)

%%
%{
% Export movies
% normal
fname_str = strrep(sprintf('N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f',...
    N, a0, Con, K, hill, initialID, p, sigma), '.','p');
i = 1;
fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
    strcat(fname_str,'-v',int2str(i),'.avi'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
        strcat(fname_str,'-v',int2str(i),'.avi'));
end

myVideo = VideoWriter(fname, 'Uncompressed AVI');
myVideo.FrameRate = 3;  % Default 30
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);
%%
% mean field 
fname_str = strrep(sprintf('N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f_MF',...
    N, a0, Con, K, hill, initialID, p, sigma), '.','p');
i = 1;
fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
    strcat(fname_str,'-v',int2str(i),'.avi'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(pwd, 'latex', '17-07-14_finite_Hill', 'interesting_trajectories',...
        strcat(fname_str,'-v',int2str(i),'.avi'));
end

myVideo = VideoWriter(fname, 'Uncompressed AVI');
myVideo.FrameRate = 15;  % Default 30
open(myVideo);
writeVideo(myVideo, framesMF);
close(myVideo);
%}