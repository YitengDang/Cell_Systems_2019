%% evaluates how well generate_I.m does its job by obtaining statistics on achieved I and running times
clear variables
close all

gridsize = 15;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;

p1 = floor(0.1*N);
dp = round(0.02*N);
p_all = (p1:dp:N-p1)/N; % initial p
%p = (30:35)/N;
dI = 0.01; % accuracy of generated I
I_all = -0.2:0.05:0.8; % initial I
ntrials = 100; % trials per initial p,I point
times = zeros(numel(I_all), numel(p_all), ntrials);
I_out = zeros(numel(I_all), numel(p_all), ntrials);

%I_min = -0.1;
%I_max = I_min+dI;
%}
%% Load saved data
%{
% round(p(end)*N))
%fname_str = strrep(sprintf('generate_I_N%d_a0_%.1f_I_all_Non_%d_%d_fracI_accept',...
%       N, a0, round(p(1)*N), 109), '.', 'p');
%fname_str = 'generate_I_N121_a0_1p0_I_all_Non_12_109_fracI_accept';
fname_str = 'generate_I_N121_a0_1p0_I_all_Non_12_109_fracI_accept_approx_new';
fname = fullfile(pwd, 'figures', 'generate_I', 'data', strcat(fname_str, '.mat'));
load(fname);

% Process wrong data to fix problem (temp)
data = [nI_accept_new(1:36, :)/10000; nI_accept_new(37:end, :)/100];
%}

%%
%{
I_all_new = 0.52:0.02:1.0;
nI_accept_all = zeros(numel(I_all)+numel(I_all_new), numel(p)); 
nI_accept_all(1:numel(I_all), :) = fracI_accept; % ntrials*fracI_accept;
idx_start = numel(I_all);
%}
%% Run
%
% Initialize parameters
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
cell_type = zeros(N, 1);
nI_accept = zeros(numel(I_all), numel(p_all));
for k=1:numel(I_all)  
    I_min = I_all(k);
    I_max = I_min+dI;
    for i=1:numel(p_all)
        iniON = round(p_all(i)*N);
        fprintf('Imin = %.2f, iniON = %d \n', I_min, iniON);
        test_aux = zeros(ntrials, 1);
        parfor j=1:ntrials
            %fprintf('Imin = %.2f, iniON = %d, trials = %d \n', I_min, iniON, j);
            cells_in = zeros(N, 1);
            cells_in(randperm(N, iniON)) = 1;
            [~, ~, ~, test_aux(j)] = generate_I_new(cells_in, I_min, I_max, dist, a0);
            %h2=figure();
            %update_cell_figure(h2, pos, 1, cells_out, cell_type, 0);
        end
        nI_accept(k,i) = sum(test_aux);
    end
end
fracI_accept = nI_accept/ntrials;
%}

% Fraction of I_out in desired range
%{
fracI_accept = zeros(numel(I_all), numel(p));
for k=1:numel(I_all)
    I_min = I_all(k); I_max = I_min+dI;
    fracI_accept(k, :) = sum((I_min < I_out(k, :, :) )&(I_out(k, :, :) < I_max), 3);
end
%}
%%
qsave = 1;
if qsave
   fname_str = strrep(sprintf('generate_I_N%d_a0_%.1f_I_all_Non_%d_%d_fracI_accept_approx',...
       N, a0, round(p_all(1)*N), round(p_all(end)*N)), '.', 'p');
   fname = fullfile(pwd, 'figures', 'generate_I', 'data', fname_str);
   save(fname);
   %save(fname, 'nI_accept', 'fracI_accept', 'N', 'a0', 'I_all', 'p', 'ntrials', 'dI', 'I_out', 'times');
end
%}


%% Compare with Imax
%{
Imax_vals = zeros(numel(p), 1);
for i=1:numel(p)
    thisp = p(i);
    [~, Imax_vals(i)] = max_moranI(a0, Rcell, N, thisp);
end
figure();
plot(p, Imax_vals);
%}
%% Plot accepted fraction
h=figure();
%I_all = [I_all I_all_new];
imagesc(p_all, I_all, fracI_accept);
%imagesc(p, I_all, nI_accept_all);
%imagesc(p_all, I_all, data);
set(gca, 'Ydir', 'normal', 'FontSize', 24);
c = colorbar;
xlabel('p');
ylabel('I');
ylabel(c, 'fraction accepted');

qsave = 1;
if qsave
   fname_str = strrep(sprintf('generate_I_times_N%d_a0_%.1f_I_%.3f_%.3f_fracI_accept_approx',...
       N, a0, I_all(1), I_all(end)), '.', 'p');
   fname = fullfile(pwd, 'figures', 'generate_I', fname_str);
   save_figure_pdf(h, 10, 8, fname);
   save_figure_eps(h, 10, 8, fname);
   %fname2 = fullfile(pwd, 'temp', fname_str);
   %save_figure_pdf(h, 10, 8, fname2);
end

%
%% Mean generated I
%{
set(0, 'defaulttextinterpreter', 'latex');
h1=figure(1);
hold on
plot(p, mean(I_out, 2));
plot([0 1], [I_min I_min], 'r-');
plot([0 1], [I_max I_max], 'r-');
xlabel('$$p$$');
ylabel('$$\langle I_{out} \rangle$$');
set(gca, 'FontSize', 24);
%ylim([0 I_max+0.05]);

qsave = 0;
if qsave
   fname_str = strrep(sprintf('generate_I_times_N%d_I_%.3f_%.3f_meanI', N, I_min, I_max),...
       '.', 'p');
   fname = fullfile(pwd, 'figures', 'generate_I', fname_str);
   save_figure_pdf(h1, 10, 8, fname);
   save_figure_eps(h1, 10, 8, fname);
end
%% Density generated I
h2=figure(2);
%dI = 0.01;
Iedges = floor(min(min(I_out))/dI)*dI:dI:ceil(max(max(I_out))/dI)*dI;
count = zeros(numel(Iedges)-1, numel(p));
for idx=1:numel(p)
    count(:, idx) = histcounts(I_out(idx, :), Iedges);
end
Ic = (Iedges(1:end-1)+Iedges(2:end))/2;
imagesc(p, Ic, count/ntrials)
c = colorbar;
xlabel('$$p$$');
ylabel('$$I_{out}$$');
ylabel(c, '$$P(I_{out}|p)$$', 'interpreter', 'latex');
set(gca, 'YDir', 'normal', 'FontSize', 24);

qsave = 1;
if qsave
   fname_str = strrep(sprintf('generate_I_times_N%d_I_%.3f_%.3f_Iprobp', N, I_min, I_max),...
       '.', 'p');
   fname = fullfile(pwd, 'figures', 'generate_I', fname_str);
   save_figure_pdf(h2, 10, 8, fname);
   save_figure_eps(h2, 10, 8, fname);
end
%% Mean equilibration times
h3=figure(3);
plot(p, mean(times, 2));
xlabel('$$p$$');
ylabel('$$\langle t \rangle$$');
set(gca, 'FontSize', 24);

qsave = 1;
if qsave
   fname_str = strrep(sprintf('generate_I_times_N%d_I_%.3f_%.3f_mean_times', N, I_min, I_max),...
       '.', 'p');
   fname = fullfile(pwd, 'figures', 'generate_I', fname_str);
   save_figure_pdf(h3, 10, 8, fname);
   save_figure_eps(h3, 10, 8, fname);
end
%}