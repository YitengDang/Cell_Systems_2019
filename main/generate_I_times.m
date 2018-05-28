%% evaluates how well generate_I.m does its job by obtaining statistics on achieved I and running times
clear variables
close all

gridsize = 11;
N = gridsize^2;

%p0 = 0.2;
p1 = floor(0.1*N);
p = (p1:N-p1)/N;
ntrials = 100;
times = zeros(numel(p), ntrials);
I_out = zeros(numel(p), ntrials);

I_min = -0.1;
dI = 0.01;
I_max = I_min+dI;

% Initialize parameters
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
cell_type = zeros(N, 1);
    
for i=1:numel(p)
    iniON = round(p(i)*N);
    for j=1:ntrials
        fprintf('iniON = %d, trials = %d \n', iniON, j);
        cells_in = zeros(N, 1);
        cells_in(randperm(N, iniON)) = 1;
        [~, I_out(i,j), times(i,j)] = generate_I_new(cells_in, I_min, I_max, dist);
        %h2=figure();
        %update_cell_figure(h2, pos, 1, cells_out, cell_type, 0);
    end
end

%% Fraction of I_out in desired range
fracI_accept = sum((I_min < I_out)&(I_out < I_max), 2);
%%
qsave = 1;
if qsave
   fname_str = strrep(sprintf('generate_I_times_stats_N%d_I_%.3f_%.3f_Non_%d_%d',...
       N, I_min, I_max, round(p(1)*N), round(p(end)*N)), '.', 'p');
   fname = fullfile(pwd, 'figures', 'generate_I', 'data', fname_str);
   save(fname);
end

%{
%% Mean generated I
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

qsave = 0;
if qsave
   fname_str = strrep(sprintf('generate_I_times_N%d_I_%.3f_%.3f_Iprobp', N, I_min, I_max),...
       '.', 'p');
   fname = fullfile(pwd, 'figures', 'generate_I', fname_str);
   save_figure_pdf(h2, 10, 8, fname);
   save_figure_eps(h2, 10, 8, fname);
end
%%
h3=figure(3);
plot(p, mean(times, 2));
xlabel('$$p$$');
ylabel('$$\langle t \rangle$$');
set(gca, 'FontSize', 24);

qsave = 0;
if qsave
   fname_str = strrep(sprintf('generate_I_times_N%d_I_%.3f_%.3f_mean_times', N, I_min, I_max),...
       '.', 'p');
   fname = fullfile(pwd, 'figures', 'generate_I', fname_str);
   save_figure_pdf(h3, 10, 8, fname);
   save_figure_eps(h3, 10, 8, fname);
end
%}