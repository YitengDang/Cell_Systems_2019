% Creates pin-pout map from (p,I) simulation 
% Previously called pin_pout_analytical, but uses MC, not analytical
% calculation
% YD: looked like the newest version for pin_pout_analytical;

clear variables
close all
qsave = 0;

% Get phenotype regions
gridsize = 15;
a0 = 1.5;
Rcell = 0.2*a0;
K = 6;
Con = 15;
alpha = 0; 

% Initial state
Ii = 0;

N = gridsize^2;
[dist, pos] = init_dist_hex(gridsize, gridsize);

dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence

% Get the signalling length
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength sphere
M = sinh(Rcell)*exp(Rcell-a0*dist)./(a0*dist);
M(1:N+1:N^2) = 0;
aux = M*M';
gN = aux(1,1);
xiN = sum(aux(1,2:end));

% Calculate the approximate dynamics of p and H
hfunc = @(p, I) -0.5*(Con-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
    -(2*p-1).*(0.5*(Con+1)*(1+fN) - K);

trials = 1000;
cntmap = zeros(N+1);
t_av = zeros(N+1, 1);

if Ii == 0
    fname = sprintf('pin_poutN%d_Con_%d_K_%d_gz_%d_a0_%d_analytic.mat', ...
            N, Con, K, gridsize, 10*a0);
else
    fname = sprintf('pin_poutN%d_Con_%d_K_%d_gz_%d_a0_%d_I%d_analytic.mat', ...
            N, Con, K, gridsize, 10*a0, 10*Ii);
end

for n = 0:N
    disp(n)
    for i = 1:trials
        p = n/N;
        I = Ii;
        theta = fN*(4*I*p*(1-p)+(2*p-1)^2);
        cont = true;
        t = 0;
        while cont
            [theta_new, p_new, pe] = update_stochastic(theta, p, N, Con, K, fN, gN, alpha);
            if rand < pe
                cont = false;
                %if pe > 0.1
                    t_av(n+1) = t_av(n+1) + t/trials;
                    nout = round(N*p);
                    cntmap(n+1, nout + 1) = cntmap(n+1, nout + 1) + 1/trials;
                %end
            else
                t = t + 1;
                p = p_new;
                theta = theta_new;
                eps = 1e-6;
                if p < eps || p > 1-eps
                    I = 0;
                else
                    I = (theta - (2*p-1)^2*fN)/4/p/(1-p)/fN;
                    %fprintf('I: %.4f, p: %.4f, theta: %.4f \n',I,p,theta);
                    %I(i) = I(i-1) + dI;
                end
            end
        end
    end
end

if qsave
    save(fullfile(pwd, 'figures', 'pin_pout', fname));
end
%% Plot and save figures
set(0, 'defaulttextinterpreter', 'latex');

if Ii == 0
    fname = strrep(sprintf('pin_pout_MC_v1_EO_N%d_a0_%.2f_Con_%d_K_%d_analytic', ...
            N, a0, Con, K), '.', 'p');
else
    fname = strrep(sprintf('pin_pout_MC_v1_EO_N%d_a0_%.2f_Con_%d_K_%d_I_%.2f_analytic', ...
            N, a0, Con, K, Ii), '.', 'p');
end

h1=figure(1);
p = (0:N)./N;
im_fig = imagesc(p,p,cntmap');
set(gca,'Ydir','normal','FontSize', 24)
set(im_fig, 'AlphaData', cntmap' > 0);
c = colorbar;
c.Label.String = 'Probability';
xlabel('$$p_{in}$$', 'FontSize', 24)
ylabel('$$p_{eq}$$', 'FontSize', 24)

qsave = 1;
if qsave
    out_file = fullfile(pwd, 'figures', 'pin_pout', strcat(fname, '_map'));
    save_figure_pdf(h1, 10, 8, out_file);
end

h2=figure(2);
plot(p, t_av)
xlabel('$$p_{in}$$')
ylabel('$$\langle t_{eq} \rangle$$')
set(gca,'FontSize', 24)

if qsave
    out_file = fullfile(pwd, 'figures', 'pin_pout', strcat(fname, '_tav'));
    save_figure_pdf(h2, 10, 8, out_file);
end