% Creates pin-pout map from (p,I) simulation (analytical = ?)
clear variables
close all

data_path = 'D:\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';

% Get phenotype regions
gridsize = 15;
a0 = 1.5;
Rcell = 0.2*a0;
Con = 11.4;
K = 5.2;
alpha = 0.26;

% These variables are to stop the running
err_thr_n = 1; % variation allowed for number of ON cells
err_thr_I = 0.01; % variation allowed for I
len_test = 10;

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
phi = @(x) exp(-0.5*x^2)/sqrt(2*pi);

n_smpl = 1000;
t_av = zeros(N+1, 1);
count = zeros(N+1,N+1);
I_av = zeros(N+1,3);

if Ii == 0
    fname = sprintf('pin_poutN%d_Con_%d_K_%d_gz_%d_a0_%d_analytic_noise%d.mat', ...
            N, round(Con), round(K), gridsize, round(10*a0), round(10*alpha));
else
    fname = sprintf('pin_poutN%d_Con_%d_K_%d_gz_%d_a0_%d_I%d_analytic_noise%d.mat', ...
            N, round(Con), round(K), gridsize, round(10*a0), round(100*Ii), round(10*alpha));
end
    
for n = 0:N
    disp(n)
    t_aux = zeros(n_smpl,1);
    k_aux = zeros(n_smpl,1);
    I_aux = zeros(n_smpl,1);
    parfor i = 1:n_smpl
        p = n/N;
        I = Ii;
        theta = fN*(4*I*p*(1-p)+(2*p-1)^2);
        cont = true;
        I = zeros(len_test,1);
        fracON = zeros(len_test,1);
        t = 0;
        while cont
            if t < len_test
                t = t+1;
                ind_aux = len_test - t+1; % populate in reverse order (older last)
                [theta, p, pe] = update_stochastic_cov(theta, p, N, ...
                    a0, Rcell, Con, K, fN, gN, alpha);
                if p == 0 || p == 1
                    I(ind_aux) = 0;
                else
                    I(ind_aux) = (theta - (2*p-1)^2*fN)/4/p/(1-p)/fN;
                end
                fracON(ind_aux) = N*p;
            else
                if std(fracON) < err_thr_n && std(I) < err_thr_I
                    cont = false;
                    t_aux(i) = t;
                    k_aux(i) = round(mean(fracON));
                    I_aux(i) = mean(I);
                else
                    t = t+1;
                    I = circshift(I,1); % rotate the list to replace the last 
                    fracON = circshift(fracON, 1);
                    [theta, p, pe] = update_stochastic_cov(theta, p, N, ...
                        a0, Rcell, Con, K, fN, gN, alpha);
                    if p == 0 || p == 1
                        I(1) = 0;
                    else
                        I(1) = (theta - (2*p-1)^2*fN)/4/p/(1-p)/fN;
                    end
                    fracON(1) = round(N*p);
                end
            end
        end
    end
    t_av(n+1) = sum(sum(t_aux))/n_smpl;
    for i = 0:N
        count(n+1,i+1) = count(n+1,i+1) + sum(sum(k_aux == i));
    end
    I_av(n+1,1) = sum(sum(I_aux))/n_smpl;
    I_av(n+1,2) = min(I_aux);
    I_av(n+1,3) = max(I_aux);
end

save(fullfile(data_path,'pin_pout',fname));

figure(1)
p = (0:N)./N;
im_fig = imagesc(p,p,count'/n_smpl);
set(gca,'Ydir','normal','FontSize', 20)
set(im_fig, 'AlphaData', count' > 0);
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{ini}', 'FontSize', 24)
ylabel('p_{eq}', 'FontSize', 24)

figure(2)
plot(p, t_av)
xlabel('p_{in}')
ylabel('Time')

