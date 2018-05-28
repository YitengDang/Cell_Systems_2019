% Creates pin-pout map from (p,I) Monte Carlo simulation with noise
clear variables
close all

data_path = 'C:\Users\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';

% Get phenotype regions
gridsize = 15;
a0 = 0.5;
Rcell = 0.2*a0;
Con = 14;
K = 15;
alpha = 1;

% Initial state
Ii = 0;

% These variables are to stop the running
err_thr_n = 1; % variation allowed for number of ON cells
err_thr_I = 0.02; % variation allowed for I
len_test = 10; % window (in steps) of the moving std

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
I_av = zeros(N+1, 3);
if Ii == 0
    fname = sprintf('pin_poutN%d_Con_%d_K_%d_gz_%d_a0_%d_noise_%d_montecarlo.mat', ...
            N, Con, K, gridsize, 10*a0, round(10*alpha));
else
    fname = sprintf('pin_poutN%d_Con_%d_K_%d_gz_%d_a0_%d_I%d_noise_%d_montecarlo.mat', ...
            N, Con, K, gridsize, 10*a0, 10*Ii, round(10*alpha));
end
    
for n = 0:N
    disp(n)
    Iaux = zeros(trials,1);
    for i = 1:trials
        p = n/N;
        I = Ii;
        I_test = zeros(len_test, 1);
        p_test = zeros(len_test, 1);
        theta = fN*(4*I*p*(1-p)+(2*p-1)^2);
        cont = true;
        t = 0;
        while cont
            if t < len_test
                t = t+1;
                ind_aux = len_test - t+1; % populate in reverse order (older last)
                [theta, p, ~] = ...
                    update_montecarlo(theta, p, N, Con, K, fN, gN, alpha);
                p_test(ind_aux) = p;
                I = calc_I(p, theta, fN);
                I_test(ind_aux) = I;
            else
                if std(p_test) < err_thr_n && std(I_test) < err_thr_I
                    cont = false;
                    t_av(n+1) = t_av(n+1) + t/trials;
                    Iaux(i) = mean(I_test);
                    nout = round(N*mean(p_test));
                    cntmap(n+1, nout + 1) = cntmap(n+1, nout + 1) + 1/trials;
                else
                    I_test = circshift(I_test,1); % rotate the list to replace the last 
                    p_test = circshift(p_test, 1);
                    [theta, p, ~] = ...
                        update_montecarlo(theta, p, N, Con, K, fN, gN, alpha);
                    t = t + 1;
                    I = calc_I(p, theta, fN);
                    p_test(1) = p;
                    I_test(1) = I;
                end
            end
        end
    end
    I_av(n+1,:) = [mean(Iaux) min(Iaux) max(Iaux)]; 
end

save(fullfile(data_path,'pin_pout_noise',fname));

figure(1)
p = (0:N)./N;
im_fig = imagesc(p,p,cntmap');
set(gca,'Ydir','normal','FontSize', 20)
set(im_fig, 'AlphaData', cntmap' > 0);
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{ini}', 'FontSize', 24)
ylabel('p_{eq}', 'FontSize', 24)

figure(2)
plot(p, t_av)
xlabel('p_{in}')
ylabel('Time')

