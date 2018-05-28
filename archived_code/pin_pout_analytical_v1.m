clear variables
close all

data_path = 'D:\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';

% Get phenotype regions
gridsize = 15;
a0 = 0.5;
Rcell = 0.2*a0;
Con = 8;
K = 16;
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
phi = @(x) exp(-0.5*x^2)/sqrt(2*pi);

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
            muon = fN*(Con*p + 1 - p + (Con-1)*(1-p)*I);
            muoff = fN*(Con*p + 1 - p - (Con-1)*p*I);
            kappa = (Con-1)^2*(gN)*p*(1-p)/ ...
                (p*muon + (1-p)*muoff);
            sigmaon = sqrt(max(kappa*muon,0)+alpha^2);
            sigmaoff = sqrt(max(kappa*muoff,0)+alpha^2);
            zoff = (K - 1 - muoff)/sigmaoff;
            zon = (K - Con - muon)/sigmaon;
            Poffoff = normcdf(zoff);
            Ponon = 1-normcdf(zon);
            ponoff = p*(1-Ponon);
            poffon = (1-p)*(1-Poffoff);
            dp = poffon-ponoff;
            vardp = (poffon*Poffoff + ponoff*Ponon)/N;
            dtheta = 8/(Con-1)*(poffon*muoff - ponoff*muon + ...
                sigmaoff*(1-p)*phi(zoff)+ sigmaon*p*phi(zon)) ...
                - 4*fN*(Con+1)/(Con-1)*dp + 4*fN*(dp^2 + vardp);
            
            vartheta = 64/(Con-1)^2*((1-p)^2*sigmaoff^2* ...
                ((1-Poffoff)^2+zoff*phi(zoff)*(1-Poffoff) - (phi(zoff))^2) ...
                + p^2*sigmaon^2* ...
                ((1-Ponon)^2-zon*phi(zon)*(1-Ponon) - (phi(zon))^2)) + ...
                16*fN^2*(4*dp^2*vardp + 2*vardp^2 + (Con+1)^2/(Con-1)^2*vardp);
            pe = (Ponon^p * Poffoff^(1-p))^N;
            if rand < pe
                cont = false;
                t_av(n+1) = t_av(n+1) + t/trials;
                nout = round(N*p);
                cntmap(n+1, nout + 1) = cntmap(n+1, nout + 1) + 1/trials;
            else
                t = t + 1;
                p = p + normrnd(dp, sqrt(vardp));
                p = max(0,min(p,1));
                theta = theta + normrnd(dtheta, sqrt(vartheta));
                theta = min(theta, fN);
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

save(fullfile(data_path,'pin_pout',fname));

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

