% Creates pin-pout map from (p,I) simulation (analytical = ?)

clear variables
close all

data_path = 'D:\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy';

% Get phenotype regions
gridsize = 25;
a0 = 0.5;
Rcell = 0.2*a0;
Con = 5;
K = 10;
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
        dx = [1 1];
        while cont
            muon = fN*(Con*p + 1 - p + (Con-1)*(1-p)*I);
            muoff = fN*(Con*p + 1 - p - (Con-1)*p*I);
            
%             kappa = (Con-1)^2*(gN)*p*(1-p)/ ...
%                 (p*muon + (1-p)*muoff);
%             sigmaon = sqrt(max(kappa*muon,0)+alpha^2);
%             sigmaoff = sqrt(max(kappa*muoff,0)+alpha^2);
            if I >=1
                sigma_aux = 0;
            else
                sigma_aux = max((Con-1)^2*p*(1-p)*(gN),0);
            end
            sigmaon = sqrt(sigma_aux + alpha^2);
            sigmaoff = sigmaon;
            zon = (K - Con - muon)/sigmaon;
            zoff = (K - 1 - muoff)/sigmaoff;

            Poffoff = normcdf(zoff);
            Ponon = 1-normcdf(zon);
            
            pe = (Ponon^p * Poffoff^(1-p))^N;
            if rand < pe
                cont = false;
                t_av(n+1) = t_av(n+1) + t/trials;
                nout = round(N*p);
                cntmap(n+1, nout + 1) = cntmap(n+1, nout + 1) + 1/trials;
            else
                t = t + 1;
                
                p_minus = p*(1-Ponon);
                p_plus = (1-p)*(1-Poffoff);
                
                dp = p_plus - p_minus; % Correct for discreteness
                
                var_p_plus = p_plus*Poffoff/N;
                var_p_minus = p_minus*Ponon/N;
                vardp = var_p_plus + var_p_minus;
                                
                if 1 - Poffoff < 1e-16
                    y_plus = muoff;
                    var_y_plus = sigmaoff^2;
                else
                    y_plus = muoff + sigmaoff*phi(zoff)/(1 - Poffoff);
                    if sigmaoff < 1e-16
                        var_y_plus = 0;
                    else
                        var_y_plus = sigmaoff^2*(1 + zoff*phi(zoff)/(1 - Poffoff) ...
                            - (phi(zoff)/(1 - Poffoff))^2);
                    end
                end
                                
                if 1 - Ponon < 1e-16
                    y_minus = muon;
                    var_y_minus = sigmaon^2;
                else
                    y_minus = muon - sigmaon*phi(zon)/(1 - Ponon);
                    if sigmaon < 1e-16
                        var_y_minus = 0;
                    else
                        var_y_minus = sigmaon^2*(1 - zon*phi(zon)/(1 - Ponon) ...
                            - (phi(zon)/(1 - Ponon))^2);
                    end
                end
                
                sigma_minus = p_minus^2*var_y_minus + ...
                    y_minus^2*var_p_minus + var_y_minus*var_p_minus;
                
                sigma_plus = p_plus^2*var_y_plus + ...
                    y_plus^2*var_p_plus + var_y_plus*var_p_plus;

                dtheta = 8/(Con-1)*(p_plus*y_plus - p_minus*y_minus) ...
                    + 4*fN*(dp^2+vardp - dp*(Con+1)/(Con-1));
                
                vardtheta = 64/(Con-1)^2*(sigma_plus + sigma_minus) + ...
                    16*fN^2*vardp*((Con+1)^2/(Con-1)^2 + 4*dp^2 + 2*vardp);
                
                covdpdtheta = 8/(Con-1)*(y_plus*(var_p_plus + p_plus*dp) ...
                    + y_minus*(var_p_minus - p_minus*dp)) ...
                    + 4*fN*(dp^3 + 3*dp*vardp - ...
                            (Con+1)/(Con-1)*(vardp + dp^2)) - dp*dtheta;
                        
                cov_d = [vardp covdpdtheta; covdpdtheta vardtheta];
                mean_d = [dp; dtheta];
                
                [c_aux, r] = chol(cov_d, 'lower');
                    
                if r>0
                    dx = [dp dtheta];
                else
                    dx = mean_d + c_aux*polar_marsaglia;                  
                end
                p = round(N*(p + dx(1)))/N;
                theta = theta + dx(2);
                p = max(0,min(p,1));
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

