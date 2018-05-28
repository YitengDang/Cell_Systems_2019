%% Estimate delta and sigma in Langevin equation
clear all
close all
warning off
%clc
digits(64); %default 32
%%
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Klist =  [10]; % 10 14 19 16 18]; %(a0=1.5) [3 6 10 15 17 20]; % (a0=0.5) [10 10 14 19 16 18]; %[3]; 
Conlist = [5]; %21 16 14 8 6]; %(a0=1.5) [24 21 21 20 14 14]; % (a0=0.5) [5 21 16 14 8 6]; %[24];
deltalist = zeros(numel(Klist), 1);
sigmalist = deltalist;

for idx=1:numel(Klist)
    %idx = 1;
    K = Klist(idx);
    Con = Conlist(idx);

    % use hexagonal lattice
    Rcell = 0.2*a0;
    [dist, pos] = init_dist_hex(gridsize, gridsize);
    dist_vec = a0*dist(1,:);
    r = dist_vec(dist_vec>0); 
    fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); 
    gN = sum(sum((sinh(Rcell)*exp(Rcell-r)./r).^2));

    % Calculate as for loop
    delta = 0;
    sigma = 0;
    I = 0;

    for n=0:N
        p = n/N;
        % Gradient term
        delh_delp = (Con-1)*2*fN*(I-1).*(2*p - 1) -(Con+1)*(fN+1)+2*K;

        % Calculate Ponon, Poffoff => delp
        muon = fN*(Con*p + 1 - p + (Con-1)*(1-p)*I);
        muoff = fN*(Con*p + 1 - p - (Con-1)*p*I);
        kappa = (Con-1)*sqrt(gN*p*(1-p));
        sigmaon = kappa; %sqrt(kappa^2+alpha^2);
        sigmaoff = kappa; %sqrt(kappa^2+alpha^2);

        zon = (K - Con - muon)/sigmaon;
        zoff = (K - 1 - muoff)/sigmaoff;
        Poffoff = normcdf(zoff);
        Ponon = 1-normcdf(zon);

        delp = (1-p)*(1-Poffoff) - p*(1-Ponon);
        prob = nchoosek(N, n)/2^N;
        
        % Update delta, sigma
        delta = delta - prob*delp/delh_delp;
        sigma = sigma + prob*sqrt(p*(1-Poffoff)*Poffoff + p*(1-Ponon)*Ponon)/sqrt(N);
    end
    % Store values
    deltalist(idx) = delta;
    sigmalist(idx) = sigma;
end
%% Old approach
%{
delta_min = zeros(30, 20);
delta_max_p = zeros(30, 20);
for i=1:30
    for j=1:20
        this_Con = i;
        this_K = j;
        p = (0:N)/N;
        I = linspace(-0.1, 0.25, 200);
        [pm, Im] = meshgrid(p, I);
        grad_p = abs( (this_Con-1)*2*fN*(1-Im).*(2*pm-1) + (this_Con+1)*(1+fN) - 2*this_K ); %abs(del h/del p)
        grad_p_inv = grad_p.^(-1);
        idx_min = find( grad_p_inv == min(min(grad_p_inv)), 1);
        %delta_min(i,j) = delta_list(idx_min)/N;
        delta_min(i,j) = median(median(grad_p_inv))/N;
        idx_max = find( grad_p_inv == max(max(grad_p_inv)), 1);
        %delta_max_p(i,j) = delta_list(idx_max);
        delta_max_p(i,j) = median(median(grad_p_inv));        
    end
end
delta_max_I = (fN*(Con-1)/2).^(-1);
delta_max = min(delta_max_p, delta_max_I);
%% Plot
figure(1);
imagesc(1:20, 1:30, delta_min);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('K');
ylabel('C_{ON}');
ylabel(c, '\delta_{min}');
%caxis([0 1]);
set(gca, 'FontSize', 24);
set(gcf, 'Units','Inches', 'Position', [0 0 10 8]);

figure(2);
imagesc(1:20, 1:30, delta_max);
set(gca, 'YDir', 'Normal');
c = colorbar;
xlabel('K');
ylabel('C_{ON}');
ylabel(c, '\delta_{max}');
%caxis([0 1]);
set(gca, 'FontSize', 24);
set(gcf, 'Units','Inches', 'Position', [0 0 10 8]);

%%
test = delta_max >= delta_min;
delta_max_av = mean(mean(delta_max));
delta_min_av = mean(mean(delta_min));
%}