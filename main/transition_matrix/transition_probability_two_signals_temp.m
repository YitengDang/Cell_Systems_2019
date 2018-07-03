% Calculate the estimated transition probability of a certain number of 
% cells activating or deactivating for a specified initial number of ON cells
% Extension: for two signals

% Speed up suggestion: save and load nchoosek data 
clear variables
close all
set(0, 'defaulttextinterpreter', 'latex');

% Parameters of the system
gridsize = 8;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;

M_int = [1 -1; 1 0];
Con = [35 15];
K = [8 25; 35 0];
lambda = [1 1.2];

% use hexagonal lattice
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);

% Calculate interaction strength
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = zeros(1, 2);
gN = zeros(1, 2);
for i=1:2
    fN(i) = sum(sinh(Rcell) * sum(exp((Rcell-r)./lambda(i)).*(lambda(i)./r)) ); % calculate signaling strength
    gN(i) = sum(sinh(Rcell)^2 * sum(exp(2*(Rcell-r)./lambda(i)).*(lambda(i)./r).^2 ) ); % calculate noise variance strength
end
%% Test 
%n = 0:N;
n = [0 0];
p = n/N;
p = [0.2 0.6];

mu_nei = fN.*(p.*Con+(1-p)); % neighbour contributions
mu_off_diag = [0 Con(2)*p(2)+(1-p(2)); Con(1)*p(1)+(1-p(1)) 0]; %off-diagonal (mean-field) self-contributions
muON_mat = diag(Con) + mu_off_diag + repmat(mu_nei, 2, 1);
muOFF_mat = diag([1 1]) + mu_off_diag + repmat(mu_nei, 2, 1);

%muON = Con + fN.*(p.*Con+(1-p));
%muOFF = 1 + fN.*(p.*Con+(1-p));
sigmap = sqrt(p.*(1-p).*gN).*(Con-1);
sigmap_mat = repmat(sigmap, 2, 1);

% long version
%{
ponon_mat = zeros(2);
ponon_mat(1,:) = M_int(1,:).*normcdf(M_int(1,:).*(-K(1,:) + muON)./sigmap) + (1-abs(M_int(1,:)));
ponon_mat(2,:) = M_int(2,:).*normcdf(M_int(2,:).*(-K(2,:) + muON)./sigmap) + (1-abs(M_int(2,:)));
ponon = prod(ponon_mat, 2);
%}

% vectorized
%ponon0 = prod(M_int.*normcdf(M_int.*(-K + muON_mat)./sigmap_mat) + (1-abs(M_int)), 2)';
%poffoff0 = prod(M_int.*normcdf(M_int.*(-K + muOFF_mat)./sigmap_mat) + (1-abs(M_int)), 2)';

ponon = prod(1 - abs(M_int).*normcdf(0, M_int.*(muON_mat - K), sigmap_mat), 2 )';
poffon = prod(1 - abs(M_int).*normcdf(0, M_int.*(muOFF_mat - K), sigmap_mat), 2)';
poffoff = 1 - poffon;

pEn0 = ponon.^(n').*poffoff.^(N-n');
pEn = ponon.^(n').*poffoff.^(N-n');

% Calculate pt, ptn
ptn = zeros(N+1);
yminv1 = 0:n(1);
yminv2 = 0:n(2);
yplusv1 = 0:(N-n(1));
yplusv2 = 0:(N-n(2));
pt = zeros(numel(yminv1), numel(yplusv1), numel(yminv2), numel(yplusv2));
for ymin1 = yminv1
    for yplus1 = yplusv1
        for ymin2 = yminv2
            for yplus2 = yplusv2
                %disp([ymin1 yplus1 ymin2 yplus2]);
                ymin = [ymin1 ymin2];
                yplus = [yplus1 yplus2];
                %tmp = zeros(2, 1);
                tmp = ponon.^(n-ymin).*(1-ponon).^ymin.*poffoff.^(N-n-yplus).*(1-poffoff).^yplus;
                tmp(1) = tmp(1)*nchoosek(n(1),ymin1)*nchoosek(N-n(1),yplus1);
                tmp(2) = tmp(2)*nchoosek(n(2),ymin2)*nchoosek(N-n(2),yplus2);   
                if sum(tmp>1)>0
                    disp(tmp)
                end
                pt(ymin1+1,yplus1+1,ymin2+1,yplus2+1) = prod(tmp);
                idx = n + yplus - ymin;
                ptn(idx(1)+1, idx(2)+1) = ptn(idx(1)+1, idx(2)+1) + prod(tmp);
            end
        end
    end
end
%}
%%
%
h = figure(1);
pv = (0:N)/N;
data = squeeze(real( ptn ));
data(data < 10^(-50)) = 0; % neglecting small values

imagesc(pv, pv, data')
c = colorbar;
set(gca, 'Ydir', 'normal', 'FontSize', 20)
xlabel('$$p_1(t+1)$$', 'FontSize', 24)
ylabel('$$p_2(t+1)$$', 'FontSize', 24)
c.Label.String = 'Probability';
%}
%% Full calculation
t_mat = zeros(N+1, N+1, N+1, N+1);
for n1 = 0:N
    for n2 = 0:N
        n = [n1 n2];
        disp(n)
        ptn = zeros(N+1);

        %--- calculate ponon, poffoff
        p = n/N;
        mu_nei = fN.*(p.*Con+(1-p)); % neighbour contributions
        mu_off_diag = [0 Con(2)*p(2)+(1-p(2)); Con(1)*p(1)+(1-p(1)) 0]; %off-diagonal (mean-field) self-contributions
        muON_mat = diag(Con) + mu_off_diag + repmat(mu_nei, 2, 1);
        muOFF_mat = diag([1 1]) + mu_off_diag + repmat(mu_nei, 2, 1);

        sigmap = sqrt(p.*(1-p).*gN).*(Con-1);
        sigmap_mat = repmat(sigmap, 2, 1);

        ponon = prod(1 - abs(M_int).*normcdf(0, M_int.*(muON_mat - K), sigmap_mat), 2 )';
        poffon = prod(1 - abs(M_int).*normcdf(0, M_int.*(muOFF_mat - K), sigmap_mat), 2)';
        poffoff = 1 - poffon;

        %------
        yminv1 = 0:n(1);
        yminv2 = 0:n(2);
        yplusv1 = 0:(N-n(1));
        yplusv2 = 0:(N-n(2));
        
        pt = zeros(numel(yminv1), numel(yplusv1), numel(yminv2), numel(yplusv2));
        for ymin1 = yminv1
            for yplus1 = yplusv1
                for ymin2 = yminv2
                    for yplus2 = yplusv2
                        %disp([ymin1 yplus1 ymin2 yplus2]);
                        ymin = [ymin1 ymin2];
                        yplus = [yplus1 yplus2];
                        %tmp = zeros(2, 1);
                        
                        tmp = ponon.^(n-ymin).*(1-ponon).^ymin.*poffoff.^(N-n-yplus).*(1-poffoff).^yplus;
                        tmp(1) = tmp(1)*nchoosek(n(1),ymin1)*nchoosek(N-n(1),yplus1);
                        tmp(2) = tmp(2)*nchoosek(n(2),ymin2)*nchoosek(N-n(2),yplus2);
                        if sum(tmp>1)>0
                            disp(tmp)
                        end
                        
                        pt(ymin1+1,yplus1+1,ymin2+1,yplus2+1) = prod(tmp);
                        idx = n + yplus - ymin;
                        
                        ptn(idx(1)+1, idx(2)+1) = ptn(idx(1)+1, idx(2)+1) + prod(tmp);
                    end
                end
            end
        end
        t_mat(n1+1, n2+1, :, :) = ptn;
        %----------------------
    end
end

%%
n = [5 2];
p=n/N; 

h=figure(2);
pv = (0:N)/N;
data = squeeze(real(t_mat(n(1)+1,n(2)+1,:,:) ))';
data(data < 10^(-50)) = 0; % neglecting small values

imagesc(pv, pv, data)
c = colorbar;
set(gca, 'Ydir', 'normal', 'FontSize', 20)
xlabel('$$p_1(t+1)$$', 'FontSize', 24)
ylabel('$$p_2(t+1)$$', 'FontSize', 24)
c.Label.String = 'Probability';
%}