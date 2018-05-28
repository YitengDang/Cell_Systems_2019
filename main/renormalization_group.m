close all
clear variables
warning off

% Parameters of the system
gridsize = 16;
N = gridsize^2;
K = 10;
Con = 30;
clesz = 2; % cluster edge size; cluster size = d^2
a0 = 1.5;
phi = 0.2; % R = phi*a0
Rcell = phi*a0;

% Get initial distances and fN
[pos,ex,ey] = init_cellpos_hex(gridsize,gridsize);
dist = dist_mat(pos,gridsize,gridsize,ex,ey);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sinh(Rcell)*sum(exp(Rcell-r)./r); % calculate signaling strength

% Get RG flowed distances
a0_RG = clesz*a0;
Rcell_RG = phi*a0_RG;
gridsize_RG = gridsize/clesz;
N_RG = N/clesz^2;
[pos_RG,ex,ey] = init_cellpos_hex(gridsize_RG, gridsize_RG);
dist_RG = dist_mat(pos_RG,gridsize_RG,gridsize_RG,ex,ey);

%% Get indices of clustered cells
s = [gridsize, gridsize];
[x,y] = ind2sub(s, (1:N)');

clusterind = zeros(gridsize/clesz, gridsize/clesz, clesz^2);
for i=1:gridsize/clesz
    for j=1:gridsize/clesz
        a = 2*i-1;
        b = 2*j-1;
        idx = find(x==a & y==b);
        %disp(idx);
        idx2 = idx+1;
        idx3 = idx+gridsize;
        idx4 = idx+gridsize+1;
        %fprintf('(x, y) = (%d, %d) \n', x(idx), y(idx) );
        %fprintf('(x, y) = (%d, %d) \n', x(idx2), y(idx2) );
        %fprintf('(x, y) = (%d, %d) \n', x(idx3), y(idx3) );
        %fprintf('(x, y) = (%d, %d) \n', x(idx4), y(idx4) );
        clusterind(i,j, :) = [idx idx2 idx3 idx4];
    end
end

%cells = zeros(N,1);
%cells(clusterind(7,6,:) ) = 1;

%% Calculate f(r_ij) matrix
% Original system
idx = dist>0;
M = zeros(size(dist)); 
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

% RG system by summing over original system
M_RG2 = zeros(size(dist_RG)); 
all_indices = [];
for c1=1:(gridsize_RG)^2
    [i1, j1] = ind2sub([gridsize_RG, gridsize_RG], c1);
    for c2=1:(gridsize_RG)^2
        [i2, j2] = ind2sub([gridsize_RG, gridsize_RG], c2);
        if c1==c2
            % self-interaction strength
            if c1==1 % calculate only once
                ind_cl1 = clusterind(i1,j1,:); %cluster 1 indices
                dist_vec = dist(ind_cl1(1), ind_cl1([2 3 4]));
                f_temp = sinh(Rcell)*exp(Rcell-a0*dist_vec)./(a0*dist_vec);
                %disp(f_temp);
                M_RG2(c1,c2) = sum(sum(f_temp));
            else
                M_RG2(c1,c2) = M_RG2(1,1); % diagonal entries the same
            end
        else
            if c2 > c1 % only calculate lower diagonal
                ind_cl1 = clusterind(i1,j1,:); %cluster 1 indices
                ind_cl2 = clusterind(i2,j2,:); %cluster 2 indices
                dist_vec = dist(ind_cl1, ind_cl2);
                
                f_temp = sinh(Rcell)*exp(Rcell-a0*dist_vec)./(a0*dist_vec);
                %disp(f_temp);
                % opt (1) average interaction strength over cells in cluster
                %M_RG2(c1,c2) = mean(sum(f_temp, 2), 1); 
                
                % opt (2) take only first cell
                f_temp = sum(f_temp, 2);
                M_RG2(c1,c2) = f_temp(1);
                
                %hin=figure(1);
                %cells = zeros(N,1);
                %cells(ind_cl1) = 1;
                %cells(ind_cl2) = 0.5;
                %cells = ones(N,1);
                %cell_type = zeros(N,1);
                %update_cell_figure_continuum(hin, pos, a0, cells, cell_type, 0);
                %k=waitforbuttonpress;
            else
                M_RG2(c1,c2) = M_RG2(c2,c1);
            end
        end
    end
end

fN_RG = sum(M_RG2, 2);

%% Estimate parameters for best fit
% Interaction matrix 
npoints = 100;
x_all = linspace(0, 5, npoints);
fN_all = zeros(npoints, 1);
for i=1:npoints
    x=x_all(i);
    idx = dist_RG>0;
    M_RG = ones(size(dist_RG)); 
    M_RG(idx) = sinh(Rcell_RG)./(x*a0_RG*dist_RG(idx)).*exp(Rcell_RG-x*a0_RG*dist_RG(idx));
    fN_all(i) = sum(M_RG(1,:));
end

figure();
hold on
plot(x_all, fN_all);
plot([0 2], [1+fN 1+fN]);
ylim([0 10]);

%% Plot summed and from dist M values together
x = 0.7;
idx = dist_RG>0;
M_RG = zeros(size(dist_RG)); 
M_RG(idx) = sinh(Rcell_RG)./(x*a0_RG*dist_RG(idx)).*exp(Rcell_RG-x*a0_RG*dist_RG(idx));
    
figure();
hold on
loglog(clesz^2*M_RG(1,:), M_RG2(1,:), 'bx');
%loglog([0 0.2], [0 0.2]);
%% Plot
%{
cells = zeros(N,1);
cells(ind_cl1) = 1;

hin=figure(1);
hold on 
%cells = ones(N,1);
cell_type = zeros(N,1);
update_cell_figure(hin, pos, a0, cells, cell_type, 0);
%%
hin=figure(2);
cells2 = zeros(N_RG,1);
cell_type2 = zeros(N_RG,1);
update_cell_figure(hin, pos_RG, a0, cells2, cell_type2, 0);
%}