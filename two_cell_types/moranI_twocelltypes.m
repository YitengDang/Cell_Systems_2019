function [theta_11, theta_12, theta_22] = moranI_twocelltypes(cells, idx1, idx2, dist)
% This function calculates the modified Moran's I for two cell types, which
% is a spatial correlation function between the cells of the types
%cells = zeros(N,1);
%cells(randperm(N, round(0.3*N)))=1;

cells_pm = 2*cells - 1; % ON cells: 1, OFF cells: -1
%cells_type1 = cells_pm(idx1);
%cells_type2 = cells_pm(idx2);

%cell_mean = [mean(cells_type1) mean(cells_type2)]; % <Xi>

% meshgrid replicates the cells states in x and y directions
[cells_matx, cells_maty] = meshgrid(cells_pm, cells_pm);

% For the Moran I we use the factor of the diffusion equation as weight
% Account for self-influence
idx = dist>0;
M = zeros(size(dist));

% Matrix of cell reading
M(idx) = exp(-dist(idx))./dist(idx);
%%
% Sum of all weights
w_sum_11 = sum(sum(M(idx1, idx1) ));
w_sum_22 = sum(sum(M(idx2, idx2) ));
w_sum_12 = sum(sum(M(idx1, idx2) ));
% w_sum_21 = sum(sum(M(idx2, idx1) )); % assume same as w_sum_12

%%
theta_11 = sum(sum(M(idx1, idx1).*cells_matx(idx1, idx1).*cells_maty(idx1, idx1)))/w_sum_11;
theta_22 = sum(sum(M(idx2, idx2).*cells_matx(idx2, idx2).*cells_maty(idx2, idx2)))/w_sum_22;
theta_12 = sum(sum(M(idx1, idx2).*cells_matx(idx1, idx2).*cells_maty(idx1, idx2)))/w_sum_12;
% theta = sum(sum(M.*cells_matx.*cells_maty))/w_sum;

%%
%{
% Multiply every weight by the cell state in a crossed fashion
tmp_11 = sum(sum(M(idx1,idx1).*(cells_matx(idx1,idx1) - cell_mean(1)).*(cells_maty(idx1,idx1)-cell_mean(1))));
%tmp_12 = sum(sum(M(idx1,idx2).*(cells_matx(idx1,idx2) - cell_mean(1)).*(cells_maty(idx1,idx2)-cell_mean(2))));
tmp_12 = sum(sum(M(idx1,idx2).*(cells_matx(idx1,idx2) - cell_mean(2)).*(cells_maty(idx1,idx2)-cell_mean(1))));
tmp_22 = sum(sum(M(idx2,idx2).*(cells_matx(idx2,idx2) - cell_mean(2)).*(cells_maty(idx2,idx2)-cell_mean(2))));

cells_var_type1 = var(cells_type1, 1);
cells_var_type2 = var(cells_type2, 1);

if tmp_11==0
    I_11 = 0;
else
    % cells state variance (using N as normalization)
    I_11 = sum(sum(tmp_11))/w_sum_11/cells_var_type1;
end

if tmp_22==0
    I_22 = 0;
else
    % cells state variance (using N as normalization)
    I_22 = sum(sum(tmp_22))/w_sum_22/cells_var_type2;
end

if tmp_12==0
    I_12 = 0;
else
    % cells state variance (using N as normalization)
    I_12 = sum(sum(tmp_12))/w_sum_12/sqrt(cells_var_type1)/sqrt(cells_var_type2);
end
%}
%if tmp==0
%    I = 0;
%else
%    % cells state variance (using N as normalization)
%    cells_var = var(cells_pm,1);
%    I = sum(sum(tmp))/w_sum/cells_var;
%end