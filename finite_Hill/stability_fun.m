function [stable, eigenvalues] = stability_fun(cells, a0, Rcell, K, Con, n, dist)
% Calculates whether a particular configuration is stable by determining
% eigenvalues of the linearized system's matrix
%clear all
%% Parameters
%Lattice parameters
%gridsize = 11;
%N = gridsize^2;
N = length(dist);
%a0 = 6;
%Rcell = 0.2*a0;

% circuit parameters
%K = 16;
%Con = 16;
%n = 2; % Hill coefficient
%% 
% calculate fN
%[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% calculate M
idx = dist>0;
M = ones(size(dist)); 
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));

%% Initiate cells
%cells = rand(N, 1);

%% Calculate eigenvalues of A = Df(x)
C0 = 1 + (Con-1).*cells; % Concentration sensed by each cell
Y = M*C0; % Reading of each cell
Yi = repmat(Y, 1, N);
A = n*Yi.^(n-1)*K^n./((K^n+Yi.^n).^2).*((Con-1)*M);
[~, evals] = eig(A);
evals = real(evals); % take only real part of complex solutions
disp('Eigenvalues |lambda| >= 1');
idx = abs(evals)>=1;
[idx2, ~] = find(idx); % linear index along diag
disp( evals(idx) ); % Show eigenvalues with |lambda| > 1.
disp('Max eigenvalue:');
disp( max(diag(evals)) ); % Show eigenvalues with |lambda| > 1.
%large_eval{i}{end+1} = evectors(:, idx2);
%fprintf('Max ev. = %.4f \nmin ev. = %.4f \n', max(diag(evals)), min(diag(evals)));%end

% --Output--
stable = isempty(idx2); %0: point unstable, 1: point stable
eigenvalues = diag(evals);
end