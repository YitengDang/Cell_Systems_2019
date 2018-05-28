% Computes the eigenvalues of the Jacobian matrix on the invariant manifold {Xi = Xm | for all i}
% This shows whether the fixed points are stable
close all
clear all
warning off
%% Parameters
%Lattice parameters
gridsize = 11;
N = gridsize^2;
a0 = 5;
Rcell = 0.2*a0;

% circuit parameters
K = 16;
Con = 12;
n = 2; % Hill coefficient
%% 
% calculate fN
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist = round(dist, 5);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% calculate M
idx = dist>0;
M = ones(size(dist)); 
M(idx) = sinh(Rcell)./(a0*dist(idx)).*exp(Rcell-a0*dist(idx));
%% Try to derive condition from matrix M only
%{
[V, D] = eig(M);
d = diag(D);
Xi = 0.5;
Yi = (1+fN)*( (Con-1)*Xi + 1);
C = n*Yi^(n-1)*K^n/((K^n+Yi^n)^2).*(Con-1);
max(d);
%}
%% Fixed points
% Plot graphically
xlist = 0:0.005:1;
ylist =  (1+fN)*( (Con-1)*xlist + 1);
figure();
hold on
plot(xlist, ylist.^n./(K^n+ylist.^n), 'LineWidth', 1.5);
plot(xlist, xlist, 'LineWidth', 1.5);
%plot(xlist, n*ylist.^(n-1)*K^n./((K^n+ylist.^n).^2));
xlabel('$$X_i(t)$$');
ylabel('$$X_i(t+1)$$');
set(gca,'FontSize', 24);
set(gcf, 'Units', 'Inches', 'Position', [0 0 10 8]);

% Calculate fixed points on invariant manifold
syms x
f = ((1+fN)*((Con-1)*x+1)).^n./(((1+fN)*((Con-1)*x+1)).^n + K^n) - x;
%start = 0.6;
sol = double(vpasolve(f, x));
disp('Fixed points:');
disp(sol)
%% Calculate eigenvalues of A = Df(x)
large_eval = cell(numel(sol),1);

%for i=1:numel(sol)
    i=1;
    fprintf('Solution %d \n', i);
    Xi = sol(i);
    Yi = (1+fN)*( (Con-1)*Xi + 1);
    A = n*Yi^(n-1)*K^n/((K^n+Yi^n)^2).*((Con-1)*M); % - diag(ones(N,1));
    [evectors, evals] = eig(A);
    disp('Eigenvalues |lambda| >= 1');
    idx = abs(evals)>=1;
    [idx2, ~] = find(idx); % linear index along diag
    disp( evals(idx) ); % Show eigenvalues with |lambda| > 1.
    large_eval{i}{end+1} = evectors(:, idx2);
    %fprintf('Max ev. = %.4f \nmin ev. = %.4f \n', max(diag(evals)), min(diag(evals)));%end
%end
%%
% Check solutions
%{
idx = 1;
x = sol(idx);
%((1+fN)*(Con-1)*x+1).^n./((1+fN)*((Con-1)*x+1).^n + K^n) - x

%}
% Check eigenvalues
%{
check = zeros(N);
for i=1:N
    check(:, i) = (A - evals(i, i)*eye(N))*evectors(:,i);
end
%}
%%
% Check matrix A at random points
%{
for i=0:10
    fprintf('x = %d/10 \n', i);
    Xi = i/10; %sol(i);
    Yi = (1+fN)*( (Con-1)*Xi + 1);
    A = n*Yi^(n-1)*K^n/((K^n+Yi^n)^2).*((Con-1)*M+1); % - diag(ones(N,1));
    [evectors, evals] = eig(A);
    disp('Eigenvalues |lambda| >= 1');
    idx = abs(evals)>=1;
    [idx2, ~] = find(idx); % linear index along diag
    disp( evals(idx) ); % Show eigenvalues with |lambda| > 1.
    %large_eval{i}{end+1} = evectors(:, idx2);
    %fprintf('Max ev. = %.4f \nmin ev. = %.4f \n', max(diag(evals)), min(diag(evals)));%end
end
%}