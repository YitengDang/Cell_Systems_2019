% Implements the Gillespie (stochastic stimulation) algorithm for diffusion
% on a hexagonal lattice.
 
% 1. Draw random number r
% 2. Propensity function alpha0 = k+gamma*X + 6*d
% 3. Reaction time 1/alpha0*log(r)
% 4. Adjust new concentration based on value of r.
    % -Need matrix of nearest neighbours
    % -Choose values for parameters
%%
clear variables
close all

% initialisation
L = 11;
%[dist, pos] = init_dist_hex(gridsize, gridsize);
xlist = []; 
ylist = [];
x = randi(L); y = randi(L);
xlist(end+1) = x;
ylist(end+1) = y;

% Define functions for moving in positive and negative directions
syms n
plus = @(n) n+1 - (n==L)*L;
minus = @(n) n-1 + (n==1)*L;

%%
% 1. Draw random number r
r = rand();
if (0 < r) && (r < 1/6)
    y = plus(y);
elseif (1/6 < r) && (r < 2/6)
    x = plus(x);
elseif (2/6 < r) && (r < 3/6)
    y = minus(y);
elseif (3/6 < r) && (r < 4/6)
    x = minus(x); 
elseif (4/6 < r) && (r < 5/6)
    x = plus(x); y = minus(y);
elseif (5/6 < r) && (r < 6/6)
    x = minus(x); y = plus(y);
end
xlist(end+1) = x;
ylist(end+1) = y;