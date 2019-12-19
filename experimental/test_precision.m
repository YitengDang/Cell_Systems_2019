clear all;

Digits_default = 32;
gz = 15;
mcsteps = 0;
rcell = 0.2;

% Generate lattice
nodisplay = 1;
%dig = 16;
%digits(16);
[positions, distances] = initial_cells_random_markov_periodic(gz, mcsteps,...
    rcell, nodisplay);

%Increase precision
%digitsOld = digits(50);
Digits = 36;
[positions_vpa, distances_vpa] = initial_cells_random_markov_periodic_vpa(gz, mcsteps,...
    rcell, nodisplay, Digits);
digits(Digits_default);

%% Compare results directly
dpos = double(positions - positions_vpa);
ddist = double(distances - distances_vpa);
max(max(dpos))
max(max(ddist))

%% Compare generating and re-rounding results
dpos = (positions - double(positions_vpa));
ddist = (distances - double(distances_vpa));
max(max(dpos))
max(max(ddist))

