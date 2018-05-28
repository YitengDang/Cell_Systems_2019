clear variables

% Plot a line of the entropy map of a saved simulation and compare to the
% analytical formula (K is changed along the line)

fname = fullfile(pwd,'data','20160518_0128_changeK_N11Son40_hexagonal.mat');

load(fname)

% Plot the result
figure(1)
plot(K,omega, 'r-')
hold on
plot(K2,omega_sim, 'ko')
hold off
xlabel('K')
ylabel('Entropy')