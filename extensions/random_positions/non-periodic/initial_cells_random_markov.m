function [pos, dist] = initial_cells_random_markov(n, L, R)

% Note: even perfect arrangement might not be accepted because distances
% are rounded off.
% Places cells randomly in a continuous space
% Parameters
%clear all
%L = 1;
%R = 0.02; 
%n = round(L/R/10); % nmax = L/R
N = n^2;
%digits(32); %default 32
 
% hexagonal placement with square packing 
delx = L/n;
dely = sqrt(3)/2*delx;
[xm, ym] = meshgrid(0:n-1, 0:n-1);
x = (xm+mod(ym,2)/2)*delx;
y = ym*dely;
xn = x;
yn = y;

% Calculate distance
[x1,x2] = meshgrid(xn, xn);
[y1,y2] = meshgrid(yn, yn);
dx = x1-x2;
dy = y1-y2;
dist = sqrt(dx.^2 + dy.^2);
    
mcsteps = 0; %10^4; % Monte Carlo steps
rejections = 0;
for i=1:mcsteps
    % Random step
    j = randi(N); %selected particle
    sigma = R/10;
    xn(j) = xn(j) + sigma*randn();
    yn(j) = yn(j) + sigma*randn();
    
    % Calculate distance
    [x1,x2] = meshgrid(xn, xn);
    [y1,y2] = meshgrid(yn, yn);
    dx = x1-x2;
    dy = y1-y2;
    dist = sqrt(dx.^2 + dy.^2);

    % Check distance
    if all(dist(dist>0) >= 2*R) 
        x = xn;
        y = yn;
    else
        rejections = rejections + 1;
        %fprintf('Rejected! \n');
    end
end
fprintf('MC steps: %d \n, Rejections: %d \n', mcsteps, rejections);
pos = [x y];

end
%% Draw configuration
%{
hin = figure();
clf(hin,'reset');
title(sprintf('N = %d, R = %.2f', N, R), ...
    'FontSize', 24);
%set(gca,'YTick',[],'XTick',[]);
set(gca,'DataAspectRatio', [1 1 1]);
axis([-R L+R -R L+R]);
box on
hold on
for i=1:N
    position = [x(i)-R y(i)-R 2*R 2*R];
    face_clr = 'k';
    curv = [1 1];
    rectangle('Position', position, 'FaceColor', face_clr, ...
        'EdgeColor', 'k', 'Curvature', curv);
end
hold off
drawnow;

%% Test whether distribution is random enough
dist2 = reshape(dist(dist>0), [N-1 N]);
nnd = min(dist2);
nnd_av = mean(nnd);

h=figure();
hold on
histogram(nnd, 'Normalization', 'pdf');
xlabel('NND');
ylabel('Probability');
title(sprintf('%s NND %s = %.2f', '\langle', '\rangle', nnd_av));
set(gca, 'FontSize', 24);
set(h, 'Position', [500 200 840 720]);

% Clark & Evans calculation
rho = @(r) 2*pi*N/L^2*r.*exp(-pi*N/L^2*r.^2);
rvals = linspace(0, max(nnd)*1.2, 1000);
plot(rvals, rho(rvals), 'LineWidth', 2);
%xlim([0 1]);
%}

% set image properties
%h = gcf;
%set(h,'Units','px');
%set(h, 'Position', [500 300 600 600]);
%}