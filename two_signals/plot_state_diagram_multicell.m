clear all;
close all;
set(0, 'defaulttextinterpreter', 'latex');
%% Input parameters
% lattice parameters
gz = 15;
N = gz^2;
a0 = 1.5;
rcell = 0.2;
Rcell = rcell*a0;

% circuit parameters
Con = [18 16];
Coff = [1 1];
M_int = [1 1; -1 -1];
K = [3 12; 13 20]; % K(i,j): sensitivity of type i to type j molecules
lambda = [1 1.2]; % diffusion length (normalize first to 1)
hill = Inf;
noise = 0;

% calculate fN
[dist, ~] = init_dist_hex(gz, gz);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = zeros(2,1);
fN(1) = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(1)).*(lambda(1)./r)) ); % calculate signaling strength
fN(2) = sum(sinh(Rcell)*sum(exp((Rcell-r)./lambda(2)).*(lambda(2)./r)) ); % calculate signaling strength


%% Calculate transition tables
% Determine phases
R1 = (repmat(1+fN', 2, 1) - K) > 0; % Everything ON
R2 = ((repmat(Con + fN', 2, 1) - K) > 0 & (repmat(1+fN', 2, 1) - K) < 0); % ON remains ON & not all ON
R3 = ((1 + repmat(fN'.*Con, 2, 1) - K) < 0 & (repmat((1+fN').*Con, 2, 1) - K) > 0) ; % OFF remains OFF & not all OFF
R4 = (repmat((1+fN').*Con, 2, 1) - K) < 0; % Everything OFF

phase = R1 + 2*R2 + 3*R3 + 4*R4;
% 0: none (activation-deactivation)
% 1: all ON(+) / OFF(-)
% 2: ON->ON (+) / ON-> OFF (-)
% 3: OFF->OFF (+) / OFF->ON (-)
% 4: all OFF(+) / ON(-)
% 5: autonomy (+) / autonomous oscillations (-)

% Map from phase to diagram
% state | activation/repression | input molecule (1/2)
g_map = cell(2, 6, 2);
% 0=OFF, 1:ON, 2:UNKNOWN
% activation 
g_map{1,1,1} = 2*ones(2);
g_map{1,1,2} = 2*ones(2);
g_map{1,2,1} = ones(2);
g_map{1,2,2} = ones(2);
g_map{1,3,1} = [2 2; 1 1];
g_map{1,3,2} = [2 1; 2 1];
g_map{1,4,1} = [0 0; 2 2];
g_map{1,4,2} = [0 2; 0 2];
g_map{1,5,1} = zeros(2);
g_map{1,5,2} = zeros(2);
g_map{1,6,1} = [0 0; 1 1];
g_map{1,6,2} = [0 1; 0 1];
% repression 
%(note: this is precisely NOT g_map{1,:,:} in the three-val
% boolean algebra with NOT 2 = 2)
g_map{2,1,1} = 2*ones(2);
g_map{2,1,2} = 2*ones(2);
g_map{2,2,1} = zeros(2);
g_map{2,2,2} = zeros(2);
g_map{2,3,1} = [2 2; 0 0];
g_map{2,3,2} = [2 0; 2 0];
g_map{2,4,1} = [1 1; 2 2];
g_map{2,4,2} = [1 2; 1 2];
g_map{2,5,1} = ones(2);
g_map{2,5,2} = ones(2);
g_map{2,6,1} = [1 1; 0 0];
g_map{2,6,2} = [1 0; 1 0];

gij = cell(2);
X_out = cell(2, 1);
for i=1:2
    for j=1:2
        if M_int(i,j)~=0
            idx = (M_int(i,j)==1) + (M_int(i,j)==-1)*2;
            gij{i,j} = g_map{idx, phase(i,j)+1, j};
        else
            gij{i,j} = ones(2); % Fix
        end
    end
    X_out{i} = gij{i,1} & gij{i,2}; % Fix 3-valued algebra
end

% Display tables
h1 = figure(1);
subplot(1, 2, 1);
imagesc([0 1], [0 1], X_out{1})
set(gca, 'YDir', 'normal');
xticks([0 1]);
yticks([0 1]);
set(gca, 'FontSize', 24);
title('$$X^{(1)}_{out}$$ ')
xlabel('$$X^{(1)}_{in}$$ ')
ylabel('$$X^{(2)}_{in}$$ ')
ncolor = numel(unique(X_out{1}));
if ncolor>1
    colormap(parula(ncolor));
    colorbar;
end

subplot(1, 2, 2);
imagesc([0 1], [0 1], X_out{2})
set(gca, 'YDir', 'normal');
xticks([0 1]);
yticks([0 1]);
set(gca, 'FontSize', 24);
title('$$X^{(2)}_{out}$$ ')
xlabel('$$X^{(1)}_{in}$$ ')
ylabel('$$X^{(2)}_{in}$$ ')

set(h1, 'Units', 'Inches', 'Position', [1 1 11 5]);
ncolor = numel(unique(X_out{1}));
if ncolor>1
    colormap(parula(ncolor));
    colorbar;
end

%% Calculate state diagram

A = zeros(4); % graph adjacency matrix
% Add deterministic transitions
for i=1:2
    for j=1:2
        state_in = i + 2*(j-1); 
        X_out_this = [X_out{1}(i,j) X_out{2}(i,j)]; % tentative 
        disp(X_out_this)
        if all(X_out_this~=2) % unambiguous out state
            state_out = X_out_this(1)+1 + 2*X_out_this(2); % (i,j) -> idx
            disp(state_out);
            A(state_in, state_out) = 1;
        end
    end
end

% Add semi-determined transitions


%% Draw state diagram
h2 = figure(2);
hold on
s = [0 1 0 1];
t = [0 0 1 1];
%A = ones(4);
Gs = digraph(A);
nLabels = {};
plot(Gs, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
    'LineWidth', 3, 'EdgeColor', 'k',...
    'Marker', 'o', 'MarkerSize', 100, 'NodeColor', [0.2 0.2 0.2], 'NodeLabel', nLabels);
text(s-0.09,t+0.015,{'(0,0)','(1,0)','(0,1)','(1,1)'}, 'Color', 'w', 'FontSize', 32)
ax = gca;
axis([-0.4 1.4 -0.4 1.4]);
ax.Visible = 'off';
h2.Color = [1 1 1];
set(ax, 'Units', 'Inches', 'Position', [0 0 9 8]);
set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);

%% function plot_state_diagram_multicell(M_int, Con, Coff, K)

% end