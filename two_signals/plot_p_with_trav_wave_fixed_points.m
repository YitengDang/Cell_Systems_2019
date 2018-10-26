%% Plot p(t) together with the travelling wave fixed points
clear all
close all
clc
%%
WxWy_all = [1 0; 0 1; 1 1; 1 2; 2 1];
nc_all = [1 1 3/2 2 5/2];
num_types = numel(nc_all);

% Parameters
gz = 16;

% Plot
WxWy_str = {};
h = figure;
cmap = colormap(copper(num_types));
markers = {'o', 'd'};
hold on
for nwaves = 1:2
    p_all = 2*nc_all/gz*nwaves;
    
    % scatter 
    sz = 150; % default=36
    for i=1:num_types
        scatter(p_all(i), p_all(i), sz, cmap(i,:), markers{nwaves}, 'filled');
    end
    
    % update legend
    WxWy_str = sprintfc('(Wx, Wy) = (%d %d)', WxWy_all);
    WxWy_str = [WxWy_str; WxWy_str];
end
legend(WxWy_str);
xlim([0 1]);
ylim([0 1]);