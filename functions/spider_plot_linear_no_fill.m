function fig_handle = spider_plot_linear_no_fill(P, P_labels, axes_interval, varargin)
% Create a spider web or radar plot with an axes specified for each column
%
% Custom version Yiteng Dang 2018

% spider_plot(P, P_labels, axes_interval, axes_precision) creates a spider
% web plot using the points specified in the array P. The column of P
% contains the data points and the rows of P contain the multiple sets of
% data points. Each point must be accompanied by a label specified in the
% cell P_labels. The number of intervals that separate the axes is
% specified by axes_interval. The number of decimal precision points is
% specified by axes_precision.
% 
% P - [vector | matrix]
% P_labels - [cell of strings]
% axes_interval - [integer] (number of spider web lines to plot)
% axes_precision - [integer] (number of decimals of axes ticks)
%
% spider_plot(P, P_labels, axes_interval, axes_precision, line_spec) works
% the same as the function above. Additional line properties can be added
% in the same format as the default "plot" function in MATLAB.
%
% line_spec - [character vector]
% 
% %%%%%%%%%%%%%%%%%%% Example of a Generic Spider Plot %%%%%%%%%%%%%%%%%%%
% % Clear workspace
% close all;
% clearvars;
% clc;
% 
% % Point properties
% num_of_points = 6;
% row_of_points = 4;
% 
% % Random data
% P = rand(row_of_points, num_of_points);
% 
% % Scale points by a factor
% P(:, 2) = P(:, 2) * 2;
% P(:, 3) = P(:, 3) * 3;
% P(:, 4) = P(:, 4) * 4;
% P(:, 5) = P(:, 5) * 5;
%
% % Make random values negative
% P(1:3, 3) = P(1:3, 3) * -1;
% P(:, 5) = P(:, 5) * -1;
% 
% % Create generic labels
% P_labels = cell(num_of_points, 1);
% 
% for ii = 1:num_of_points
%     P_labels{ii} = sprintf('Label %i', ii);
% end
% 
% % Figure properties
% figure('units', 'normalized', 'outerposition', [0 0.05 1 0.95]);
% 
% % Axes properties
% axes_interval = 2;
% axes_precision = 1;
% 
% % Spider plot
% spider_plot(P, P_labels, axes_interval, axes_precision,...
%     'Marker', 'o',...
%     'LineStyle', '-',...
%     'LineWidth', 2,...
%     'MarkerSize', 5);
% 
% % Title properties
% title('Sample Spider Plot',...
%     'Fontweight', 'bold',...
%     'FontSize', 12);
% 
% % Legend properties
% legend('show', 'Location', 'southoutside');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Point Properties %%%
% Number of points
[row_of_points, num_of_points] = size(P);

%%% Error Check %%%
% Check if axes properties are an integer
if floor(axes_interval) ~= axes_interval %|| floor(axes_precision) ~= axes_precision
    error('Error: Please enter in an integer for the axes properties.');
end

% Check if axes properties are positive
if axes_interval < 1 %|| axes_precision < 1
    error('Error: Please enter value greater than one for the axes properties.');
end

% Check if the labels are the same number as the number of points
if length(P_labels) ~= num_of_points
    error('Error: Please make sure the number of labels is the same as the number of points.');
end

% Pre-allocation
max_values = zeros(1, num_of_points);
min_values = zeros(1, num_of_points);
axis_increment = zeros(1, num_of_points);

% Normalized axis increment
normalized_axis_increment = 1/(axes_interval+1);

% Iterate through number of points
for ii = 1:num_of_points
    %{
    % Group of points
    group_points = P(:, ii);
    
    % Max and min value of each group
    max_values(ii) = max(group_points);
    min_values(ii) = min(group_points);
    range = max_values(ii) - min_values(ii);
    
    % Axis increment
    axis_increment(ii) = range/axes_interval;
    
    % Normalize points to range from [0, 1]
    P(:, ii) = (P(:, ii)-min(group_points))/range;
    
    % Shift points by one axis increment
    P(:, ii) = P(:, ii) + normalized_axis_increment;
    %}
    
    % set custom range and normalize
    % Max and min value of each group
    max_values(ii) = 1000;
    min_values(ii) = 1;
    
    % Axis increment
    range = max_values(ii) - min_values(ii);
    axis_increment(ii) = range/(axes_interval+1);
    
    P(:, ii) = (P(:, ii)-min_values(ii))/(max_values(ii)-min_values(ii));    
end


%%% Figure Properties %%%
colors = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.9290, 0.6940, 0.1250;...
    0.4940, 0.1840, 0.5560;...
    0.4660, 0.6740, 0.1880;...
    0.3010, 0.7450, 0.9330;...
    0.6350, 0.0780, 0.1840];
%}
% Repeat colors is necessary
repeat_colors = fix(row_of_points/size(colors, 1))+1;
colors = repmat(colors, repeat_colors, 1);

%%% Polar Axes %%%
% Polar increments
polar_increments = 2*pi/num_of_points;

% Normalized max limit of axes
axes_limit = 1;

% Shift axes limit by one axis increment
% axes_limit = axes_limit + normalized_axis_increment;

% Polar points
radius = [0; axes_limit];
theta = 0:polar_increments:2*pi;

% Convert polar to cartesian coordinates
[x_axes, y_axes] = pol2cart(theta, radius);
%%
%%% Data Points %%%
%fig_handle = figure;
hold on
% Iterate through all the rows
%{
for ii = 1:row_of_points
    % Convert polar to cartesian coordinates
    [x_points, y_points] = pol2cart(theta(1:end-1), P(ii, :));
    
    % Make points circular
    x_circular = [x_points, x_points(1)];
    y_circular = [y_points, y_points(1)];
    
    % Plot data points
    p = plot(x_circular, y_circular,...
        'Color', colors(ii, :),...
        'MarkerFaceColor', colors(ii, :),...
        varargin{:});
    p.Color(4) = 0.1; % transparency
end
%}

% (1) Plot lines
%
for ii = 1:row_of_points
    % Convert polar to cartesian coordinates
    [x_points, y_points] = pol2cart(theta(1:end-1), P(ii, :));
    
    % Make points circular
    x_circular = [x_points, x_points(1)];
    y_circular = [y_points, y_points(1)];
    
    % Plot data points
    hold on
    p = plot(x_circular, y_circular,...
        varargin{:});
    p.Color(4) = 0.5; % transparency
end
%}

% (1b) Plot colored areas
%{
P_min = min(P, [], 1);
P_max = max(P, [], 1);

% fill 
[x_points, y_points] = pol2cart(theta(1:end-1), P_max);
p = fill(x_points, y_points, 'r', 'FaceAlpha', 0.2); % 'EdgeColor', 'r',
% empty inside
[x_points, y_points] = pol2cart(theta(1:end-1), P_min);
p = fill(x_points, y_points, 'w', 'FaceAlpha', 1); % 'EdgeColor', 'r',
%}

% (1c) Plot mean values
%{
P_mean = mean(P, 1);
[x_points, y_points] = pol2cart(theta(1:end-1), P_mean);
x_circular = [x_points, x_points(1)];
y_circular = [y_points, y_points(1)];
p = plot(x_circular, y_circular,...
        'Marker', 'none',...
        'Color', [1 0 0],...
        'LineStyle', '-',...
        'LineWidth', 3);
%}

% (2) plot markers
%{
for ii = 1:row_of_points
    % Convert polar to cartesian coordinates
    [x_points, y_points] = pol2cart(theta(1:end-1), P(ii, :));
    
    % Make points circular
    x_circular = [x_points, x_points(1)];
    y_circular = [y_points, y_points(1)];
    
    % Plot data points
    hold on
    p = plot(x_circular, y_circular,...
        'Marker', 'o',...
        'Color', [1 0 0],...
        'MarkerFaceColor', [1 0 0],...
        'LineStyle', 'none',...
        'MarkerSize', 3);
    p.Color(4) = 0.1; % transparency
end
%}
%% Plot axes after data
%
hold on
% Plot polar axes
%{
grey = [1, 1, 1] * 0.5;
h = line(x_axes, y_axes,...
    'LineWidth', 1,...
    'Color', grey);
%}
% Plot polar axes
%
for ii=1:size(x_axes, 2)
    %grey = [1, 1, 1] * 0.5;
    h = line(x_axes(:,ii), y_axes(:,ii),...
        'LineWidth', 1,...
        'Color', 'k');
    h.Color(4) = 0.5; % set each of the axes to transparent
end
%}
% Iterate through all the line handles
for ii = 1:length(h)
    % Remove polar axes from legend
    h(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
end

%%% Polar Isocurves %%%
% Shifted axes interval
shifted_axes_interval = axes_interval+1;

% Incremental radius
radius = (0:axes_limit/shifted_axes_interval:axes_limit)';

% Convert polar to cartesian coordinates
[x_isocurves, y_isocurves] = pol2cart(theta, radius);

% Plot polar isocurves
%
for ii=1:size(x_isocurves, 1)
    h = plot(x_isocurves(ii, :)', y_isocurves(ii, :)',...
        'LineWidth', 1,...
        'Color', 'k');
    h.Color(4) = 0.5; % set each of the axes to transparent
end
%{
hold on;
h = plot(x_isocurves', y_isocurves',...
    'LineWidth', 1,...
    'Color', grey);
%}
%alpha(h, 0.5);

% Iterate through all the plot handles
for ii = 1:length(h)
    % Remove polar isocurves from legend
    h(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
end

%}

%%
%%% Axis Properties %%%
% Figure background
fig = gcf;
fig.Color = 'white';

% Plot text with an offset
r_offset = 0.02; % default 0.05;
[x_isocurves_text, y_isocurves_text] = pol2cart(theta, radius+r_offset);
%{
% Iterate through all the number of points
for hh = 1:num_of_points
    % Shifted min value
    shifted_min_value = min_values(hh)-axis_increment(hh);
    
    % Axis label for each row
    %row_axis_labels = (shifted_min_value:axis_increment(hh):max_values(hh))';
    
    % Don't shift min value:
    row_axis_labels = (min_values(hh):axis_increment(hh):max_values(hh))';
    
    % Iterate through all the isocurve radius
    for ii = length(radius) %2:length(radius)
        % Display axis text for each isocurve
        %text(x_isocurves_text(ii, hh), y_isocurves_text(ii, hh), sprintf(sprintf('$10^{%%.%if}$', axes_precision), row_axis_labels(ii)),...
        text(x_isocurves_text(ii, hh), y_isocurves_text(ii, hh), sprintf('$10^{%d}$', row_axis_labels(ii)),...
            'Units', 'Data',...
            'Color', 'k',...
            'FontSize', 16,...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'middle');
    end
end

% Plot text at center
text(0,0, '1',...
    'Units', 'Data',...
    'Color', 'k',...
    'FontSize', 16,...
    'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'middle');
%}
% Label points
x_label = x_isocurves_text(end, :);
y_label = y_isocurves_text(end, :);

% Shift axis label
shift_pos = 0.02; % default 0.07;

% Iterate through each label
for ii = 1:num_of_points
    % Angle of point in radians
    theta_point = theta(ii);
    
    % Find out which quadrant the point is in
    if theta_point == 0
        quadrant = 0;
    elseif theta_point == pi/2
        quadrant = 1.5;
    elseif theta_point == pi
        quadrant = 2.5;
    elseif theta_point == 3*pi/2
        quadrant = 3.5;
    elseif theta_point == 2*pi
        quadrant = 0;
    elseif theta_point > 0 && theta_point < pi/2
        quadrant = 1;
    elseif theta_point > pi/2 && theta_point < pi
        quadrant = 2;
    elseif theta_point > pi && theta_point < 3*pi/2
        quadrant = 3;
    elseif theta_point > 3*pi/2 && theta_point < 2*pi
        quadrant = 4;
    end
    
    % Adjust text alignment information depending on quadrant
    switch quadrant
        case 0
            horz_align = 'left';
            vert_align = 'middle';
            x_pos = shift_pos;
            y_pos = 0;
        case 1
            horz_align = 'left';
            vert_align = 'bottom';
            x_pos = shift_pos;
            y_pos = shift_pos;
        case 1.5
            horz_align = 'center';
            vert_align = 'bottom';
            x_pos = 0;
            y_pos = shift_pos;
        case 2
            horz_align = 'right';
            vert_align = 'bottom';
            x_pos = -shift_pos;
            y_pos = shift_pos;
        case 2.5
            horz_align = 'right';
            vert_align = 'middle';
            x_pos = -shift_pos;
            y_pos = 0;
        case 3
            horz_align = 'right';
            vert_align = 'top';
            x_pos = -shift_pos;
            y_pos = -shift_pos;
        case 3.5
            horz_align = 'center';
            vert_align = 'top';
            x_pos = 0;
            y_pos = -shift_pos;
        case 4
            horz_align = 'left';
            vert_align = 'top';
            x_pos = shift_pos;
            y_pos = -shift_pos;
    end
    
    % Display text label
    text(x_label(ii)+x_pos, y_label(ii)+y_pos, P_labels{ii},...
        'Units', 'Data',...
        'HorizontalAlignment', horz_align,...
        'VerticalAlignment', vert_align,...
        'BackgroundColor', 'w',...
        'FontSize', 32); 
        %'EdgeColor', 'k',...
end

% Axis limits
axis square;
axis([-axes_limit, axes_limit, -axes_limit, axes_limit]);
axis off;


end