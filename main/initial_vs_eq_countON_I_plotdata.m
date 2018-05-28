clear variables
close all
warning off

% Plot saved data of a simulated pin_pout map and save as PDF. This is done
% for the output of initial_vs_eq_countON_I_batch_save

% Asks the user to browse the .mat file
data_path = '~/Dropbox/Matlab codes/data_onecelltype_entropy';
[fnamein, path, ~] = uigetfile(fullfile(data_path,'pin_pout_noise'));
noise = 0;
load(fullfile(path,fnamein));

% Save Figure?
save_fig = 1;

% This is to be compatible with new saved data that used Con instead of Son
if exist('Con', 'var')
    Son = Con;
end


[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = dist(1,:);
r = a0*dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

[path, name, ~] = fileparts(fnamein);
disp(name);

p = (0:N)./N;
prob = transpose(count./repmat(sum(count,2),1,N+1));

% Plot the map
h1 = figure(1);
im_fig = imagesc(p,(0:N)/N,prob);
% Plot lines
p0 = 0.5-0.25*(0.5*(Son+1)*(1+fN) - K)/fN;
vline_p = [p0];
hold on
for a = vline_p
plot([a a],[0 1], '--r', 'Linewidth', 1.5)
end
hold off
set(gca,'Ydir','normal','FontSize', 20)
set(im_fig, 'AlphaData', count' > 0);
c = colorbar;
c.Label.String = 'Probability';
xlabel('p_{ini}', 'FontSize', 24)
ylabel('p_{eq}', 'FontSize', 24)
% Check whether there is noise or not to figure the data properly
if noise > 0
    title(sprintf('C_{ON}: %.1f, K: %.1f, a0: %.1f, noise: %.2f', ...
        Son, K, a0, noise), 'FontSize', 24)    
else
    title(sprintf('C_{ON}: %.1f, K: %.1f, a0: %.1f', Son, K, a0), ...
        'FontSize', 24)
end
% Organize and save
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout', strcat(name,'_map'));
    save_figure_svg(h1, 10, 6, out_file);
end

% Plot the average number of steps
h2 = figure(2);
plot(p, t_av, 'r-o')
set(gca,'FontSize', 20)
xlabel('p_{ini}', 'FontSize', 24)
ylabel('n_{steps}', 'FontSize', 24)
if noise > 0
    title(sprintf('C_{ON}: %.1f, K: %.1f, a0: %.1f, noise: %.1f', ...
        Son, K, a0, noise), 'FontSize', 24)    
else
    title(sprintf('C_{ON}: %.1f, K: %.1f, a0: %.1f', Son, K, a0), ...
        'FontSize', 24)
end
% Organize and save
if save_fig > 0
    out_file = fullfile(pwd, 'figures', 'pin_pout', strcat(name,'_time'));
    save_figure_pdf(h2, 10, 6, out_file);
end
