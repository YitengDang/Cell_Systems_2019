% Plot all parameters of a saved dynamics data for all the files saved that
% have the set of parameters established, including noise
clear variables
close all
warning off

% Parameters of the system
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;
iniON = 48;
noise = 0.26;
radius = 0.25;
% parameters
Son = 11.4;
K = 5.2;

path = 'D:\eduardopavinat\Dropbox\Matlab codes\data_onecelltype_entropy\dynamics_noise\circle';
straux = '(\d+)';
fpattern = sprintf('N%d_n%d_neq_%s_a0%d_K_%d_Son_%d_r_%d_noise_%d_t_%s-v%s', ...
    N, iniON, straux, round(10*a0), round(K), round(Son), round(10*radius), ...
    round(10*noise), straux, straux);

% Get all file names
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
count = 0;
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        count = count + 1;
        names{count} = name;
    end
end

h1 = figure(1);
hold on
h2 = figure(2);
hold on
first = true;
for i = 1:numel(names)
    % first get the filename given by the student
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        disp(names{i})
        load(fullfile(path,strcat(names{i},'.mat')), 'cells_hist', 'fN', 'I', 'mom', 't');
        p = zeros(numel(cells_hist),1);
        for k = 1:numel(cells_hist)
            p(k) = sum(cells_hist{k})/N;
        end
        figure(h1)
        % plot energy contour
        if first
            first = false;
            E = @(p, I) -0.5*(Son-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
                -(2*p-1).*(0.5*(Son+1)*(1+fN) - K);
            piv = (0:N)/N;
            Iv = -0.05:0.05:1;
            [p_i,Imesh] = meshgrid(piv, Iv);
            contourf(piv, Iv, E(p_i, Imesh),'LineStyle', 'none')
            colormap('summer')
            c = colorbar;
        end
        plot(p, I, 'Color', 'r')
        scatter(p(end), I(end), 'kx')
        figure(h2)
        plot(0:t, -mom/N, 'Color', 'r')
        scatter(t, -mom(end)/N, 'kx')
    end
end

figure(h1)
hold off

set(gca,'FontSize', 20)
xlabel('p', 'FontSize', 24)
ylabel('I', 'FontSize', 24)
ylabel(c, 'h=H/N', 'FontSize', 24)
%ylim([-0.05 0.3])

name = sprintf('N%d_n%d_a0%d_K_%d_Son_%d_radius_%d_noise_%d', ...
    N, iniON, round(10*a0), round(K), round(Son), round(10*radius),round(10*noise));
out_file = fullfile(pwd, 'figures', 'dynamics_noise', 'circle', strcat(name,'_h_path'));
save_figure_pdf(h1, 10, 6, out_file);

figure(h2)
hold off

set(gca,'FontSize', 20)
xlabel('Time (steps)', 'FontSize', 24)
ylabel('h', 'FontSize', 24)
%ylim([-0.05 0.3])
box on
out_file = fullfile(pwd, 'figures', 'dynamics_noise','circle', strcat(name,'_allH'));
save_figure_pdf(h2, 10, 6, out_file);
