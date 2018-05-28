% Calculates the energy changes for simulated runs and plots them in a
% histogram

% CAUTION: mom does not always give the right energy. Calculate directly
% from p, I.
%
for B=0:100:1200
%clearvars
clearvars -except B
close all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
p_initial = 0.65;
iniON = round(p_initial*N);
% parameters of the circuit
K = 16;
Son = 8;
%B = 1200;
qsave = 1;

% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = fullfile(pwd, 'data','dynamics', 'burst', '2017-06-06 autonomy'); 
straux = '(\d+)';
fpattern = sprintf('N%d_n%d_neq%s_a0%d_K%d_Son%d_t%s_B%d-v%s', ...
    N, iniON, straux, 10*a0, K, Son, straux, B, straux);

% Get all file names in the directory
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
%}
%%
% Extract starting and final (p,I) values and compute hamiltonian change
pini = zeros(numel(names),1);
%p2 = pini; %p(t=2) 
pend = pini;
Iini = pini;
%I2 = pini; %I(t=2)
Iend = pini;
tend = pini;
count = 1;
for i = 1:numel(names)
    % first get the filename given by the student
    [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
    if numel(tokens) > 0
        disp(names{i}) % displays the file name
        % load the data
        load(fullfile(path,strcat(names{i},'.mat')), 'Non', 'fN', 'I', 'mom', 't');
        % Get the sequence of p
        p = Non/N;
        % Save the initial and end points for plotting the marks later
        pend(count) = p(end);
        Iend(count) = I(end);
        pini(count) = p(1);
        Iini(count) = I(1);
        %p2(count) = p(2);
        %I2(count) = I(2);
        tend(count) = t;
        count = count+1;
    else %delete entries that do not match parameters
        pend(count) = [];
        Iend(count) = [];
        pini(count) = [];
        Iini(count) = [];
        %p2(count) = [];
        %I2(count) = [];
        tend(count) = [];
    end
end

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength
E = @(p, I) -0.5*(Son-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
            -(2*p-1).*(0.5*(Son+1)*(1+fN) - K);
delH = E(pend,Iend)-E(pini,Iini);
%W = E(p2, I2)-E(pini,Iini);

% Check whether energy from (p,I) is the same as mom
%{
disp(delh==hend-hini);
figure(2);
scatter(delh, hend-hini)
%}
%%
% Plot delH histogram
h1 = figure();
nbins=20;
histogram(delH/N,nbins,'Normalization','count');
set(gca,'FontSize', 24)
xlabel('$$\Delta h$$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Count', 'Interpreter', 'latex', 'FontSize', 24);
box on

%{
% Plot W histogram
h2 = figure();
nbins=20;
histogram(W/N,nbins,'Normalization','count');
set(gca,'FontSize', 24);
xlabel('$$W$$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Count', 'Interpreter', 'latex', 'FontSize', 24);
box on

% Compare delH and W in scatter plot
Hlim = linspace(min(delH/N), max(delH/N), 100);
h3 = figure();
hold on
scatter(delH/N, W/N);
plot(Hlim, Hlim);
set(gca,'FontSize', 24);
xlabel('$$\Delta H/N$$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$$W/N$$', 'Interpreter', 'latex', 'FontSize', 24);
hold off
%}

% Store mean and std energy change 
delh_av = mean(delH/N);
delh_std = std(delH/N);
delh_std_mean = std(delH/N)/sqrt(count-1); %number of runs = count-1
%W_av = mean(W/N);
%%
% save the histograms

if qsave
    name = sprintf('N%d_n%d_a0%d_K_%d_Son_%d_bins_%d_B%d', ...
        N, iniON, 10*a0, K, Son, nbins, B);
    out_fig = fullfile(pwd,'figures', 'energy_change', 'burst', ...
        strcat(name,'_delh_hist')); %filename
    save_figure_pdf(h1, 10, 6, out_fig);
    %out_fig2 = fullfile(pwd,'figures', 'hamiltonian_map', 'burst energy change', ...
    %    strcat(name,'_W_hist')); %filename
    %save_figure_pdf(h2, 10, 6, out_fig2);
    % save scatter plot
    %out_fig3 = fullfile(pwd,'figures', 'hamiltonian_map', 'burst energy change', ...
    %    strcat(name,'_delh_W_scatter')); %filename
    %save_figure_pdf(h3, 10, 6, out_fig3);
    % save file variables
    close all;
    out_name = fullfile(pwd, 'figures', 'energy_change', 'burst', ...
            strcat(name,'.mat'));  
    save(out_name);
end

end