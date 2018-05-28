% Calculates the energy changes for simulated runs and plots them in a
% histogram
clear variables
close all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
p_initial = 0.6;
iniON = round(p_initial*N);
% parameters of the circuit
Son = 8;
K = 16;
B = 1000;

% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = 'H:\My Documents\Multicellular automaton\data\dynamics\burst';
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

% Extract starting and final (p,I) values and compute hamiltonian change
pini = zeros(numel(names),1);
pend = pini;
Iini = pini;
Iend = pini;
hini = pini;
hend = pini;
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
        hini(count) = mom(1);
        hend(count) = mom(end);
        tend(count) = t;
        count = count+1;
    else %delete entries that do not match parameters
        pend(count) = [];
        Iend(count) = [];
        pini(count) = [];
        Iini(count) = [];
        hini(count) = [];
        hend(count) = [];
        tend(count) = [];
    end
end

% CAUTION: mom does not always give the right energy. Calculate directly
% from p, I.
E = @(p, I) -0.5*(Son-1)*(1 + 4*fN.*p.*(1-p).*I + fN*(2*p-1).^2) ...
            -(2*p-1).*(0.5*(Son+1)*(1+fN) - K);
delh = E(pend,Iend)-E(pini,Iini);

% Check whether energy from (p,I) is the same as mom
%{
disp(delh==hend-hini);
figure(2);
scatter(delh, hend-hini)
%}

% Plot histogram
h1 = figure(1);
nbins=20;
histogram(delh/N,nbins,'Normalization','count');
set(gca,'FontSize', 20)
xlabel('$$\Delta h$$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Count', 'Interpreter', 'latex', 'FontSize', 24);
box on

%%
% save the histogram
name = sprintf('N%d_n%d_a0%d_K_%d_Son_%d_bins_%d_B%d', ...
    N, iniON, 10*a0, K, Son, nbins, B);
out_file = fullfile(pwd,'figures', 'energy_change', 'burst', ...
    strcat(name,'_h_hist')); %filename
save_figure_pdf(h1, 8, 6, out_file);

% save file variables
%fname_str = sprintf('N%d_n%d_a0%d_K_%d_Son_%d_B%d', ...
%    N, iniON, 10*a0, K, Son);
fname = fullfile(pwd, 'figures', 'energy_change', 'burst', ...
        strcat(name,'.mat'));  
save(fname)