% Plots average values of delta h with burst against the burst size
clear variables
close all
warning off

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0 = 0.5;
Rcell = 0.2*a0;
p_initial = [0.5 0.65 0.8];
iniON = round(p_initial*N);

% parameters of the circuit
K = 16;
Son = 8;
Bvals = 0:100:1200;
nbins = 20;
delhvals = zeros(length(Bvals), length(p_initial));
delhstd = zeros(length(Bvals), length(p_initial));

% Supply list of folders in which to search, in the order of the p_initial
% values of the data
folderlist = {'2017-06-06 K16 Son8 deactivation', '2017-06-06 K16 Son8 autonomy', '2017-06-06 K16 Son8 activation'};

for j=1:length(folderlist)
    % Path to search for the saved data. It searchs by the name, defined by the
    % parameters chosen
    path = fullfile(pwd, 'figures', 'energy_change', 'burst', folderlist{j}); 
    straux = '(\d+)';
    %name = sprintf('N%d_n%d_a0%d_K_%d_Son_%d_bins_%d_B%d', ...
    %    N, iniON, 10*a0, K, Son, nbins, B);
    fpattern = sprintf('N%d_n%d_a0%d_K_%d_Son_%d_bins_%d_B%s',...
        N, iniON(j), 10*a0, K, Son, nbins, straux);
    
    % Get all file names in the directory
    listing = dir(path);
    num_files = numel(listing)-2; %first two entries are not useful
    count = 0;
    for i=1:num_files
        filename = listing(i+2).name;
        % remove extension and do not include txt files
        [~,name,ext] = fileparts(filename);
        if strcmp(ext, '.mat')
            count = count + 1;
            names{count} = name;
        end
    end

    for i=1:numel(names)
        [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
        if numel(tokens) > 0
            disp(names{i}) % displays the file name
            thisB = str2double(tokens{1,1}{1,1});
            idx = find(Bvals == thisB, 1); %find the right index for B
            % load the data
            load(fullfile(path,strcat(names{i},'.mat')), 'delh_av', 'delh_std_mean');
            delhvals(idx,j) = delh_av;
            delhstd(idx,j) = delh_std_mean;
        end
    end
end
%%
set(0, 'defaulttextinterpreter', 'latex');
h1=figure();
errorbar(repmat(Bvals, length(p_initial),1)', delhvals, delhstd, 'o--');
xlabel('$$B$$');
ylabel('$$\langle \Delta h \rangle$$');
set(gca, 'FontSize', 24);
xlim([-50 1250]);
legend('p=0.5', 'p=0.65', 'p=0.8');

%% Save figure
qsave=1;
if qsave
    fname = sprintf('N%d_a0%d_K_%d_Son_%d_bins_%d_B%d_to_%d', ...
        N, 10*a0, K, Son, nbins, Bvals(1), Bvals(end));
    out_fig = fullfile(pwd,'figures', 'energy_change', 'burst', ...
        strcat('Delh_vs_B_all_pini_', fname)); %filename
    save_figure_pdf(h1, 10, 8, out_fig);
end