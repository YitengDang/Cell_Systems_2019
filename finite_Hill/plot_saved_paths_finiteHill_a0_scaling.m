%% Study the spatial order of patterns in the finite Hill model
% Plot how mean final p, I change with a0 for fixed K, Con

% Parameters of the system
gridsize = 11;
N = gridsize^2;
a0_list = 5.2:0.05:5.7;
initialID = 'uniform';
K = 8;
Con = 16;
hill = 2;

% save?
qsave = 0;

%% Load data
% Path to search for the saved data. It searchs by the name, defined by the
% parameters chosen
path = 'H:\My Documents\Multicellular automaton\data\dynamics\finiteHill\2017-08-07'; 
straux = '(\d+)';
straux2 = '(\w+)';
%fpattern = strrep(sprintf('N%d_a0%d_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_%s-v%s',...
%            N, 10*a0, K, Con, hill, straux, straux2, initialID, straux), '.','p');
fpattern = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_%s-v%s',...
            N, a0, K, Con, hill, straux, straux2, initialID, straux), '.','p');

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

pend = zeros(numel(names),numel(a0_list));
Iend = pend;
a0_count = 0;
for a0=a0_list
    a0_count = a0_count + 1;
    fpattern = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_t%s_xmeanf_%s_%s-v%s',...
            N, a0, K, Con, hill, straux, straux2, initialID, straux), '.', 'p');
    for i = 1:numel(names)
        % first get the filename given by the student
        [tokens, ~] = regexp(names{i},fpattern,'tokens','match');
        if numel(tokens) > 0
            disp(names{i}) % displays the file name
            nruns = nruns + 1;
            % load the data
            load(fullfile(path,strcat(names{i},'.mat')), 'cells_hist', 'fN', 'I', 'mom', 't');
            % Get the sequence of p
            p = zeros(numel(cells_hist),1);
            for k = 1:numel(cells_hist)
                p(k) = sum(cells_hist{k})/N;
            end
            pend(i, a0_count) = p(end);
            Iend(i, a0_count) = I(end);
        end
    end
end
%% Parse data
% Only works if each a0 has an equal number of simualtions
nruns0 = 200; % number of sims per a0 value
idx = (pend ~= 0);
pend2 = reshape(pend(idx), [nruns0, numel(a0_list)]);
Iend2 = reshape(Iend(idx), [nruns0, numel(a0_list)]);

%% Plot final p and I against a0
set(0, 'defaulttextinterpreter', 'latex');
h1=figure();
plot(a0_list, mean(pend2, 1), 'o-', 'LineWidth', 1.5);
xlim([5.2 5.7]);
xlabel('$$a_0$$', 'FontSize', 24);
ylabel('$$\langle p_{final} \rangle$$', 'FontSize', 24);
set(gca,'FontSize', 20);

h2=figure();
plot(a0_list, mean(Iend2, 1), 'o-', 'LineWidth', 1.5);
xlim([5.2 5.7]);
xlabel('$$a_0$$', 'FontSize', 24);
ylabel('$$\langle I_{final} \rangle$$', 'FontSize', 24);
set(gca,'FontSize', 20);

qsave = 0;
if qsave
    name = strrep(sprintf('N%d_a0_%.2fto%.2f_K%d_Con%d_hill%.2f_%s_runs%d',...
                N, a0(1), a0(end), K, Con, hill, initialID, nruns0), '.','p');
    out_file = fullfile(pwd, 'figures',...
    'finite_Hill_final_pI_distribution', strcat(name,'_p_scaling'));
    out_file2 = fullfile(pwd, 'figures',...
        'finite_Hill_final_pI_distribution', strcat(name,'_I_scaling'));
    save_figure_pdf(h1, 10, 8, out_file);
    save_figure_pdf(h2, 10, 8, out_file2);
end