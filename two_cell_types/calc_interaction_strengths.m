function [f_mat, g_mat, f_std, g_std] = calc_interaction_strengths(N, frac_type1, a0, rcell, dist, Mcomm)
% calculates the interaction strengths f_lm by averaging over many lattices
%close all
%clear all
%warning off
set(0, 'defaulttextinterpreter', 'latex');
%--------------------------------------------------------------------------
% Set parameters of the system
%gridsize = 15;
%N = gridsize^2;

% Parameters 
%frac_type1 = 0.5;
N1 = round(frac_type1*N);
N2 = N - N1;
%a0 = 0.5;
%rcell = [0.2 0.2];
Rcell = rcell*a0;

ntrials = 10^3;
%Mcomm = [1 1; 1 1]; % Communication matrix, type i reacts to type j iff M_ij=1
%Mcomm = logical(Mcomm);
%--------------------------------------------------------------------------
% Distance and position
%[dist, pos] = init_dist_hex(gridsize, gridsize);

%% Check whether result already exists
skip = 0;
subfolder = sprintf('Mcomm_%d_%d_%d_%d', Mcomm(1,1), Mcomm(1,2), Mcomm(2,1), Mcomm(2,2));
fname_pattern = strrep(...
    sprintf('N1_%d_N2_%d_frac1_%.2f_a0_%.2f_Rcell1_%.2fa0_Rcell2_%.2fa0_ntrials_%s_%s-v%s',...
    N1, N2, frac_type1, a0, rcell(1), rcell(2), '(\d+)', subfolder, '(\d+)'), '.', 'p');
%fname_pattern = 'N1_113_N2_112_frac1_0p50_a0_0p50_Rcell1_0p20a0_Rcell2_0p20a0_ntrials_(\d+)-v(\d+)';
path = fullfile(pwd, 'data', 'interaction_strengths');
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
names = {};
for i=1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        names{end+1} = name;
        [tokens, ~] = regexp(names{i}, fname_pattern, 'tokens', 'match');
        if numel(tokens)>0 % load first version of found file
            load(fullfile(path, filename));
            skip = 1;
            break
        end
    end
end


%% calculate single interaction strengths
if ~skip
    f11 = zeros(ntrials, 1); 
    f12 = f11; f21 = f11; f22 = f11;
    g11 = f11; g12 = f11; g21 = f11; g22 = f11;
    for i=1:ntrials
        cell_type = zeros(N,1);
        idx1 = randperm(N,N1);
        idx2 = setdiff(1:N, idx1);
        cell_type(idx2) = 1;
        % extra randomization of cell positions by requiring that cell_type also
        % has a Moran's I value of around 0
        [cell_type, ~] = generate_I_new(cell_type, 0, 0.01, dist, a0);
        idx2 = find(cell_type);
        idx1 = setdiff(1:N, idx2);

        % Matrix of cell reading
        M = zeros(size(dist));
            % same type
        M(idx1, idx1) = Mcomm(1,1)*sinh(Rcell(1))./(a0*dist(idx1, idx1)).*exp(Rcell(1)-a0*dist(idx1, idx1));
        M(sub2ind(size(M), idx1, idx1)) = Mcomm(1,1)*1; % 1: normal, 0: no self communication

        M(idx2, idx2) = Mcomm(2,2)*sinh(Rcell(2))./(a0*dist(idx2, idx2)).*exp(Rcell(2)-a0*dist(idx2, idx2));
        M(sub2ind(size(M), idx2, idx2)) = Mcomm(2,2)*1; % 1: normal, 0: no self communication
            % different type
        M(idx1, idx2) = Mcomm(1,2)*Rcell(2)/(Rcell(1))*sinh(Rcell(1))./(a0*dist(idx1, idx2)).*exp(Rcell(2)-a0*dist(idx1, idx2)); % conc. cell of type 1 senses due to cells of type 2
        M(idx2, idx1) = Mcomm(2,1)*Rcell(1)/(Rcell(2))*sinh(Rcell(2))./(a0*dist(idx2, idx1)).*exp(Rcell(1)-a0*dist(idx2, idx1));
            % self-communication
        %M(sub2ind(size(M), 1:N, 1:N)) = 1; % 1: normal, 0: no self communication

        % Interaction strengths
        f11(i) = sum(sum(M(idx1, idx1)))/N;
        f12(i) = sum(sum(M(idx1, idx2)))/N; 
        f21(i) = sum(sum(M(idx2, idx1)))/N; %not equal to f12 if cell radii different
        f22(i) = sum(sum(M(idx2, idx2)))/N;

        g11(i) = sum(sum(M(idx1, idx1).^2))/N;
        g12(i) = sum(sum(M(idx1, idx2).^2))/N; 
        g21(i) = sum(sum(M(idx2, idx1).^2))/N; %not equal to f12 if cell radii different
        g22(i) = sum(sum(M(idx2, idx2).^2))/N;
    end

    qsave = 1;
    if qsave
        fname_str = strrep(...
            sprintf('N1_%d_N2_%d_frac1_%.2f_a0_%.2f_Rcell1_%.2fa0_Rcell2_%.2fa0_ntrials_%d_%s',...
            N1, N2, frac_type1, a0, rcell(1), rcell(2), ntrials, subfolder), '.', 'p');
        i = 1;
        fname = fullfile(pwd, 'data', 'interaction_strengths',...
            strcat(fname_str,'-v',int2str(i),'.mat'));
        while exist(fname, 'file') == 2
            i=i+1;
            fname = fullfile(pwd, 'data', 'interaction_strengths',...
                strcat(fname_str,'-v',int2str(i),'.mat'));
        end

        % save data
        save(fname, 'N1', 'N2', 'Rcell', 'a0', 'ntrials', 'f11', 'f12', 'f21', 'f22',...
            'g11', 'g12', 'g21', 'g22');
    end
end

%% Package end result
f_mat = [mean(f11) mean(f12); mean(f21) mean(f22)];
g_mat = [mean(g11) mean(g12); mean(g21) mean(g22)];
f_std = [std(f11) std(f12); std(f21) std(f22)];
g_std = [std(g11) std(g12); std(g21) std(g22)];

%%
% Plot average and std of interaction strengths
% Confirmed that std is always small
%figure();
%hold on
%c = categorical({'f11', 'f12', 'f21', 'f22'});
%f_avg = [mean(f11) mean(f12) mean(f21) mean(f22)];
%f_std = [std(f11) std(f12) std(f21) std(f22)];
%bar(c, f_mat(:))
%errorbar(1:4, f_mat(:), f_std(:), '.');
%}

%}
