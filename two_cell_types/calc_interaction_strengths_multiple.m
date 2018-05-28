% calculates the interaction strengths f_lm by averaging over many lattices
%close all
%clear all
%warning off
set(0, 'defaulttextinterpreter', 'latex');
%--------------------------------------------------------------------------
% Set parameters of the system
gridsize = 15;
N = gridsize^2;

% Parameters 
frac_type1 = 0.5;
N1 = round(frac_type1*N);
N2 = N - N1;
a0 = 0.5;
rcell = [0.2 0.2];
Rcell = rcell*a0;

Mcomm = [1 1; 1 1]; % Communication matrix, type i reacts to type j iff M_ij=1
Mcomm = logical(Mcomm);
%--------------------------------------------------------------------------
% Distance and position
[dist, pos] = init_dist_hex(gridsize, gridsize);

%% Plot for a range of a0 values
frac_type1_all = 0.5; %0.0:0.1:0.5;
%for var=1:numel(frac_type1_all)
I_all = -0.1:0.05:0.5;
a0 = 1;
%for var = 1:numel(I_all)
    frac_type1 = 0.5; 
    %frac_type1_all(var);
    N1 = round(frac_type1*N);
    N2 = N - N1;

    a0_range = 0.1:0.1:2;
    ntrials = 10^2;
    %f11 = zeros(numel(a0_range), ntrials); 
    f11 = zeros(numel(I_all), ntrials); 
    f12 = f11; f21 = f11; f22 = f11;
    g11 = f11; g12 = f11; g21 = f11; g22 = f11;

    %for i=1:numel(a0_range)
    for i=1:numel(I_all)   
        
        I_12 = I_all(i);
        fprintf('I12 = %.1f', I_12);
        %a0 = a0_range(i);
        rcell = [0.2 0.2];
        Rcell = rcell*a0;
        for j=1:ntrials
            disp(j);
            cell_type = zeros(N,1);
            idx1 = randperm(N,N1);
            idx2 = setdiff(1:N, idx1);
            cell_type(idx2) = 1;
            
            % --No extra randomization here!--
            % generate lattices with spatial order in distribution of cells
            % of type 1, 2 (fix I_12).
            [cell_type, ~] = generate_I_new(cell_type, I_12, I_12+0.01, dist, a0);
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
            M1 = M; 
            M1(sub2ind(size(M1), 1:N, 1:N)) = 0; % remove self-interaction terms
            f11(i,j) = sum(sum(M1(idx1, idx1)))/N;
            f12(i,j) = sum(sum(M1(idx1, idx2)))/N; 
            f21(i,j) = sum(sum(M1(idx2, idx1)))/N; %not equal to f12 if cell radii different
            f22(i,j) = sum(sum(M1(idx2, idx2)))/N;

            % Interaction strengths
            g11(i,j) = sum(sum(M1(idx1, idx1).^2))/N;
            g12(i,j) = sum(sum(M1(idx1, idx2).^2))/N; 
            g21(i,j) = sum(sum(M1(idx2, idx1).^2))/N; %not equal to f12 if cell radii different
            g22(i,j) = sum(sum(M1(idx2, idx2).^2))/N;
        end
    end
    %% Plot against a0
    h1=figure(1);
    errorbar(I_all, mean(f11, 2), std(f11, 1, 2), '-', 'LineWidth', 1.5);
    hold on
    errorbar(I_all, mean(f12, 2), std(f12, 1, 2), '-', 'LineWidth', 1.5);
    errorbar(I_all, mean(f21, 2), std(f21, 1, 2),'-', 'LineWidth', 1.5);
    errorbar(I_all, mean(f22, 2), std(f22, 1, 2),'-', 'LineWidth', 1.5);
    legend({'f11', 'f12', 'f21', 'f22'});
    %xlabel('$$a_0$$');
    xlabel('$I_{1,2}$');
    ylabel('$$f_{lm}$$');
    set(gca,'FontSize', 24);
    %%
    h2=figure(2);
    errorbar(I_all, mean(g11, 2), std(g11, 1, 2), '-', 'LineWidth', 1.5);
    hold on
    errorbar(I_all, mean(g12, 2), std(g12, 1, 2), '-', 'LineWidth', 1.5);
    errorbar(I_all, mean(g21, 2), std(g21, 1, 2),'-', 'LineWidth', 1.5);
    errorbar(I_all, mean(g22, 2), std(g22, 1, 2),'-', 'LineWidth', 1.5);
    legend({'g11', 'g12', 'g21', 'g22'});
    %xlabel('$$a_0$$');
    xlabel('$I_{1,2}$');
    ylabel('$$g_{lm}$$');
    set(gca,'FontSize', 24);
    %% Save data and fig
    f11_avg = mean(f11, 2); f11_std = std(f11, 1, 2);
    f12_avg = mean(f12, 2); f12_std = std(f12, 1, 2);
    f21_avg = mean(f21, 2); f21_std = std(f21, 1, 2);
    f22_avg = mean(f22, 2); f22_std = std(f22, 1, 2);
    g11_avg = mean(g11, 2); g11_std = std(g11, 1, 2);
    g12_avg = mean(g12, 2); g12_std = std(g12, 1, 2);
    g21_avg = mean(g21, 2); g21_std = std(g21, 1, 2);
    g22_avg = mean(g22, 2); g22_std = std(g22, 1, 2);
    
    qsave = 0;
    if qsave
        fname_str = strrep(...
            sprintf('scaling_I12_N1_%d_N2_%d_frac1_%.2f_a0_%.2fto%.2f_Rcell1_%.2fa0_Rcell2_%.2fa0_ntrials_%d',...
            N1, N2, frac_type1, a0_range(1), a0_range(end), rcell(1), rcell(2), ntrials), '.', 'p');
        i = 1;
        fname = fullfile(pwd, 'data', 'interaction_strengths',...
            strcat(fname_str,'-v',int2str(i),'.mat'));
        while exist(fname, 'file') == 2
            i=i+1;
            fname = fullfile(pwd, 'data', 'interaction_strengths',...
                strcat(fname_str,'-v',int2str(i),'.mat'));
        end

        % save figure
        fname_str = strrep(...
            sprintf('scaling_f_lm_I12_N1_%d_N2_%d_frac1_%.2f_a0_%.2fto%.2f_Rcell1_%.2fa0_Rcell2_%.2fa0_ntrials_%d',...
            N1, N2, frac_type1, a0_range(1), a0_range(end), rcell(1), rcell(2), ntrials), '.', 'p');
        fname1 = fullfile(pwd, 'figures', 'interaction_strengths',...
            strcat(fname_str,'-v',int2str(i)));
        save_figure_pdf(h1, 10, 8, fname1)
        
        fname_str = strrep(...
            sprintf('scaling_g_lm_I12_N1_%d_N2_%d_frac1_%.2f_a0_%.2fto%.2f_Rcell1_%.2fa0_Rcell2_%.2fa0_ntrials_%d',...
            N1, N2, frac_type1, a0_range(1), a0_range(end), rcell(1), rcell(2), ntrials), '.', 'p');       
        fname2 = fullfile(pwd, 'figures', 'interaction_strengths',...
            strcat('g_lm_', fname_str,'-v',int2str(i)));
        save_figure_pdf(h2, 10, 8, fname2)
        close all;
        % save data
        save(fname, 'a0_range', 'f11_avg', 'f12_avg', 'f21_avg', 'f22_avg',...
            'f11_std', 'f12_std', 'f21_std', 'f22_std', ...
            'g11_avg', 'g12_avg', 'g21_avg', 'g22_avg',...
            'g11_std', 'g12_std', 'g21_std', 'g22_std');
    end
%end

%% Core code
%{
%--------------------------------------------------------------------------
% Set parameters of the system
gridsize = 15;
N = gridsize^2;
frac_type1 = 0.5;
N1 = round(frac_type1*N);
N2 = N - N1;
a0 = 1.5;

% Parameters 
%Con = [35 5]; %column j = type j 
%K = [5 5];
%p = [0.5 0.5];
%Rcell = [0.2 0.2]*a0;

Mcomm = [1 1; 1 1]; % Communication matrix, type i reacts to type j iff M_ij=1
Mcomm = logical(Mcomm);
%--------------------------------------------------------------------------
% Distance and position
[dist, pos] = init_dist_hex(gridsize, gridsize);

rcell = [0.2 0.2];
Rcell = rcell*a0;

% cell types
cell_type = zeros(N,1);
idx1 = randperm(N,N1);
idx2 = setdiff(1:N, idx1);
cell_type(idx2) = 1;

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
M1 = M; 
M1(sub2ind(size(M1), 1:N, 1:N)) = 0; % remove self-interaction terms
f11 = sum(sum(M1(idx1, idx1)))/N1;
f12 = sum(sum(M1(idx1, idx2)))/N1; 
f21 = sum(sum(M1(idx2, idx1)))/N2; %not equal to f12 if cell radii different
f22 = sum(sum(M1(idx2, idx2)))/N2;

% Interaction strengths
g11 = sum(sum(M1(idx1, idx1).^2))/N1;
g12 = sum(sum(M1(idx1, idx2).^2))/N1; 
g21 = sum(sum(M1(idx2, idx1).^2))/N2; %not equal to f12 if cell radii different
g22 = sum(sum(M1(idx2, idx2).^2))/N2;
%}