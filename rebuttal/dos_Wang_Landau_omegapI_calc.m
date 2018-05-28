function [log_omegapI, pv, Iv] = dos_Wang_Landau_omegapI_calc(N, a0, path)
%% Load files
%L = 11; % size of lattice
%N = L^2;
%a0 = 0.5;
n = 11:61;
nbins = 30; % number of histogram bins
c1=1; f0 = exp(c1); % initial modification factor 
c2=4; ffinal = exp(10^(-c2)); % final modification factor
p_flat=0.8; % flatness parameter

g_all = zeros(nbins, numel(n));
for i=1:numel(n)
    this_n = n(i);
    fileid = strrep(sprintf('WL_norm_N%d_n0_%d_a0_%.2f_f0exp%d_ffin_e10e-%d_pflat%.1f_%dbins_fixed%.2f',...
        N, this_n, a0, c1, c2, p_flat, nbins), '.', 'p');
    fname = fullfile(path, strcat(fileid, '.mat') ); % filename
    load(fname, 'lgn2', 'edges');
    g_all(:, i) = exp(lgn2)/sum(exp(lgn2));
end

% Get binomial coefficients (possibly not exact)
fileid = sprintf('nchoosek_n%d', N);
fname = fullfile(pwd,'data','nchoosek', strcat(fileid, '.mat'));
if exist(fname)==2
    %disp('exists!');
    load(fname);
else
    %disp('Doesnt exist!');
    binom = zeros(1,N+1);
    for j=1:N+1  
        disp(j);
        % if not, calculate
        binom(j) = nchoosek(N,j-1);
    end
    save(fname, 'binom');
end

omegap = binom(n);

% calculate omegapI
pv = unique(sort([n N-n]))/N; % remove duplicates
Iv = (edges(1:end-1)+edges(2:end))/2;
omegapI_half = g_all.*repmat(omegap, nbins, 1); % up to values of p=0.5
log_omegapI = log([omegapI_half fliplr(omegapI_half(:, 1:end-mod(N,2)-1 ) )]);

end