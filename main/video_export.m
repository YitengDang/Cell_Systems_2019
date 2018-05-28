% makes and exports video from saved frames
clear variables
close all

%% Give parameters for generating filename
% fixed lattice parameters
gridsize = 15;
N = gridsize^2;
a0 = 1.5;
Rcell = 0.2*a0;

% use hexagonal lattice
[dist, pos] = init_dist_hex(gridsize, gridsize);
dist_vec = a0*dist(1,:);
r = dist_vec(dist_vec>0); % exclude self influence
fN = sum(sinh(Rcell)*sum(exp(Rcell-r)./r)); % calculate signaling strength

% fixed circuit parameters
Son = 18;
K = 14;
hill = 2;

% initial conditions
p0 = 0.27;
iniON = round(p0*N);
sigma = 0.1;

% variable lattice parameters
%{
Son_0 = 5; 
K_0 = (Son_0+1)/2*(1+fN);

tsteps = 10; % complete procedure in tsteps steps
%Son_all=[5 5.24528812001061 5.53124349331728 5.86746718837490 6.26660291113958 6.74553754586907 7.32716227705940 8.04299585757536 8.93715627322097 10.0724797646909 11.5401266664718];
Son_all = linspace(5, 15, tsteps+1);
K_all = linspace(K_0, K_0, tsteps+1);
B_all = (Son_all+1)/2*(1+fN) - K_all;
%}
%% Import data
straux = '(\d+)';
%{
fpattern = strrep(sprintf('N%d_n%d_a0_%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d_frame%s-v%s', ...
            N, iniON, a0, Son_all(1), Son_all(end), B_all(1),...
            round(B_all(end),2), tsteps, straux, straux), '.', 'p');
fpattern = strrep(sprintf('N%d_inirand_a0_%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d_frame%s-v%s', ...
            N, a0, Son_all(1), Son_all(end), B_all(1),...
            round(B_all(end),2), tsteps, straux, straux), '.', 'p');
%}
fpattern = strrep(sprintf('N%d_a0%.1f_Con%.2f_K%.2f_hill%.2f_%s_pini%.2f_sigma%.2f-v%s',...
    N, a0, Con, K, hill, initialID, p, sigma, straux), '.','p');

path = fullfile(pwd, 'figures', 'time_evolution');
listing = dir(path);
num_files = numel(listing)-2; %first two entries are not useful
names = {};
for i = 1:num_files
    filename = listing(i+2).name;
    % remove extension and do not include txt files
    [~,name,ext] = fileparts(filename);
    if strcmp(ext, '.fig')
        names{end+1} = name;
    end
end

frames = struct('cdata',[],'colormap',[]);
for i=1:numel(names)
    [tokens, ~] = regexp(names{i}, fpattern, 'tokens','match');
    if numel(tokens)>0
        disp(names{i});
        h=openfig(fullfile(path, strcat(names{i}, '.fig'))) ;
        %disp(str2double(tokens{1}{1}));
        frames(str2double(tokens{1}{1})) = getframe(gcf);
        close(h);
    end
end
%% display MATLAB movie (ignore axes that are off)
movie(frames, 1, 2)

%% export video
%myVideo = VideoWriter('myfile.avi');
out_file = strrep(sprintf('N%d_n%d_a0_%.1f_Son_%.fto%.f_B_%.2fto%.2f_tsteps%d', ...
            N, iniON, a0, Son_all(1), Son_all(end), B_all(1),...
            round(B_all(end),2), tsteps), '.', 'p');   
myVideo = VideoWriter(fullfile(pwd, 'figures', 'time_evolution', strcat(out_file, '.avi')), 'Uncompressed AVI');
myVideo.FrameRate = 1;  % Default 30
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);