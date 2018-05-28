%% Makes a movie of a saved trajectory
close all
clear all
warning off

%% Load trajectory
% Parameters
%gridsize = 11;
%N = gridsize^2;
%a0 = 3.5;
%K = 7;
%Con = 13;
%hill = 2; % Hill coefficient
%initialID = 'uniform'; %'uniform' 'allON' 'singleOFF' 'fixedp'
% If fixedp, assume initial states ~ Normal(p,sigma)
%p = 0.4;
%sigma = 0.3; 
qsave = 0;
version = 1;
%tf = 28;
%xmf = 0.65;
%fname_str = strrep(sprintf('N%d_a0_%.2f_K%d_Con%d_hill%.2f_t%d_xmeanf_%.2f_%s-v%d',...
%    N, a0, K, Con, hill, tf, xmf, initialID, version), '.','p');
%fname_str = 'N121_a0_5p20_K8_Con16_hill2p00_t293_xmeanf_0p66_uniform-v1';
%fname = fullfile(pwd, 'data', 'dynamics', 'finiteHill', '2017-08-07', fname_str);
fname_str = 'N121_n61_neq_20_a015_K_17_Son_14_t_7-v1';
fname = fullfile(pwd, 'temp', fname_str);
load(fname);

%% video for saving
frames = struct('cdata',[],'colormap',[]);
hin = figure();
for i=1:t+1
    cells = cells_hist{i};
    %update_cell_figure_continuum(hin, pos, a0, cells, cell_type, i);
    update_cell_figure(hin, pos, a0, cells, cell_type, i);
    frames(i) = getframe(gcf);
end

%% Check video
%
figure(6);
fps = 1;
movie(frames, 1, fps)

%% Save video
i = 1;
fname = fullfile(pwd, 'videos', ...%'finiteHill',...
    strcat(fname_str,'-v',int2str(i),'.avi'));
while exist(fname, 'file') == 2
    i=i+1;
    fname = fullfile(pwd, 'videos', ...%'finiteHill',...
        strcat(fname_str,'-v',int2str(i),'.avi'));
end

myVideo = VideoWriter(fname, 'Uncompressed AVI');
myVideo.FrameRate = fps;  % Default 30
open(myVideo);
writeVideo(myVideo, frames);
close(myVideo);
%}