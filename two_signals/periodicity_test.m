function period_out = periodicity_test(cells_hist, save_consts_struct)

N = save_consts_struct.N;
gz = sqrt(N);
t_out = numel(cells_hist)-1;
%% Plot loaded trajectory
%{
hin = figure(1);
[pos,~,~] = init_cellpos_hex(gz,gz);
%a0 = 1;
disp_mol = 1;
cell_type = zeros(N,1);
cells = cells_hist{1};
update_cell_figure_continuum(hin, pos, cells, cell_type, 0, disp_mol);
for t=1:numel(cells_hist)-1
    waitforbuttonpress;
    %pause(1);
    cells = cells_hist{t+1};
    update_cell_figure_continuum(hin, pos, cells, cell_type, t, disp_mol);
    %update_cell_figure(hin, pos, a0, cells, cell_type, i-1);
end
%}
%% Periodicity test
% Scan over all initial frames
period = zeros(t_out+1, 1);
for t1=0:t_out-2
    %disp(t1);
    %t1 = t_out-10; % initial frame
    cells_ini = cells_hist{t1+1};
    for t=t1+1:t_out
        %disp(t);
        cells = cells_hist{t+1};
        if all(all(cells==cells_ini))
            period(t1+1) = t-t1;
            fprintf('t1=%d, period %d \n', t1, t-t1);
            break
        end
    end
end
disp('finished');
%% return variables
% return first found period
idx = find(period~=0, 1);
period_out = period(idx);
if isempty(period_out)
    period_out = 0;
end
%% Check found periods
%{
% check frame at t1, compare with frame at t2
t1 = 323; 
t2 = t1 + period(t1+1); % using found period
cells_1 = cells_hist{t1+1};
cells_2 = cells_hist{t2+1};
disp_mol = 1;

h1 = figure(1);
h2 = figure(2);
update_cell_figure_continuum(h1, pos, cells_1, cell_type, t1, disp_mol);
update_cell_figure_continuum(h2, pos, cells_2, cell_type, t2, disp_mol);
%}
end