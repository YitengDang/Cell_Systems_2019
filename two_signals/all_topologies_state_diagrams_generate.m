%% Generate and save state diagrams of all possible topologies (up to equivalence) 
clear variables
close all
clc
%% Estimate # of topologies
% # phases per interaction
n_phases_all = [6 2];

% Single cell
%{
disp('Single cell:');
n_phases = n_phases_all(2); % number of phases for each interaction
Ts = [1 0 3 0 10];
Ta = n_phases.^(0:4);
ns = [1 0 4 0 4]; % self-similar
na = [0 8 20 32 12]./2; % remaining

N_phases = sum(Ts.*ns+Ta.*na); % total # phases
N_phases_all = sum( (ns+2*na).*Ta );
fprintf('Total # single-cell phases (before symmetries): %d \n', N_phases_all)
fprintf('Total # single-cell phases (after symmetries): %d \n', N_phases)

% multiple cells
disp('Multiple cells:');
n_phases = n_phases_all(1);
Ta = n_phases.^(0:4);
n = zeros(1,5);
for k=0:4
    n(k+1) = nchoosek(4,k)*2^k;
end
Ns1 = sum(n.*Ta);

% account for symmetry
ns2 = [1 0 4 0 4]; % self-similar
na2 = [0 8 20 32 12]./2; % remaining
Ts2 = [1 0 21 0 666]; % #states for self-similar topologies
Ta2 = n_phases.^(0:4);
%Ns2a = sum(Ta.*na); % states without self-similarity
Ns2 = sum(Ts2.*ns2+Ta2.*na2);

% exclude trivial topologies
ns3 = [0 0 2 0 4];
na3 = [0 0 10 32 12]./2; % remaining
Ns3 = sum(Ts2.*ns3+Ta2.*na3);

% account for symmetries
fprintf('Total # multi-cell phases (before symmetries): %d \n', Ns1)
fprintf('Total # multi-cell phases (after symmetries): %d \n', Ns2)
%fprintf('States without self-similarity: %d \n', Ns2a)
fprintf('Excluding trivial topologies: %d \n', Ns3)

% Convention:
% phase 1: (Con - K) > 0
% phase 2: (Con - K) < 0
%}
%% Loop over topologies and phases
% N.B. some of the phases are not possible, e.g. when autonomy and
% activation-deactivation are both present for the same gene

% Settings
single_cell = 0;
draw_diagram = 0; % draw state diagram?
sym = 0; % include symmetries? 0: symmetric diagrams are excluded, 1: everything included
n_phases = n_phases_all(single_cell+1);

% tracking variables
done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
M = [0 1 -1]; % index to interaction
count = 0; %count topologies
countP = zeros(3^4, 1); % count phases

% main output
M_int_all = {}; %cell(3^4, 1);
M_int_by_topology = {}; 
phases_all = {}; % store all phases, stored as matrix with phase of each interaction (1-6), 0 if no interaction (M_int(i,j)==0)
phases_all_by_topology = {};
state_diagrams = {}; % cell(Ns2, 1); % store all state diagrams (graph transition matrices)
state_diagrams_by_topology = {};
steady_states = {}; %cell(Ns2, 1); % store all steady states
steady_states_by_topology = {};
cycles_all = {}; %cell(Ns2, 1); % store all loop structures
cycles_all_by_topology = {};

for k=1:3^4
    [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], k);
    % matrix associated with indices
    M_int = [M(i11) M(i12); M(i21) M(i22)];
    
    if all(M_int(:)==0) % skip the trivial network
        continue
    elseif done(i11,i12,i21,i22)
        continue
    %elseif (i11==1 && i12==1) || (i21==1 && i22==1) || (i12==1 && i21==1)
        % skip: no input to gene 1 || 2 || no cross-talk
        % between gene 1 and 2
        
        %disp(M_int);
        
        %done(i11,i12,i21,i22)=1;
        % also symmetry partner is done:
        %gM = g([i11 i12; i21 i22]);
        %done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
    else
        %disp(M_int);
        
        % update loop variables
        done(i11,i12,i21,i22)=1;
        count = count+1;
        
        % prepare variables
        M_int_by_topology{count} = M_int;
        phases_all_by_topology{count} = {};
        state_diagrams_by_topology{count} = {};
        steady_states_by_topology{count} = {};
        cycles_all_by_topology{count} = {};

        % get symmetric partners
        gM = g([i11 i12; i21 i22]);
        if ~sym
            done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
        end
        
        % loop over all phases
        ni = sum(abs(M_int(:))==1);
        sz = n_phases*abs(M_int)+(1-abs(M_int)); % size matrix, for correct index conversion
        doneP = zeros(sz(1,1), sz(1,2), sz(2,1), sz(2,2));
        
        % distinguish between self-similar topologies and others
        if all(all(g(M_int) == M_int)) % self-similar topology
            % topology symmetry 1<->2
            %disp(M_int);
            for k1=1:n_phases^ni
                [i11b, i12b, i21b, i22b] = ind2sub([sz(1,1), sz(1,2), sz(2,1), sz(2,2)], k1);
                if doneP(i11b,i12b,i21b,i22b)
                    continue
                elseif (i11b==1)&&(i21b==6)&&(M_int(1,1)~=0)&&(M_int(2,1)~=0) ||...
                        (i11b==6)&&(i21b==1)&&(M_int(1,1)~=0)&&(M_int(2,1)~=0) ||...
                        (i21b==1)&&(i22b==6)&&(M_int(2,1)~=0)&&(M_int(2,2)~=0) ||...
                        (i21b==1)&&(i22b==6)&&(M_int(2,1)~=0)&&(M_int(2,2)~=0)
                    % incompatible phases
                    doneP(i11b,i12b,i21b,i22b) = 1;
                    continue
                else
                    %disp(sum(countP));
                    P = [i11b i12b; i21b i22b]; % matrix with phases
                    
                    % calculate result
                    [A, ss, cycles] = all_topologies_analyze(single_cell, P, M_int, draw_diagram);
                    
                    % store results
                    phases_all{end+1} = abs(M_int).*P; % store phase matrix 
                    M_int_all{end+1} = M_int; %store interaction matrix
                    state_diagrams{end+1} = A;
                    steady_states{end+1} = ss;
                    cycles_all{end+1} = cycles;
                    
                    phases_all_by_topology{count}{end+1} = abs(M_int).*P; % store phase matrix 
                    state_diagrams_by_topology{count}{end+1} = A;
                    steady_states_by_topology{count}{end+1} = ss;
                    cycles_all_by_topology{count}{end+1} = cycles;
                    %}
                    
                    % update tracking variables: also consider P symmetries
                    doneP(i11b,i12b,i21b,i22b) = 1;
                    gP = g(P);
                    if ~sym
                        doneP(gP(1,1),gP(1,2),gP(2,1),gP(2,2))=1;
                    end
                    countP(k) = countP(k) + 1;
                    
                    if draw_diagram
                        pause(1);
                        close all
                    end
                end
            end
        else
            % topology asymmetric
            for k1=1:n_phases^ni
                [i11b, i12b, i21b, i22b] = ind2sub([sz(1,1), sz(1,2), sz(2,1), sz(2,2)], k1);
                if doneP(i11b,i12b,i21b,i22b)
                    continue
                elseif (i11b==1)&&(i21b==6)&&(M_int(1,1)~=0)&&(M_int(2,1)~=0) ||...
                        (i11b==6)&&(i21b==1)&&(M_int(1,1)~=0)&&(M_int(2,1)~=0) ||...
                        (i21b==1)&&(i22b==6)&&(M_int(2,1)~=0)&&(M_int(2,2)~=0) ||...
                        (i21b==1)&&(i22b==6)&&(M_int(2,1)~=0)&&(M_int(2,2)~=0)
                    % incompatible phases
                    doneP(i11b,i12b,i21b,i22b) = 1;
                    continue
                else
                    %
                    disp(sum(countP));
                    P = [i11b i12b; i21b i22b]; % matrix with phases
                    
                    % calculate result
                    [A, ss, cycles] = all_topologies_analyze(single_cell, P, M_int, draw_diagram);
                    
                    % store results
                    phases_all{end+1} = abs(M_int).*P; % store phase matrix 
                    M_int_all{end+1} = M_int; %store interaction matrix
                    state_diagrams{end+1} = A;
                    steady_states{end+1} = ss;
                    cycles_all{end+1} = cycles;
                    
                    phases_all_by_topology{count}{end+1} = abs(M_int).*P; % store phase matrix 
                    state_diagrams_by_topology{count}{end+1} = A;
                    steady_states_by_topology{count}{end+1} = ss;
                    cycles_all_by_topology{count}{end+1} = cycles;
                    
                    % update tracking variables: no need to consider P symmetries
                    doneP(i11b,i12b,i21b,i22b) = 1;
                    countP(k) = countP(k) + 1;
                    
                    if draw_diagram
                        pause(1);
                        close all
                    end
                end
            end
        end
        %}
        
    end
end

fprintf('Total # phases considered: %d \n', sum(countP))
%%
% save data
save_path = 'H:\My Documents\Multicellular automaton\data\two_signals\all_topologies';
%save_path = 'D:\Multicellularity\data\two_signals\all_topologies';
labels = {'multi_cell', 'single_cell', 'multi_cell_all_incl', 'single_cell_all_incl'};
qsave = 1;
if qsave
    label = labels{single_cell+1 + 2*sym};
    fname_str = sprintf('all_topologies_data_%s', label);
    fname = fullfile(save_path, fname_str);
    save(fname, 'single_cell', 'M_int_all', 'phases_all', 'state_diagrams', 'steady_states', 'cycles_all',...
        'M_int_by_topology', 'phases_all_by_topology', 'state_diagrams_by_topology',...
        'steady_states_by_topology', 'cycles_all_by_topology');
end

%% For given topology (M_int), loop over all phases
%{
countP = 0;
M_int = [0 -1; 1 1];

% Settings
single_cell = 1;
draw_diagram = 1; % draw state diagram?
save_diagram = 1;
labels = {'multi_cell', 'single_cell'};
label = labels{single_cell+1};

% loop over all phases
n_phases = n_phases_all(single_cell+1);
ni = sum(abs(M_int(:))==1);
sz = n_phases*abs(M_int)+(1-abs(M_int));
doneP = zeros(sz(1,1), sz(1,2), sz(2,1), sz(2,2));

phases_all = {}; % store all phases, stored as matrix with phase of each interaction (1-6), 0 if no interaction (M_int(i,j)==0)
state_diagrams = {}; % cell(Ns2, 1); % store all state diagrams (graph transition matrices)
steady_states = {}; %cell(Ns2, 1); % store all steady states
cycles_all = {}; %cell(Ns2, 1); % store all loop structures

count = 0;
for k1=1:n_phases^ni
    disp(k1)
    [i11b, i12b, i21b, i22b] = ind2sub([sz(1,1), sz(1,2), sz(2,1), sz(2,2)], k1);
    if doneP(i11b,i12b,i21b,i22b)
        continue
    elseif (i11b==1)&&(i12b==6) || (i11b==6)&&(i12b==1) 
        doneP(i11b,i12b,i21b,i22b) = 1;
        continue
    elseif (i21b==1)&&(i22b==6) || (i21b==6)&&(i22b==1) 
        doneP(i11b,i12b,i21b,i22b) = 1;
        continue    
    else
        P = [i11b i12b; i21b i22b]; % matrix with phases
        phases_all{end+1} = abs(M_int).*P;
        [A, ss, cycles, h] = all_topologies_analyze(single_cell, P, M_int, draw_diagram);
        state_diagrams{end+1} = A;
        steady_states{end+1} = ss;
        cycles_all{end+1} = cycles;
        
        count = count+1;
        %save_folder = fullfile('H:\My Documents\Multicellular automaton\figures\two_signals\all_topologies\temp');
        save_folder = 'H:\My Documents\Multicellular automaton\figures\two_signals\all_topologies\single_cell_by_topology';
        phase = abs(M_int).*P;
        fname_str = sprintf('state_diagram_%s_phase_%d_%d_%d_%d', label,...
            phase(1,1), phase(1,2), phase(2,1), phase(2,2));
        fname = fullfile(save_folder, fname_str);
        save_figure(h, 7, 6, fname, '.pdf', save_diagram)
            
        % update tracking variables: also consider P symmetries
        doneP(i11b,i12b,i21b,i22b) = 1;
        %gP = g(P);
        %doneP(gP(1,1),gP(1,2),gP(2,1),gP(2,2))=1;
        countP = countP + 1;
        pause(1);
        close all
    end
end
fprintf('Total # phases considered: %d \n', sum(countP))
%}
%% Functions

% symmetry operation 1<->2, finds the topology that is equivalent up to
% symmetry
function gM = g(M)
    gM = [M(2,2) M(2,1); M(1,2) M(1,1)];
end
