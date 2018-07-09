%% Obtain statistics on all possible topologies (up to equivalence)
clear variables
close all
clc
%% Count # of topologies
n_phases_all = [6 2];

% Single cell
disp('Single cell:');
n_phases = n_phases_all(2); % number of phases for each interaction
Ts = [1 0 3 0 10];
Ta = n_phases.^(0:4);
ns = [1 0 4 0 4]; % self-similar
na = [0 8 20 32 12]./2; % remaining

N_phases = sum(Ts.*ns+Ta.*na); % total # phases
fprintf('Total # single-cell phases: %d \n', N_phases)

% multiple cells

disp('Multiple cells:');
n_phases = n_phases_all(1);
n = zeros(1,5);
for k=0:4
    n(k+1) = nchoosek(4,k)*2^k;
end
na = n_phases.^(0:4);
Ns1 = sum(n.*na);

% account for symmetry
Ts = [1 0 4 0 4]; % self-similar
Ta = [0 8 20 32 12]./2; % remaining
ns = [1 0 21 0 666]; % #states for self-similar topologies
%Ns2a = sum(Ta.*na); % states without self-similarity
Ns2 = sum(Ts.*ns+Ta.*na);

% exclude trivial topologies
Ts = [0 0 2 0 4];
Ta = [0 0 10 32 12]./2; % remaining
Ns3 = sum(Ts.*ns+Ta.*na);

% account for symmetries
fprintf('Total # multi-cell phases: %d \n', Ns1)
fprintf('Accounting for symmetry: %d \n', Ns2)
%fprintf('States without self-similarity: %d \n', Ns2a)
fprintf('Excluding trivial topologies: %d \n', Ns3)

% Convention:
% phase 1: (Con - K) > 0
% phase 2: (Con - K) < 0

%% Loop over topologies and phases
% N.B. some of the phases are not possible, e.g. when autonomy and
% activation-deactivation are both present for the same gene

% Settings
single_cell = 0;
draw_diagram = 0; % draw state diagram?
n_phases = n_phases_all(single_cell+1);

% tracking variables
done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
M = [0 1 -1]; % index to interaction
count = 0; %count topologies
countP = zeros(3^4, 1); % count phases

% main output 
M_int_all = {}; %cell(3^4, 1);
phases_all = {}; % store all phases, stored as matrix with phase of each interaction (1-6), 0 if no interaction (M_int(i,j)==0)
state_diagrams = {}; % cell(Ns2, 1); % store all state diagrams (graph transition matrices)
steady_states = {}; %cell(Ns2, 1); % store all steady states
cycles_all = {}; %cell(Ns2, 1); % store all loop structures

for k=1:3^4
    [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], k);
    % matrix associated with indices
    M_int = [M(i11) M(i12); M(i21) M(i22)];
    
    if done(i11,i12,i21,i22)
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
        
        done(i11,i12,i21,i22)=1;
        gM = g([i11 i12; i21 i22]);
        done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
        
        % loop over all phases
        ni = sum(abs(M_int(:))==1);
        sz = n_phases*abs(M_int)+(1-abs(M_int)); % size matrix, for correct index conversion
        doneP = zeros(sz(1,1), sz(1,2), sz(2,1), sz(2,2));
        % distinguish between self-similar topologies and others
        if all(all(g(M_int) == M_int))
            % special case: topology symmetry 1<->2
            %disp(M_int);
            for k1=1:n_phases^ni
                [i11b, i12b, i21b, i22b] = ind2sub([sz(1,1), sz(1,2), sz(2,1), sz(2,2)], k1);
                if doneP(i11b,i12b,i21b,i22b)
                    continue
                else
                    % 
                    disp(sum(countP));
                    P = [i11b i12b; i21b i22b]; % matrix with phases
                    phases_all{end+1} = abs(M_int).*P; % store phase matrix 
                    M_int_all{end+1} = M_int; %store interaction matrix
                    
                    [A, ss, cycles] = all_topologies_analyze(single_cell, P, M_int, draw_diagram);
                    state_diagrams{end+1} = A;
                    steady_states{end+1} = ss;
                    cycles_all{end+1} = cycles;
                    
                    % update tracking variables: also consider P symmetries
                    doneP(i11b,i12b,i21b,i22b) = 1;
                    gP = g(P);
                    doneP(gP(1,1),gP(1,2),gP(2,1),gP(2,2))=1;
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
                else
                    % 
                    disp(sum(countP));
                    P = [i11b i12b; i21b i22b]; % matrix with phases
                    phases_all{end+1} = abs(M_int).*P; % store phase matrix
                    M_int_all{end+1} = M_int; %store interaction matrix
                    
                    [A, ss, cycles] = all_topologies_analyze(single_cell, P, M_int, draw_diagram);
                    state_diagrams{end+1} = A;
                    steady_states{end+1} = ss;
                    cycles_all{end+1} = cycles;
                    
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
        count = count+1;
    end
end

fprintf('Total # phases considered: %d \n', sum(countP))
%%
% save data
%save_path = 'H:\My Documents\Multicellular automaton\data\two_signals\all_topologies';
save_path = 'D:\Multicellularity\data\two_signals\all_topologies';
labels = {'multi_cell', 'single_cell'};
qsave = 1;
if qsave
    label = labels{single_cell+1};
    fname_str = sprintf('all_topologies_data_%s', label);
    fname = fullfile(save_path, fname_str);
    save(fname, 'M_int_all', 'phases_all', 'state_diagrams', 'steady_states', 'cycles_all', 'single_cell');
end

%% For given topology (M_int), loop over all phases
countP = 0;
M_int = [1 0; 0 1];

% Settings
single_cell = 1;
draw_diagram = 1; % draw state diagram?

% loop over all phases
n_phases = n_phases_all(single_cell+1);
ni = sum(abs(M_int(:))==1);
sz = n_phases*abs(M_int)+(1-abs(M_int));
doneP = zeros(sz(1,1), sz(1,2), sz(2,1), sz(2,2));

phases_all = {}; % store all phases, stored as matrix with phase of each interaction (1-6), 0 if no interaction (M_int(i,j)==0)
state_diagrams = {}; % cell(Ns2, 1); % store all state diagrams (graph transition matrices)
steady_states = {}; %cell(Ns2, 1); % store all steady states
cycles_all = {}; %cell(Ns2, 1); % store all loop structures

for k1=1:n_phases^ni
    disp(k1)
    [i11b, i12b, i21b, i22b] = ind2sub([sz(1,1), sz(1,2), sz(2,1), sz(2,2)], k1);
    if doneP(i11b,i12b,i21b,i22b)
        continue
    else
        P = [i11b i12b; i21b i22b]; % matrix with phases
        phases_all{end+1} = abs(M_int).*P;
        [A, ss, cycles] = all_topologies_analyze(single_cell, P, M_int, draw_diagram);
        state_diagrams{end+1} = A;
        steady_states{end+1} = ss;
        cycles_all{end+1} = cycles;
        
        % update tracking variables: also consider P symmetries
        doneP(i11b,i12b,i21b,i22b) = 1;
        gP = g(P);
        doneP(gP(1,1),gP(1,2),gP(2,1),gP(2,2))=1;
        countP = countP + 1;
        pause(1);
        close all
    end
end
fprintf('Total # phases considered: %d \n', sum(countP))

%% Functions

% symmetry operation 1<->2, finds the topology that is equivalent up to
% symmetry
function gM = g(M)
    gM = [M(2,2) M(2,1); M(1,2) M(1,1)];
end
