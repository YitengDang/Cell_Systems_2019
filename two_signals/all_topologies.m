%% Obtain statistics on all possible topologies (up to equivalence)
clear variables
close all

%% Count # of topologies
n = zeros(1,5);
n_phases = 6; % number of phases to distinguish
for k=0:4
    n(k+1) = nchoosek(4,k)*2^k;
end
p = n_phases.^(0:4);
Ns1 = sum(n.*p);

% account for symmetry
ns = [1 0 4 0 4]; % self-similar
ps = [1 0 21 0 666]; % #states for self-similar topologies
nr = [0 8 20 32 12]./2; % remaining
Ns2a = sum(nr.*p); % states without self-similarity
Ns2 = sum(ns.*ps+nr.*p);

% exclude trivial topologies
ns = [0 0 2 0 4];
nr = [0 0 10 32 12]./2; % remaining
Ns3 = sum(ns.*ps+nr.*p);

% account for symmetries
fprintf('Total # phases: %d \n', Ns1)
fprintf('Accounting for symmetry: %d \n', Ns2)
fprintf('States without self-similarity: %d \n', Ns2a)
fprintf('Excluding trivial topologies: %d \n', Ns3)

%% Loop over topologies and phases
% N.B. some of the phases are not possible, e.g. when autonomy and
% activation-deactivation are both present for the same gene

done = zeros(3,3,3,3); % keeps track of which topologies have been found already (up to symmetry)
M = [0 1 -1]; % index to interaction
count = 0;  %count topologies
countP = zeros(3^4, 1); % count phases
M_int_all = cell(3^4, 1);

% main output 
phases_all = {}; % store all phases, stored as matrix with phase of each interaction (1-6), 0 if no interaction (M_int(i,j)==0)
state_diagrams = {}; % cell(Ns2, 1); % store all state diagrams (graph transition matrices)
steady_states = {}; %cell(Ns2, 1); % store all steady states
cycles_all = {}; %cell(Ns2, 1); % store all loop structures

for k=1:3^4
    [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], k);
    % matrix associated with indices
    M_int = [M(i11) M(i12); M(i21) M(i22)];
    M_int_all{k} = M_int;
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
                    P = [i11b i12b; i21b i22b]; % matrix with phases
                    phases_all{end+1} = abs(M_int).*P;
                    
                    % update tracking variables: also consider P symmetries
                    doneP(i11b,i12b,i21b,i22b) = 1;
                    gP = g(P);
                    doneP(gP(1,1),gP(1,2),gP(2,1),gP(2,2))=1;
                    countP(k) = countP(k) + 1;
                end
            end
        else
            for k1=1:n_phases^ni
                [i11b, i12b, i21b, i22b] = ind2sub([sz(1,1), sz(1,2), sz(2,1), sz(2,2)], k1);
                if doneP(i11b,i12b,i21b,i22b)
                    continue
                else
                    % update tracking variables: no need to consider P symmetries
                    doneP(i11b,i12b,i21b,i22b) = 1;
                    countP(k) = countP(k) + 1;
                end
            end
        end
        %}
        count = count+1;
    end
end

fprintf('Total # phases considered: %d \n', sum(countP))


% save data

%% For given topology (M_int), loop over all phases
countP = 0;
M_int = [1 0; 1 1];
% loop over all phases
ni = sum(abs(M_int(:))==1);
sz = 6*abs(M_int)+(1-abs(M_int));
doneP = zeros(sz(1,1), sz(1,2), sz(2,1), sz(2,2));

phases_all = {}; % store all phases, stored as matrix with phase of each interaction (1-6), 0 if no interaction (M_int(i,j)==0)
state_diagrams = {}; % cell(Ns2, 1); % store all state diagrams (graph transition matrices)
steady_states = {}; %cell(Ns2, 1); % store all steady states
cycles_all = {}; %cell(Ns2, 1); % store all loop structures

for k1=1:n_phases^ni
    [i11b, i12b, i21b, i22b] = ind2sub([sz(1,1), sz(1,2), sz(2,1), sz(2,2)], k1);
    if doneP(i11b,i12b,i21b,i22b)
        continue
    else
        % 
        P = [i11b i12b; i21b i22b]; % matrix with phases
        phases_all{end+1} = abs(M_int).*P;
        [A, ss, cycles] = all_topologies_analyze(P);
        state_diagrams{end+1} = A;
        steady_states{end+1} = ss;
        cycles_all{end+1} = cycles;
        
        % update tracking variables: also consider P symmetries
        doneP(i11b,i12b,i21b,i22b) = 1;
        gP = g(P);
        doneP(gP(1,1),gP(1,2),gP(2,1),gP(2,2))=1;
        countP = countP + 1;
    end
end
fprintf('Total # phases considered: %d \n', sum(countP))

%% Transition table
g_map = cell(2, 6, 2);
% 0=OFF, 1:ON, 2:UNKNOWN
% activation 
g_map{1,1,1} = 2*ones(2);
g_map{1,1,2} = 2*ones(2);
g_map{1,2,1} = ones(2);
g_map{1,2,2} = ones(2);
g_map{1,3,1} = [2 2; 1 1];
g_map{1,3,2} = [2 1; 2 1];
g_map{1,4,1} = [0 0; 2 2];
g_map{1,4,2} = [0 2; 0 2];
g_map{1,5,1} = zeros(2);
g_map{1,5,2} = zeros(2);
g_map{1,6,1} = [0 0; 1 1];
g_map{1,6,2} = [0 1; 0 1];
% repression 
%(note: this is precisely NOT g_map{1,:,:} in the three-val
% boolean algebra with NOT 2 = 2)
g_map{2,1,1} = 2*ones(2);
g_map{2,1,2} = 2*ones(2);
g_map{2,2,1} = zeros(2);
g_map{2,2,2} = zeros(2);
g_map{2,3,1} = [2 2; 0 0];
g_map{2,3,2} = [2 0; 2 0];
g_map{2,4,1} = [1 1; 2 2];
g_map{2,4,2} = [1 2; 1 2];
g_map{2,5,1} = ones(2);
g_map{2,5,2} = ones(2);
g_map{2,6,1} = [1 1; 0 0];
g_map{2,6,2} = [1 0; 1 0];


%% Functions

% symmetry operation 1<->2, finds the topology that is equivalent up to
% symmetry
function gM = g(M)
    gM = [M(2,2) M(2,1); M(1,2) M(1,1)];
end
