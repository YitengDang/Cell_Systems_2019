%% Obtain statistics on all possible topologies (up to equivalence)
clear variables
close all

%% Count # of topologies
n = zeros(1,5);
for k=0:4
    n(k+1) = nchoosek(4,k)*2^k;
end
p = 6.^(0:4);
Ns1 = sum(n.*p);

% account for symmetry
ns = [1 0 4 0 4]; % self-similar
ps = [1 0 21 0 610]; % #states for self-similar topologies
nr = [0 8 20 32 12]./2; % remaining
Ns2 = sum(ns.*ps+nr.*p);

% exclude trivial topologies
ns = [0 0 2 0 4];
ps = [1 0 21 0 610];
nr = [0 0 10 32 12]./2; % remaining
Ns3 = sum(ns.*ps+nr.*p);

% account for symmetries
fprintf('Total # topologies: %d \n', Ns1)
fprintf('After accounting for symmetry: %d \n', Ns2)
fprintf('After excluding trivial topologies: %d \n', Ns3)

%% Loop over topologies
done = zeros(3,3,3,3);
M = [0 1 -1];
count = 0;
countP = zeros(3^4, 1);

for k=1:3^4
    [i11, i12, i21, i22] = ind2sub([3, 3, 3, 3], k);
    if done(i11,i12,i21,i22)
        continue
    elseif (i11==1 && i12==1) || (i21==1 && i22==1) || (i12==1 && i21==1)
        % skip: no input to gene 1 || 2 || no cross-talk
        % between gene 1 and 2
        done(i11,i12,i21,i22)=1;
        gM = g([i11 i12; i21 i22]);
        done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
        %M_int = [M(i11) M(i12); M(i21) M(i22)];
        %disp(M_int);
    else 
        M_int = [M(i11) M(i12); M(i21) M(i22)];
        disp(M_int);
        
        done(i11,i12,i21,i22)=1;
        gM = g([i11 i12; i21 i22]);
        done(gM(1,1),gM(1,2),gM(2,1),gM(2,2))=1;
        
        % loop over all phases
        ni = sum(abs(M_int(:))==1);
        sz = 6*abs(M_int)+(1-abs(M_int));
        doneP = zeros(sz(1,1), sz(1,2), sz(2,1), sz(2,2));
        for k1=1:6^ni
            [i11b, i12b, i21b, i22b] = ind2sub([sz(1,1), sz(1,2), sz(2,1), sz(2,2)], k1);
            if doneP(i11b,i12b,i21b,i22b)
                continue
            else
                doneP(i11b,i12b,i21b,i22b) = 1;
                gP = g([i11b i12b; i21b i22b]);
                doneP(gP(1,1),gP(1,2),gP(2,1),gP(2,2))=1;
                countP(k) = countP(k) + 1;
            end
        end
        disp(k1);
        disp(sum(doneP(:)==1));
        %{
        if all(all(g(M_int) == M_int))
            disp(M_int);
            % special case -> to work out
        else
            ni = sum(abs(M_int(:))==1); %number of interactions
            countP = countP + 6^ni;
            for k=1:6^ni                
                sz = 6*abs(M_int)+(1-abs(M_int));
                [i11, i12, i21, i22] = ind2sub([sz(1,1), sz(1,2), sz(2,1), sz(2,2)], k);
                %countP = countP + 1;
            end
        end
        %}
        count = count+1;
    end
end
%%
countP = 0;
M_int = [1 0; 1 0];
% loop over all phases
ni = sum(abs(M_int(:))==1);
sz = 6*abs(M_int)+(1-abs(M_int));
doneP = zeros(sz(1,1), sz(1,2), sz(2,1), sz(2,2));
for k1=1:6^ni
    [i11, i12, i21, i22] = ind2sub([sz(1,1), sz(1,2), sz(2,1), sz(2,2)], k1);
    if doneP(i11,i12,i21,i22)
        continue
    else
        doneP(i11, i12, i21, i22) = 1;
        gP = g([i11 i12; i21 i22]);
        doneP(gP(1,1),gP(1,2),gP(2,1),gP(2,2))=1;
        countP = countP + 1;
    end
end

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

% symmetry operation
function gM = g(M)
    gM = [M(2,2) M(2,1); M(1,2) M(1,1)];
end
