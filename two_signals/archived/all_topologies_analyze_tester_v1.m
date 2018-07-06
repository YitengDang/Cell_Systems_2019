%function [A, ss, cycles] = all_topologies_analyze(P)
    % Analyzes the dynamics of a system with given interaction topology and
    % phases of the interactions. 
    % Input:    P, phase matrix (2x2) with phases (1-6) for each interaction.
    %           If an interaction i<-j is absent, P(i,j) = 0.
    % Output:   A, 4x4 matrix, state diagram of the system
    %           ss, vector, all steady states
    %           cycles, cell, all cycles of the system
    
    % Initial conditions
    clear variables
    close all
    
    M_int = [1 -1; 1 1];
    phase = [1 5; 5 6];
    
    [A, ss, cycles] = all_topologies_analyze(phase, M_int);
    draw_state_diagram(A)
    
    states = {'(0,0', '(1,0)', '(0,1)', '(1,1)'};
    disp('Steady states:');
    disp(states(ss));
    
    disp('Found cycles:');
    for i=1:numel(cycles)
        disp(cycles{i});
    end
    
    %% Calculate state diagram
    %---------Modification of plot_state_diagram_multicell.m---------
    % 1: none (activation-deactivation)
    % 2: all ON(+) / OFF(-)
    % 3: A1: ON->ON (+) / ON-> OFF (-)
    % 4: A0: OFF->OFF (+) / OFF->ON (-)
    % 5: all OFF(+) / ON(-)
    % 6: A01: autonomy (+) / autonomous oscillations (-)

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
    
    %% Calculate transitions
    gij = cell(2);
    X_out = cell(2, 1);
    for i=1:2
        if all(M_int(i,:)==0)
            fprintf('No input for gene %d \n', i);
            % no input => output=initial state
            X1_in = [0 0; 1 1]; 
            X2_in = [0 1; 0 1];
            X_in = (i==1).*X1_in + (i==2).*X2_in; 
            X_out{i} = X_in;
        else
            % normal case
            for j=1:2
                if M_int(i,j)~=0
                    idx = (M_int(i,j)==1) + (M_int(i,j)==-1)*2;
                    gij{i,j} = g_map{idx, phase(i,j), j};
                else
                    gij{i,j} = ones(2); % Fixed ambiguous inputs
                end
            end
            X_out{i} = and3(gij{i,1}, gij{i,2}); % three-valued logic
        end
    end
    %% Calculate state diagram
    A = zeros(4); % graph adjacency matrix
    for i=1:2
        for j=1:2
            state_in = i + 2*(j-1); 
            X_out_this = [X_out{1}(i,j) X_out{2}(i,j)]; % tentative 
            %disp(X_out_this)
            if all(X_out_this~=2) % unambiguous out state
                state_out = X_out_this(1)+1 + 2*X_out_this(2); % (i,j) -> idx
                A(state_in, state_out) = 1;
                %disp(state_out);
            elseif sum(X_out_this==2)==1 % semi-definite
                if (X_out_this(1)==2)
                    X_out_both = [0 X_out_this(2); 1 X_out_this(2)];
                elseif (X_out_this(2)==2)
                    X_out_both = [X_out_this(1) 0; X_out_this(1) 1];
                end
                state_out = X_out_both*[1; 2]+1;
                %[X_out_both(1,1)+1 + 2*X_out_both(1,2);...
                %    X_out_both(2,1)+1 + 2*X_out_both(2,2)];
                A(state_in, state_out) = 1;
                %disp(state_out);
            elseif sum(X_out_this==2)==2 
                A(state_in, :) = 1;
            end
        end
    end
    
    %% Draw state diagram
    draw_state_diagram(A)
    
    %% Find steady states
    ss = find(diag(A)==1);
    
    %% Find cycles
    cycles = all_topologies_cycle_finder(A);
    
    disp('Found cycles:');
    for i=1:numel(cycles)
        disp(cycles{i});
    end
    
    %% Functions
    function out = and3(x,y)
        % out = ceil((x.*y)/2); % this is for 1=UNKNOWN, 2=TRUE (not used)
        out = min(x.*y, 2);
    end
    
    function draw_state_diagram(A)
        h10 = figure(1);
        hold on
        s = [0 1 0 1];
        t = [0 0 1 1];
        %A = ones(4);
        Gs = digraph(A);
        nLabels = {};
        g=plot(Gs, 'XData', s, 'YData', t, 'ArrowSize', 20, 'EdgeAlpha', 1, ...
            'LineWidth', 3, 'EdgeColor', 'k',...
            'Marker', 'o', 'MarkerSize', 100, 'NodeColor', [0.2 0.2 0.2], 'NodeLabel', nLabels);
        % Make edges dashed if state has two outgoing edges
        for i=1:4
            if sum(A(i,:))==2
                idx = find(A(i,:));
                highlight(g, i, idx, 'LineStyle', '--', 'LineWidth', 2);
            elseif sum(A(i,:))==4
                highlight(g, i, 1:4, 'LineStyle', ':', 'LineWidth', 2);
            end
        end
        text(s-0.11,t+0.019,{'(0,0)','(1,0)','(0,1)','(1,1)'}, 'Color', 'w', 'FontSize', 32)
        ax = gca;
        axis([-0.4 1.4 -0.4 1.4]);
        ax.Visible = 'off';
        h10.Color = [1 1 1];
        %set(ax, 'Units', 'Inches', 'Position', [0 0 9 8]);
        %set(h2, 'Units', 'Inches', 'Position', [1 1 9 8]);
        set(ax, 'Units', 'Inches', 'Position', [0 0 7 6]);
        set(h10, 'Units', 'Inches', 'Position', [0.2 0.2 7 6]);
    end
%end