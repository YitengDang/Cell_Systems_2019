clear variables 
close all

    %% Testing code: periodicity_test_cycle
    % Works
    clc
    test_cycles = {};
    found_cycle = [];
    test_cycles{end+1} = [1 6 8 5 9 0 1 2 3];
    test_cycles{end+1} = [randperm(10,5) randperm(10,5)];
    test_cycles{end+1} = [randperm(10,5) randperm(10,5)];
    test_cycles{end+1} = [randperm(10,5) randperm(10,5)];
    
    for i=1:numel(test_cycles)
        [periodic, found_cycle] = periodicity_test_cycle(test_cycles{i});
        disp(test_cycles{i});
        fprintf('Periodic? %d, found_cycle =  \n', periodic);
        disp(found_cycle)
    end
    
    %% Testing code: trim_cycles(cycles)
    % Works
    clc
    [found_cycles, remaining_cycles] = trim_cycles(test_cycles);
    
    disp('Found cycles:');
    for i=1:numel(found_cycles)
        disp(found_cycles{i});
    end
    disp('Remaining cycles:');
    for i=1:numel(remaining_cycles)
        disp(remaining_cycles{i});
    end
    

    %% Testing code: shift_cycle
    % Passed
    clc
    
    % shift cycle
    cycle = [1 6 8 5 9 0 1 2 3];
    cycle_trimmed = trim_cycles({cycle});
    cycle_trimmed = cycle_trimmed{1};
    
    shifted_cycle = cycle_trimmed;
    disp('Original cycle:');
    disp(cycle_trimmed);
    disp('Shifted cycles:');
    for i=1:numel(cycle_trimmed)
        shifted_cycle = shift_cycle(shifted_cycle);
        disp(shifted_cycle);
    end
    
    %% Testing: compare_cycles(cycle1, cycle2)
    % 
    clc
    cycle1 = [1 2 3 4 1];
    cycle2 = [2 3 4 1 2];
    same = compare_cycles(cycle1, cycle2)
    
    cycle1 = [1 2 3 4 1];
    cycle2 = [4 3 2 1 4];
    same = compare_cycles(cycle1, cycle2)
    
    cycle1 = [1 3 4 5 3 1];
    cycle2 = [3 4 5 3 1 3];
    same = compare_cycles(cycle1, cycle2)
    % Wrong
    %% Testing code: update_cycle
    % Works (but double check)
    % Next: implement while loop instead of fixed # trials
    clc;
    
    A = [1 1 0 0; 0 0 0 1; 1 1 0 0; 0 0 1 1];
    %A = randi(2, 4)-1;
    Gs = digraph(A);
    g=plot(Gs);
    
    start_node = 1;
    cycles = {start_node};
    net_cycles = {}; % "net" cycles
    
    for trial=1:15
        [updated_cycles] = update_cycle(A, cycles);
        [found_cycles, remaining_cycles] = trim_cycles(updated_cycles);
        for i=1:numel(found_cycles)
            net_cycles{end+1} = found_cycles{i};
        end
        cycles = remaining_cycles;
    end
    %% remove equivalent cycles
    net_cycles2 = {net_cycles{1}};
    for i1=1:numel(net_cycles)
        for i2=1:numel(net_cycles2)
            cycle1_trimmed = trim_cycles({net_cycles{i1}}); cycle1_trimmed = cycle1_trimmed{1};
            cycle2_trimmed = trim_cycles({net_cycles{i2}}); cycle2_trimmed = cycle2_trimmed{1};
            same = compare_cycles(cycle1_trimmed, cycle2_trimmed);
            if same
                break % do nothing, go to next cycle
            end
        end
        if ~same
            net_cycles2{end+1} = net_cycles{i1};
        end
    end
    
    disp('Found cycles:');
    for i=1:numel(net_cycles)
        disp(net_cycles{i});
    end
    
    disp('Found cycles (after removing equivalent ones):');
    for i=1:numel(net_cycles2)
        disp(net_cycles2{i});
    end
    %% Functions
    function [found_cycles, remaining_cycles] = trim_cycles(cycles)
        remaining_cycles = {};
        found_cycles = {};
        for i=1:numel(cycles)
            cycle = cycles{i};
            
            % test for periodicity
            [periodic, found_cycle] = periodicity_test_cycle(cycle);
            % if not periodic, keep cycle
            if periodic
                % register trimmed cycle
                found_cycles{end+1} = found_cycle;
            else
                % keep original cyle
                remaining_cycles{end+1} = cycle;
            end
        end
        
        % return number array instead of cell array if there is only 1
        % element
        %if length(found_cycles)==1
        %    found_cycles = found_cycles{1};
        %end
        
    end
    
    function [periodic, found_cycle] = periodicity_test_cycle(cycle)
        if numel(unique(cycle))==cycle
            periodic = 0;
            found_cycle = [];
            return
        else
            for t1=1:length(cycle)
                for t2=t1+1:length(cycle)
                    if (cycle(t1)==cycle(t2))
                        periodic = 1;
                        found_cycle = cycle(t1:t2);
                        return
                    end
                    
                end
            end
            periodic = 0;
            found_cycle = [];
        end
    end
    
    function same = compare_cycles(cycle1, cycle2)
        % compares whether two cycles are equivalent
        % input: cycles in form of number arrays
        same = 0;
        if length(cycle1)==length(cycle2)
            shifted_cycle = cycle2;
            for i=1:numel(shifted_cycle)
                shifted_cycle = shift_cycle(shifted_cycle);
                if all(shifted_cycle==cycle1)
                    same = 1;
                    return
                end
            end
            
        end
    end
    
    function shifted_cycle = shift_cycle(cycle)
        % shifts a cycle by 1 place
        if cycle(1)~=cycle(end)
            warning('Shift_cycle: Not a cycle!');
        else
            shifted_cycle = cycle(1:end-1);
            shifted_cycle = [shifted_cycle(end) shifted_cycle];
        end
    end
    function updated_cycles = update_cycle(A, cycles)
        updated_cycles = {};
        for i1=1:numel(cycles)
            % find children
            cycle = cycles{i1};
            last_node = cycle(end);
            children_nodes = children(A, last_node);
            
            % add children nodes to cycles
            if numel(children_nodes)==0 
                % if node has no children, return same cycle
                updated_cycles{end+1} = cycle;
            else
                % update to new cycles
                %new_cycles = cell(numel(children_nodes), 1);
                for i2=1:numel(children_nodes)
                    new_cycle = cycle;
                    new_cycle(end+1) = children_nodes(i2);
                    %new_cycles{i2} = new_cycle;
                    updated_cycles{end+1} = new_cycle;
                end
            end
            
            % add new cycles to input cycles
            % updated_cycles{i1} = new_cycles;
        end
    end
    
    function out = children(A, node)
        idx = find(A(node, :)==1);
        out = idx(idx~=node);
    end