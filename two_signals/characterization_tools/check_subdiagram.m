function match = check_subdiagram(A, subdiagram)
    % checks whether a subdiagram is present in a given state diagram
    % A: input state diagram
    % subdiagram: n x 2 matrix with transitions (i,j) which need to be
    % present
    
    % subdiagrams for travelling wave propagation
    %{
    subdiagrams = {};
    % anti-clockwise loop
    subdiagrams{1} = [1 2; 2 4; 3 1; 4 3]; % sets of (i,j) indices of transitions that need to be present
    % clockwise loop
    subdiagrams{2} = [2 1; 4 2; 1 3; 3 4]; % sets of (i,j) indices of transitions that need to be present
    %}
    
    if sum(diag(A))==0
        % no stable state
        match = 0;
        return
    elseif all(A(:)==1)
        % disregard unconstrained diagram
        match = 0;
        return
    else
        match = 1;
        % check all conditions
        for i=1:size(subdiagram, 1)
            if ~A(subdiagram(i, 1), subdiagram(i, 2))
                match = 0;
                return
            end
        end
    end
end
