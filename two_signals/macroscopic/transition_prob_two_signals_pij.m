function W = transition_prob_two_signals_pij(M_int, Con, Coff, K, fN, gN, pij_in, I_in)
    % Calculate transition matrix
    % pij: Background p^{ij}
    if nargin<8     
        I_in = [0 0];
    end
    
    W = zeros(4, 4); %transition matrix
    % loop over p_in, p_out
    for i=1:4
        [i1, i2] = ind2sub([2 2], i);
        p_in = [i1-1 i2-1]; % order: (0,0), (1,0), (0,1), (1,1)
        %disp(p_in);
        for j=1:4
            [i1, i2] = ind2sub([2 2], j);
            p_out = [i1-1 i2-1];
            
            p_nei = [sum(pij_in(2, :)) sum(pij_in(:, 2))];
            %p_in = [0 0];
            %p_out = [0 0];
            
            % -------calculate W1 (prob. p_out^{s} = 1)--------------------
            %W1 = zeros(1, 2);

            % self-contributions to Y
            Y_self = Con.*p_in + Coff.*(1-p_in);

            % neighbour contributions to Y
            Y_nei_avg = fN.*(p_nei.*Con + (1-p_nei).*Coff + (Con-1).*(1-p_nei).*I_in); % neighbour contributions

            % total sensed concentration
            Y = Y_self + Y_nei_avg;
            Y_mat = repmat(Y, 2, 1);
            %disp(Y_mat);
            
            % std
            sigma = real(sqrt(p_nei.*(1-p_nei).*gN).*(Con-1)); 
            sigma_mat = repmat(sigma, 2, 1);

            % Note that W1 = Ponon
            %disp( M_int.*(Y_mat - K) );
            %disp(sigma_mat);
            W1 = prod(1 - abs(M_int).*normcdf(0, M_int.*(Y_mat - K), sigma_mat), 2 )';
            % disp(W1);
            % -------------------------------------------------------------

            P1 = (1-W1).*(1-p_out) + W1.*p_out;
            W(i,j) = prod(P1);
        end
    end
    
    % Check whether the probabilities add up to 1
    if ~all(round(sum(W, 2), 5)==1)
        warning('Probabilities do not add up correctly!');
        disp('Probabilities add up to:');
        disp(sum(W, 2));
    end
end