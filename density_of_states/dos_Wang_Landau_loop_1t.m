function lg = dos_Wang_Landau_loop_1t(cells_in, f_in, edges, ffinal, p_flat, nbins, dist, a0)
    % Implements 1/t algorithm from "Wang-Landau algorithm: A theoretical analysis of the saturation
    % of the error", https://arxiv.org/pdf/cond-mat/0702414
    iteration=1;
    time = 0; % keep track of simulation steps in units of 10^5
    cells = cells_in;
    N = sum(cells);
    f = f_in;
    lg = zeros(1,nbins);
    tc_reached = 0; % tC: time when F(t) ~ 1/t should start
    while f > ffinal
        fprintf('N=%d \n', N);
        fprintf('Iteration: %d \n', iteration);
        fprintf('f=%.5f, f_final = %.5f \n', f, ffinal);
        fprintf('time = %d * 10^5 \n', time);
        H = zeros(1,nbins); % initialize/reset histogram  
        % --inner loop
        %for i=1:1000
        flat=0;
        while ~flat
            %disp('check');
            for i=1:10000 % Do 10000 MC steps before checking
                % flip two random spins: 1 ON and 1 OFF
                cellsp = cells;
                ONcells = find(cells == 1); x1=randi(length(ONcells)); cellsp(ONcells(x1))=0; % turn ON cell OFF
                OFFcells = find(cells == 0); x2=randi(length(OFFcells)); cellsp(OFFcells(x2))=1; % turn OFF cell ON
                I = moranI(cells, a0*dist);
                Ip = moranI(cellsp, a0*dist); 
                bin = dos_Wang_Landau_binning(I, edges); % bin of original system
                binp = dos_Wang_Landau_binning(Ip, edges); % bin of flipped system
                %disp(g(E)); disp(g(Ep));
                %disp(lg(E)); disp(lg(Ep));
                prob = min(exp(lg(bin)-lg(binp)), 1); %acceptance probability min(g()/g(E2), 1)
                %p = min(g(E)/g(Ep), 1); %acceptance probability min(g()/g(E2), 1)
                if rand() < prob
                    cells = cellsp;
                    %I = Ip;
                    bin = binp;
                end
            %disp(lattice);
            %disp(strcat('new h=', num2str(h)));

            %-update ge and energies
            %g(E) = f_in*g(E);
            lg(bin) = lg(bin) + log(f);
            H(bin) = H(bin) + 1;
            %disp(g);
            %disp(H);
            end
            % update time
            time = time + 1;
            
            %-check flatness criterion (all H(E) > p*<H(E)>)
            flat = all(H > mean(H)*p_flat); %flatness criterion
            disp(flat);
            %end
            % --end flat loop
        end
        % update f
        % f = f^(1/2); % original rule
        if ~tc_reached 
            tc_reached = log(f) < 1/(time*10^5);
            f = f^(1/2);
        else
            f = exp(1/(time*10^5));
        end
        iteration = iteration+1;
        
    end
    %f_out = f;
end