function lg = dos_Wang_Landau_loop(cells_in, f_in, edges, ffinal, p_flat, nbins, dist, a0)
    iteration=1;
    cells = cells_in;
    f = f_in;
    lg = zeros(1,nbins);
    while f > ffinal
        fprintf('Iteration: %d \n', iteration);
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

            %-check flatness criterion (all H(E) > p*<H(E)>)
            flat = all(H > mean(H)*p_flat); %flatness criterion
            disp(flat);
            %end
            % --end flat loop
        end
        % update f
        f = f^(1/2);
        fprintf('N=%d \n', sum(cells));
        fprintf('f=%.5f \n', f);
        iteration = iteration+1;
    end
    %f_out = f;
end