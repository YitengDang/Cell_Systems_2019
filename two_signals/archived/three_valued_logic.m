%% Implement 3-valued logic in MATALB
s = [0 1 2]; % states

and_table = zeros(3);
or_table = zeros(3);
for i=1:3
    for j=1:3
        and_table(i,j) = and(s(i), s(j));
        or_table(i,j) = or(s(i), s(j));
    end
end
disp('AND table');
disp(and_table);

disp('OR table');
disp(or_table);

function out = not(x)
    switch x
        case 0
            out = 1;
        case 1 
            out = 0;
        case 2
            out = 2;
    end
end

function out = and(x,y)
    out = min(x*y,2);
end

function out = or(x,y)
    if x==1 || y==1
        out = 1;
    elseif x==2 || y==2
        out = 2;
    else
        out = 0;
    end
end