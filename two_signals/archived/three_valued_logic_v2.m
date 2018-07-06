%% Implement 3-valued logic in MATALB
clear all
close all
%%
s = [0 1 2]; % states
s1 = repmat(s, 3, 1);
s2 = repmat(s', 1, 3);
% 0: False
% 1: True
% 2: UNKNOWN

and_table = zeros(3);
or_table = zeros(3);
for i=1:3
    for j=1:3
        and_table(i,j) = and3(s(i), s(j));
        or_table(i,j) = or(s(i), s(j));
    end
end

disp('NOT test');
test = randi(3, 5, 5)-1;
disp(test);
disp(not3(test));

disp('AND table');
%disp(and_table);
disp(and3(s1,s2));

disp('OR table');
%disp(or_table);
disp(or3(s1,s2));

%%
function out = and3(x,y)
    out = min(x.*y, 2);
end

function out = or3(x,y)
    out = zeros(size(x));
    
    idx2 =  (x==2)|(y==2);
    out(idx2) = 2;
    
    idx1 = (x==1)|(y==1);
    out(idx1) = 1;
end

function out = not3(x)
    idx0 = (x==0);
    idx1 = (x==1);
    
    out = x;
    out(idx0) = 1;
    out(idx1) = 0;
end