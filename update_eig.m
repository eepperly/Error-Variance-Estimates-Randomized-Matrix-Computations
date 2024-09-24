function [Q,D] = update_eig(d,v,c)
%% Negate, sort, and normalize
negative = c < 0;
if negative; d = -d; c = -c; end
[d,p] = sort(d,'ascend');
v = v(p);
c = c * norm(v)^2;
v = v / norm(v);

%% Compute update
try 
    [Q,D] = update_eig_mex(d,v,c);
catch ME
    if strcmp(ME.identifier, "MATLAB:UndefinedFunction")
        error(['Need to build MEX function '...
            '"mex -v update_eig_mex.c -lmwlapack"'])
    else
        rethrow(ME)
    end
end

%% Undo negation and sorting
Q(p,:) = Q;
if negative; D = -D; end
D = diag(D);
end