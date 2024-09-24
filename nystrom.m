function [V,Lambda,est] = nystrom(A,s,varargin)
%NYSTROM Nystrom approximation to matrix A along with jackknife variance
%estimate or leave-one-out error estimate
%
%   Jackknife variance estimates can be applied to either the plain
%   approximation X = V*Lambda*V' or transformed versions of the
%   approximation like spectral projectors. For transformed versions, the
%   user must input a function transform such if the Nystrom replicates are
%   X = V*M*V', then the transformed replicates are V*transform(Q,D)*V'
%   where [Q,D] = eig(M)
%   
%   Arguments:
%   - A: psd matrix
%   - s: approximation rank
%   Optional arguments:
%   - transform: function handle taking [Q,D] = eig(M) and the
%   tranformed replicate is V*transform(Q,D)*V'
%   - err_est_type: either 'jack' or 'leave-one-out'). Defaults to
%   leave-one-out unless transform is set, in which case 'jack' is the only
%   available option

[transform, err_est_type, q] = read_varargin(varargin{:});

Omega = randn(size(A,1),s);
if q > 0
    Omega0 = Omega; Z = A*Omega; Omega = Z;
end
for i = 2:q
    Omega = A*Omega;
end
Y = A*Omega;
nu = eps()*norm(Y);
Y = Y + nu*Omega;
[Q,R] = qr(Y,0);
H = Omega'*Y;
C = chol((H+H')/2);
[U,Sigma,~] = svd(R/C,'econ');
Lambda = max(Sigma^2 - nu*eye(s), 0);
V = Q*U;

warning('off','MATLAB:nearlySingularMatrix');
diagHinv = diag(inv(H))';
warning('on','MATLAB:nearlySingularMatrix');

if isempty(transform)
    if strcmp(err_est_type, 'jack')
        warning('off','MATLAB:nearlySingularMatrix');
        T = ((U'*R/C)/(C')) .* diagHinv.^(-0.5);
        warning('on','MATLAB:nearlySingularMatrix');
        avg = T*T'/s;
        est = 0;
        for i = 1:s
            est = est + norm(avg - T(:,i)*T(:,i)','fro')^2;
        end
        est = sqrt(est);
    else
        if q ~= 0
            warning('off','MATLAB:nearlySingularMatrix');
            T = ((U'*R/C)/(C')) .* diagHinv.^(-0.5);
            warning('on','MATLAB:nearlySingularMatrix');
            VOm = V'*Omega0;
            est = norm(Z-V*(Lambda*VOm)+V*(T.*diag(T'*VOm).'),'fro')/sqrt(s);
        else
            est = norm(((R/C)/(C')) ./ diagHinv, 'fro')/sqrt(s);
        end
    end
else
    T = ((U'*R/C)/(C')) .* diagHinv.^(-0.5);
    replicates = zeros(s,s,s);
    lambda = diag(Lambda); 
    for i = 1:s
        [Q,D] = update_eig(lambda, T(:,i), -1);
        replicates(:,:,i) = transform(Q,D);
    end
    est = norm(replicates - mean(replicates, 3), 'fro');
end

end

