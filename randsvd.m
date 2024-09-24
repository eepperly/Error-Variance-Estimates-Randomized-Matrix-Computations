function [U,S,V,est] = randsvd(A,s,varargin)
%RANDSVD Randomized SVD approximation for matrix A along with jackknife
%variance estimate or leave-one-out error estimate
%
%   Jackknife variance estimates can be applied to either the plain
%   approximation X = U*S*V' or transformed versions of the approximation
%   like left or right singular projectors. For transformed versions, the
%   user must input a function transform such if the RSVD replicates are
%   X = U*M*V', then the transformed replicates are U*transform(u,s,v)*V'
%   where [u,s,v] = svd(M)
%   
%   Arguments:
%   - A: matrix
%   - s: approximation rank
%   Optional arguments:
%   - transform: function handle taking [u,s,v] = svd(M) and the
%   tranformed replicate is U*transform(u,s,v)*V'
%   - err_est_type: either 'jack' or 'leave-one-out'. Defaults to
%   leave-one-out unless transform is set, in which case 'jack' is the only
%   available option

[transform, err_est_type, q] = read_varargin(varargin{:});

Omega = randn(size(A,2),s);
Z = A*Omega; Y = Z;
for iter = 1:q
    Y = A*(A'*Y);
end
[Q,R] = qr(Y,0);
C = Q'*A;

[W,S,V] = svd(C,'econ');
U = Q * W;

if isempty(transform)
    if strcmp(err_est_type, 'jack')
        T = cnormc(inv(R'));
        avg = T*T'*W*S/s;
        est = 0;
        for i = 1:s
            est = est + norm(avg - T(:,i)*((T(:,i)'*W)*S),'fro')^2;
        end
        est = sqrt(est);
    else
        if q ~= 0
            warning('off','MATLAB:nearlySingularMatrix');
            T = cnormc(inv(R'));
            warning('on','MATLAB:nearlySingularMatrix');
            QZ = Q'*Z;
            est = norm(Z-Q*QZ+Q*T.*(diag(T'*QZ).'),'fro')/sqrt(s);
        else
            est = norm(1./vecnorm(inv(R')))/sqrt(s);
        end
    end
else
    warning('off','MATLAB:nearlySingularMatrix');
    T = cnormc(inv(R'));
    warning('on','MATLAB:nearlySingularMatrix');
    WT = W'*T;
    replicates = zeros(s,s,s);
    for i = 1:s
        [UU,SS,VV] = svd(S - WT(:,i)*(WT(:,i)'*S));
        replicates(:,:,i) = transform(UU,SS,VV);
    end
    est = norm(replicates - mean(replicates, 3), 'fro');
end
end

