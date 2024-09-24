function [clusters, jack, coords, V, lambda] = spectral_clustering(pts,kernel,l,s,...
    varargin)

%% Set up
if isnumeric(kernel) && size(kernel,1) == size(pts,1)
    K = kernel;
else
    if isnumeric(kernel)
        kernel = @(x,y) exp(-pdist2(x,y,"euclidean").^2 / (2*kernel^2));
    end
    K = kernel(pts, pts);
end
d = sum(K, 2);
A = d.^(-0.5) .* (K .* (d').^(-0.5));

if ~isempty(varargin)
    q = varargin{1};
else
    q = 3;
end

%% Nystrom approximation
Omega = randn(size(A,1),s);
for i = 1:q
    Omega = A*Omega;
end
Y = A*Omega;
nu = eps()*norm(Y);
Y = Y + nu*Omega;
[Q,R] = qr(Y,0);
H = Omega'*Y;
C = chol((H+H')/2);
[U,Sigma,~] = svd(R/C,'econ');
lambda = max(diag(Sigma).^2 - nu, 0);
V = Q*U;

%% Spectral clustering
coords = d.^(-0.5).*V(:,1:l);
clusters = kmeans(coords, l);

%% Jackknife variance estiamte
G = V' * (V ./ d);
diagHinv = diag(inv(H))';
T = (U'*R/H) .* diagHinv.^(-0.5);
replicates = zeros(s,l,s); avg = zeros(s);
for i = 1:s
    Q = update_eig(lambda,T(:,i),-1);
    Q = Q(:,1:l);
    Q = Q / sqrt(norm(Q'*G*Q,'fro'));
    replicates(:,:,i) = Q;
    avg = avg + Q*Q'/s;
end
GavgG = G*avg*G;
jack = s * trace(GavgG*avg);
for i = 1:s
    Q = replicates(:,:,i);
    jack = jack - 2*trace(Q'*GavgG*Q) + norm(Q'*G*Q, 'fro')^2;
end
jack = sqrt(jack);

%% Fix clusters
idxes = cell(l,1);
x_coords = zeros(l,1);
for i = 1:l
    indices = find(clusters == i);
    idxes{i} = indices;
    x_coords(i) = pts(indices(1),1);
end
[~,p] = sort(x_coords);
for i = 1:l
    clusters(idxes{p(i)}) = i;
end

end