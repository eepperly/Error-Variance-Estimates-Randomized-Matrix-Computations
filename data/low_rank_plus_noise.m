function A = low_rank_plus_noise(n,R,xi)
G = randn(n,n);
A = diag([ones(1,R) zeros(1,n-R)]) + xi/n * (G*G');
end

