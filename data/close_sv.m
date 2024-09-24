function A = close_sv(n,r,eps,small_sv)
A = diag([ones(1,r) 1-eps small_sv*ones(1,n-1-r)]);
end

