function A = exp_decay(n,R,q)
A = spdiags([ones(1,R) 10.^(-q*(1:(n-R)))].',0,n,n);
end

