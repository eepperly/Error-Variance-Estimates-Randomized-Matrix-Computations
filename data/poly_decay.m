function A = poly_decay(n,R,p)
A = spdiags([ones(1,R) (2:(n-R+1)).^(-p)].',0,n,n);
end

