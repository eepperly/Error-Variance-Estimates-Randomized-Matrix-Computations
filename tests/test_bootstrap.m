rng(42389324)
d = 1000;
s = 100;
B = 1000;
trials = 1000;

true_svals = [linspace(1,0.26,75) 0.25*(1:(d-75)).^(-2)]';
A = spdiags(true_svals,0,d,d);
Om = randn(d,s);

Y = A*Om;
Q = orth(Y);
sval = norm(Q'*A);

svals = zeros(s,1);
fprintf('jack')
for i = 1:s
    if mod(i,s/20) == 0
        fprintf('.');
    end
    Q = orth(Y(:,[1:(i-1) (i+1):s]));
    svals(i) = norm(Q'*A);
end
fprintf('\n')
jack = norm(svals - mean(svals))

svals = zeros(B,1);
fprintf('boot')
for i = 1:B
    if mod(i,B/20) == 0
        fprintf('.');
    end
    Q = orth(Y(:,unique(randsample(s,s,true))));
    svals(i) = norm(Q'*A);
end
fprintf('\n')
boot = norm(svals - mean(svals))/sqrt(B-1)

svals = zeros(trials,1);
fprintf('true')
for i = 1:trials
    if mod(i,trials/20) == 0
        fprintf('.');
    end
    Y = A*randn(d,s);
    Q = orth(Y);
    svals(i) = norm(Q'*A);
end
fprintf('\n')
true_std = norm(svals - mean(svals))/sqrt(trials-1)