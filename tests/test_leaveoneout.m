addpath('..')
load('../data/A_kernel.mat')
A = A_kernel;
num_trials = 1000;
s_vals = 5:5:200;
ests = zeros(length(s_vals),1);
est_stds = zeros(length(s_vals),1);
gh_ests = zeros(length(s_vals),1);
gh_est_stds = zeros(length(s_vals),1);
errs = zeros(length(s_vals),1);
est_errs = zeros(length(s_vals),1);
gh_est_errs = zeros(length(s_vals),1);
rng(42);

for i = 1:length(s_vals)
    s = s_vals(i);
    fprintf('s=%d\t', s)
    avg = zeros(size(A));
    est_copies = zeros(num_trials,1);
    normest_copies = zeros(num_trials,1);
    normA = norm(A,'fro');
    for trial = 1:num_trials
        if mod(trial, ceil(num_trials/10)) == 0
            fprintf('.');
        end
        [V, D, est_copies(trial)] = nystrom(A,s);
        Ahat = V*D*V';
        normest_copies(trial) = norm((A - Ahat)*randn(size(A,1),10),'fro')/sqrt(10);
        avg = avg + V*D*V'/num_trials;
        myerr = norm(A - Ahat,'fro');
        errs(i) = errs(i) + myerr/num_trials/normA;
        est_errs(i) = est_errs(i) + abs(myerr - est_copies(trial))/myerr/num_trials;
        gh_est_errs(i) = gh_est_errs(i) + abs(myerr - normest_copies(trial))/myerr/num_trials;
    end
    ests(i) = mean(est_copies) / normA;
    est_stds(i) = std(est_copies) / normA;
    gh_ests(i) = mean(normest_copies) / normA;
    gh_est_stds(i) = std(normest_copies) / normA;
    fprintf('\n')
end

%% Plots
close all
figure(1)
semilogy(s_vals, errs); hold on
errorbar(s_vals, ests, est_stds,'--','LineWidth',2)
axis([0 150 1e-3 3e-2])
xlabel('Approximation rank $s$')
ylabel('Average relative error')
legend({'Error $\|\mbox{\boldmath $A$}-\mbox{\boldmath $X$}\|_{\rm F}/\|\mbox{\boldmath $A$}\|_{\rm F}$','Estimate $\mathrm{Err}(\mbox{\boldmath $X$})/\|\mbox{\boldmath $A$}\|_{\rm F}$'})
saveas(gcf,'../figs/leaveoneout.fig')
saveas(gcf,'../figs/leaveoneout.png')

figure(2)
semilogy(s_vals, gh_est_errs, '-.', 'Color', "#EDB120",'LineWidth',2);
hold on
semilogy(s_vals, est_errs, '--', 'Color', "#D95319",'LineWidth',2);
axis([0 150 -Inf Inf])
xlabel('Approximation rank $s$')
ylabel('Relative error of norm estimate')
legend({'Girard--Hutchinson','Leave-one-out'})
saveas(gcf,'../figs/leaveoneout_compare.fig')
saveas(gcf,'../figs/leaveoneout_compare.png')
