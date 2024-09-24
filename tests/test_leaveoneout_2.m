addpath('..')
addpath('../data')
method = 'nystrom';

n = 1000;
R = 5;
As = { exp_decay(n,R,0.1), exp_decay(n,R,0.05),...
       poly_decay(n,R,2), poly_decay(n,R,1),...
       low_rank_plus_noise(n,R,1e-4), low_rank_plus_noise(n,R,1e-2) };
names = {'expfast','expslow','polyfast','polyslow','nlrfast','nlrslow'};
s_vals = 10:10:100;
trials = 2;
close all

if strcmp(method, 'nystrom')
    % load('A_kernel.mat')
    % As{end+1} = A_kernel;
    % names{end+1} = 'kernel';
elseif strcmp(method, 'rsvd') || strcmp(method, 'randsvd')
    % load('A_vel.mat')
    % As{end+1} = A_vel;
    % names{end+1} = 'vel';
else
    error("Method '%s' not recognized", method)
end

for run = 1:length(As)
    A = As{run};
    name = names{run};
    norm_A = norm(A,'fro');

    errs = zeros(length(s_vals),1);
    leaveoneouts = zeros(length(s_vals),1);
    leaveoneout_stds = zeros(length(s_vals),1);

    for s_idx = 1:length(s_vals)
        s = s_vals(s_idx);
        leaveoneout_trials = zeros(trials,1);
        avg = zeros(size(A));
        rng(42)
        for trial = 1:trials
            if strcmp(method, 'rsvd') || strcmp(method, 'randsvd')
                [U,S,V,leaveoneout_trials(trial)] = randsvd(A,s);
                result = U*S*V';
            elseif strcmp(method, 'nystrom')
                [V,D,leaveoneout_trials(trial)] = nystrom(A,s);
                result = V*D*V';
            else
                error("Method '%s' not recognized", method)
            end
            avg = avg + result / trials;
            errs(s_idx) = errs(s_idx) + norm(result - A,'fro')...
                /trials/norm_A;
        end
        leaveoneouts(s_idx) = mean(leaveoneout_trials) / norm_A;
        leaveoneout_stds(s_idx) = std(leaveoneout_trials) / norm_A;
    end

    figure
    errorbar(s_vals,leaveoneouts,leaveoneout_stds,'--','LineWidth',2,'Color',"#D95319")
    hold on
    semilogy(s_vals, errs,'Color',"#0072BD")
    set(gca, 'YScale', 'log')
    if run == 1
        legend({'$\mathrm{Err}(\mbox{\boldmath $X$}) / \|\mbox{\boldmath $A$}\|_{\rm F}$',...
            '$\|\mbox{\boldmath $X$}-\mbox{\boldmath $A$}\|_{\rm F} / \|\mbox{\boldmath $A$}\|_{\rm F}$'},...
            "Location","best")
    end
    xlabel('Approximation rank $s$')
    ylabel('Error metric')
    drawnow
    saveas(gcf, sprintf('../figs/leaveoneout_%s_%s.fig', method, name))
    saveas(gcf, sprintf('../figs/leaveoneout_%s_%s.png', method, name))
end