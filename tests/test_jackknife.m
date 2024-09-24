addpath('..')
addpath('../data')
method = 'randsvd';

n = 1000;
R = 5;
q = 0;
As = { exp_decay(n,R,0.1), exp_decay(n,R,0.05),...
       poly_decay(n,R,2), poly_decay(n,R,1),...
       low_rank_plus_noise(n,R,1e-4), low_rank_plus_noise(n,R,1e-2) };
names = {'expfast','expslow','polyfast','polyslow','nlrfast','nlrslow'};
s_vals = 10:10:100;
trials = 1000;
close all

if strcmp(method, 'nystrom')
    % load('A_kernel.mat')
    % As{end+1} = A_kernel;
    % names{end+1} = 'kernel';
    As(1) = [];
    names(1) = [];
    transform = @(VV,DD) VV(:,1:5)*VV(:,1:5)';
elseif strcmp(method, 'rsvd') || strcmp(method, 'randsvd')
    % load('A_vel.mat')
    % As{end+1} = A_vel;
    % names{end+1} = 'vel';
    transform = @(UU,SS,VV) VV(:,1:5)*VV(:,1:5)';
else
    error("Method '%s' not recognized", method)
end

for run = 1:length(As)
    A = As{run};
    name = names{run};
    if strcmp(method, 'rsvd') || strcmp(method, 'randsvd')
        [U,S,V] = svd(full(A), 'econ');
        target = transform(U,S,V);
    elseif strcmp(method, 'nystrom')
        [V,D] = eig(full(A));
        V = V(:,end:-1:1); D = D(end:-1:1,end:-1:1);
        target = transform(V,D);
    else
        error("Method '%s' not recognized", method)
    end

    stds = zeros(length(s_vals),1);
    errs = zeros(length(s_vals),1);
    jacks = zeros(length(s_vals),1);
    jack_stds = zeros(length(s_vals),1);

    for s_idx = 1:length(s_vals)
        s = s_vals(s_idx);
        fprintf('s=%d\t',s)
        jack_trials = zeros(trials,1);
        avg = zeros(size(target));
        rng(42)
        for trial = 1:trials
            if mod(trial, ceil(trials/10)) == 0
                fprintf('.')
            end
            if strcmp(method, 'rsvd') || strcmp(method, 'randsvd')
                [U,S,V,jack_trials(trial)] = randsvd(A,s,q,transform);
                result = transform(U,S,V);
            elseif strcmp(method, 'nystrom')
                [V,D,jack_trials(trial)] = nystrom(A,s,q,transform);
                result = transform(V,D);
            else
                error("Method '%s' not recognized", method)
            end
            avg = avg + result / trials;
            errs(s_idx) = errs(s_idx) + norm(result - target,'fro')/trials;
        end
        jacks(s_idx) = mean(jack_trials);
        jack_stds(s_idx) = std(jack_trials);

        rng(42)
        for trial = 1:trials
            if mod(trial, ceil(trials/10)) == 0
                fprintf('.')
            end
            if strcmp(method, 'rsvd') || strcmp(method, 'randsvd')
                [U,S,V,jack_trials(trial)] = randsvd(A,s,q,transform);
                result = transform(U,S,V);
            elseif strcmp(method, 'nystrom')
                [V,D,jack_trials(trial)] = nystrom(A,s,q,transform);
                result = transform(V,D);
            else
                error("Method '%s' not recognized", method)
            end
            stds(s_idx) = stds(s_idx) + norm(result - avg,'fro')^2/trials;
        end
        stds(s_idx) = sqrt(stds(s_idx));
        fprintf('\n')
    end

    figure
    errorbar(s_vals,jacks,jack_stds,'--','LineWidth',2,'Color',"#D95319")
    hold on
    semilogy(s_vals, errs,'Color',"#0072BD")
    semilogy(s_vals, stds,':o','Color',"#EDB120",...
        'MarkerFaceColor',"#EDB120")
    set(gca, 'YScale', 'log')
    if run == 1
        legend({'$\mathrm{Jack}(\mbox{\boldmath $X$})$',...
            '$\mathrm{Std}(\mbox{\boldmath $X$})$',...
            '$\|\mbox{\boldmath $X$}-\mbox{\boldmath $\Pi$}\|_{\rm F}$'},...
            "Location","best")
    end
    xlabel('Approximation rank $s$')
    ylabel('Error metric')
    drawnow
    saveas(gcf, sprintf('../figs/jack_%s_%s.fig', method, name))
    saveas(gcf, sprintf('../figs/jack_%s_%s.png', method, name))
end