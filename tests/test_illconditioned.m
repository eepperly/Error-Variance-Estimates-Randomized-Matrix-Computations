%% Setup
clear 
close all
addpath('..')
load('../data/A_clustering.mat')
A = A_clustering;
clear A_clustering
[V_full,D_full] = eig(A);

s_vals = 50:5:150;
q = 3;
trials = 1000;

ks = [3 4];
errs_full = cell(length(ks),1);
stds_full = cell(length(ks),1);
jacks_full = cell(length(ks),1);
jack_stds_full = cell(length(ks),1);

%% Run
for k_idx = 1:length(ks)
    k = ks(k_idx);
    fprintf('k=%d\n',k)
    errs= zeros(size(s_vals));
    stds = zeros(size(s_vals));
    jacks = zeros(size(s_vals));
    jack_stds = zeros(size(s_vals));

    transform = @(VV,DD) VV(:,1:k) * VV(:,1:k)';
    target = transform(V_full, D_full);
    norm_target = norm(target,'fro');

    for s_idx = 1:length(s_vals)
        s = s_vals(s_idx);
        jack_trials = zeros(trials,1);
        avg = zeros(size(target));
        rng(42)
        fprintf('s=%d\t',s)
        for trial = 1:trials
            if mod(trial, ceil(trials/10)) == 0
                fprintf('.');
            end
            [V,D,jack_trials(trial)] = nystrom(A,s,q,transform);
            result = transform(V,D);
            avg = avg + result / trials;
            errs(s_idx) = errs(s_idx) + norm(result - target,'fro')...
                / trials / norm_target;
        end
        jacks(s_idx) = mean(jack_trials) / norm_target;
        jack_stds(s_idx) = std(jack_trials) / norm_target;

        rng(42)
        for trial = 1:trials
            if mod(trial, ceil(trials/10)) == 0
                fprintf('.');
            end
            [V,D,jack_trials(trial)] = nystrom(A,s,q,transform);
            result = transform(V,D);
            stds(s_idx) = stds(s_idx) + norm(result - avg,'fro')^2/trials;
        end
        stds(s_idx) = sqrt(stds(s_idx)) / norm_target;

        % save('illconditioned_test.mat')
        fprintf('\n')
    end

    errs_full{k_idx} = errs;
    stds_full{k_idx} = stds;
    jacks_full{k_idx} = jacks;
    jack_stds_full{k_idx} = jack_stds;
end

%% Plots
figure
errorbar(s_vals,jacks_full{1},jack_stds_full{1},'--','LineWidth',2,'Color',"#0072BD")
hold on
semilogy(s_vals, stds_full{1},'Color',"#4DBEEE")
errorbar(s_vals,jacks_full{2},jack_stds_full{2},':','LineWidth',2,'Color',"#A2142F")
semilogy(s_vals, stds_full{2},'-.','Color',"#D95319")
set(gca,'YScale','log')
legend({'$\mathrm{Jack}(\mbox{\boldmath $X$})$ ($n_{\rm dim}=3$)',...
    '$\mathrm{Std}(\mbox{\boldmath $X$})$ ($n_{\rm dim}=3$)',...
    '$\mathrm{Jack}(\mbox{\boldmath $X$})$ ($n_{\rm dim}=4$)',...
    '$\mathrm{Std}(\mbox{\boldmath $X$})$ ($n_{\rm dim}=4$)'},...
    'Location','southwest')
axis([50 150 4e-4 4e0])
xlabel('Approximation rank $s$')
ylabel('Standard deviation')
saveas(gcf, '../figs/illconditioned.fig')
saveas(gcf, '../figs/illconditioned.png')