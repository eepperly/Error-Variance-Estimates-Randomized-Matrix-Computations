addpath('..')
load('../data/clustering.mat')
num_trials = 1000;
s_vals = 5:5:150;
correct = zeros(size(s_vals));
variances = zeros(size(s_vals));
jacks = zeros(size(s_vals));
jack_stds = zeros(size(s_vals));
errs = zeros(size(s_vals));
bandwidth = 0.05;

kernel = @(x,y) exp(-pdist2(x,y,"euclidean").^2 / (2*bandwidth^2));
K = kernel(X,X);
d = sum(K, 2);
A = d.^(-0.5) .* (K .* (d').^(-0.5));

for i = 1:length(s_vals)
    s = s_vals(i);
    fprintf('s=%d:\t', s)
    jack_trials = zeros(num_trials, 1);
    avg = zeros(size(X,1), size(X,1));
    rng(42)
    for trial = 1:num_trials
        if mod(trial, ceil(num_trials/10)) == 0
            fprintf('.');
        end
        [clusters, jack_trials(trial), coords, V, lambda] ...
            = spectral_clustering(X,K,4,s);
        Y = coords * coords'; Y = Y / norm(Y,'fro');
        correct(i) = correct(i) + (norm(clusters - correct_clustering)...
            ==0)/num_trials;
        avg = avg + Y/num_trials;
        errs(i) = errs(i) + norm(A - V*diag(lambda)*V', 'fro')/norm(A,'fro')/num_trials;
    end
    jacks(i) = mean(jack_trials);
    jack_stds(i) = std(jack_trials);
    rng(42)
    for trial = 1:num_trials
        if mod(trial, ceil(num_trials/10)) == 0
            fprintf('.');
        end
        [~, ~, coords] = spectral_clustering(X,K,4,s);
        Y = coords * coords'; Y = Y / norm(Y,'fro');
        variances(i) = variances(i) + norm(Y - avg,'fro')^2 / num_trials;
    end
    variances(i) = sqrt(variances(i));
    fprintf('\n%e\t%e\t%e\t%f\n', variances(i), jacks(i), errs(i), correct(i));
    save('spectral_clustering_test.mat')
end

%% Plot
close all

yyaxis left 
errorbar(s_vals, jacks, jack_stds, '--', 'LineWidth', 2)
hold on
semilogy(s_vals, variances, '-', 'LineWidth', 2)
semilogy(s_vals, errs, '*-.', 'LineWidth', 2, 'Color',"#00008B",'MarkerFaceColor',"#00008B")
set(gca,'YScale','log')
axis([50 150 6e-4 3e0])

xlabel('Approximation rank $s$')
ylabel('Quality metric')

yyaxis right
plot(s_vals, correct, ':')
ylabel('Probability of correct recovery')
axis([50 150 0 1.1])

legend({'$\mathrm{Jack}(\mbox{\boldmath $X$})$',...
    '$\mathrm{Std}(\mbox{\boldmath $X$})$',...
    '$\mathrm{Err}(\mbox{\boldmath $\hat{A}$},\mbox{\boldmath $A$})/\|\mbox{\boldmath $A$}\|_{\rm F}$',...
    'Success chance'}, "Location","best")

saveas(gcf, '../figs/spectral_clustering.png')
saveas(gcf, '../figs/spectral_clustering.fig')

%% Clustering plots

clusters = spectral_clustering(X,0.05,4,50);
figure 
scatter(X(:,1),X(:,2),10,clusters,'filled')
axis equal
axis([0 4 0 1])
axis off
colormap([100,143,255;120,94,240;220,38,127;254,97,0]/255)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

saveas(gcf, '../figs/clustering_incorrect.png')
saveas(gcf, '../figs/clustering_incorrect.fig')

clusters = spectral_clustering(X,0.05,4,150);
figure 
scatter(X(:,1),X(:,2),10,clusters,'filled')
axis equal
axis([0 4 0 1])
axis off
colormap([100,143,255;120,94,240;220,38,127;254,97,0]/255)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

saveas(gcf, '../figs/clustering_correct.png')
saveas(gcf, '../figs/clustering_correct.fig')