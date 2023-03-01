load('figSX_globalSignal_data');
%%
lbl = {'saline','iGluR anatg','mAChR antag','D1/2R antag'};

[~,p] = ttest(global_da(:,1), global_da(:,2));
[~,p(2)] = ttest(global_da(:,1), global_da(:,3));
[~,p(3)] = ttest(global_da(:,1), global_da(:,4));
fig = figure; fig.Position(3) = 1000;
subplot(1,2,2); hold on
plot(global_da','.m', 'MarkerSize', 20); 
errorbar([1.2:1:4.2], nanmean(global_da), SEM(global_da,1), '.k', 'MarkerSize', 20); 
xlim([0.5 4.5]); xticks([1:5]); xticklabels(lbl);
ylabel('Global DA signal'); yticks([0:0.25:1.25]); ylim([0 1.25]);
title(sprintf('DA: sal/iglur = %1.3f | sal/machr = %1.3f | sal/d1d2 = %1.5f',p(1),p(2),p(3)))
title(sprintf('DA: sal/iglur = %1.3f | sal/machr = %1.3f | sal/d1d2 = 1.1e-04',p(1),p(2)))

[~,p] = ttest(global_ach(:,1), global_ach(:,2));
[~,p(2)] = ttest(global_ach(:,1), global_ach(:,3));
[~,p(3)] = ttest(global_ach(:,1), global_ach(:,4));
subplot(1,2,1); hold on
plot(global_ach','.g', 'MarkerSize', 20); 
errorbar([1.2:1:4.2], nanmean(global_ach), SEM(global_ach,1), '.k', 'MarkerSize', 20); 
xlim([0.5 4.5]); xticks([1:5]); xticklabels(lbl);
ylabel('Global ACh signal'); yticks([0:0.25:1.25]); ylim([0 1.25]);
title(sprintf('ACh: sal/iglur = %1.3f | sal/machr = %1.3f | sal/d1d2 = %1.3f',p(1),p(2),p(3)))
title(sprintf('ACh: sal/iglur = %1.3f | sal/machr = 6.9e-06 | sal/d1d2 = %1.3f',p(1),p(3)))